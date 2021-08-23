include("../src/MimiNICE.jl")
using Mimi
using CSV, CSVFiles, DataFrames
function scenario_compare()
    nice = Main.MimiNICE.create_nice()
    prtps = [0, 0.01, 0.02]
    etas = [1.0, 2.0]
    whose_moneys = ["independent", "average", "poorest"]

    nice = Main.MimiNICE.create_nice()
    run(nice)
    variant_nice = Main.MimiNICE.create_nice()
    run(variant_nice)
    base_emissions = nice[:emissions, :EIND]
    scenarios_file = DataFrame(load(joinpath(@__DIR__, "..", "data", "SSP_IAM_V2_201811.csv")))
    scenarios_of_interest = [
        ("IMAGE", "SSP1-19"), ("GCAM4", "SSP4-34"), ("MESSAGE-GLOBIOM", "SSP2-45"), 
        ("AIM/CGE", "SSP3-60"), ("REMIND-MAGPIE", "SSP5-Baseline")
    ]
    results_table = DataFrame(Model=String[], Scenario=String[], whose=String[], prtp=Float64, eta=Float64, scc=Float64)
    # TODO: correct for inflation
    for item in scenarios_of_interest
        model = item[1]
        scenario = item[2]
        variant_emissions = scenarios_file[
            (scenarios_file.MODEL .== model) .& (scenarios_file.SCENARIO .== scenario) .& 
            (scenarios_file.REGION .== "World") .& (scenarios_file.VARIABLE .== "Emissions|CO2|Fossil Fuels and Industry"), :]
        variant_treemissions = scenarios_file[
            (scenarios_file.MODEL .== model) .& (scenarios_file.SCENARIO .== scenario) .& 
            (scenarios_file.REGION .== "World") .& (scenarios_file.VARIABLE .== "Emissions|CO2|Land Use"), :]
        emissions_ratios = DataFrame(year=Int64[], value=Float64[])
        treemissions = DataFrame(year=Int64[], value=Float64[])
        years_to_mod = (2015:10:2095)
        # We start at the second year index, 2015
        year_ind = 2
        for year in years_to_mod
            # Convert interpolated average scenario emissions in MtCO2/year into factor for GtC/year
            emissions_ratio = (variant_emissions[1, string(year-5)] + variant_emissions[1, string(year+5)]) / 2 / sum(
                base_emissions[year_ind, :]) * 12/44 / 1000
            push!(emissions_ratios, [year, 1-emissions_ratio])
            treemission = (variant_treemissions[1, string(year-5)] + variant_treemissions[1, string(year+5)]) / 2 * 12/44 / 1000
            push!(treemissions, [year, treemission])
            year_ind += 1
        end
        # Do the remaining years to mimic 2100
        emissions_ratio = variant_emissions[1, string(2100)]/ sum(base_emissions[year_ind, :]) * 12/44 / 1000
        push!(emissions_ratios, [2105, 1-emissions_ratio])
        treemission = (variant_treemissions[1, string(2100)]) * 12/44 / 1000
        push!(treemissions, [2105, treemission])
        # The code cannot deal with cases more extreme than the baseline, so we remove emissions values below 0. We also limit 
        # Negative emissions after 2100 to only tree-based emissions. 
        emissions_ratios.value[emissions_ratios.value .< 0 ] .= 0
        new_miu = vcat(vcat(zeros(1, 12), repeat(emissions_ratios.value, 1, 12)), min(emissions_ratios.value[end], 1) * ones(49, 12))
        update_param!(variant_nice, :MIU, new_miu)
        # Our new construction requires the original starting value, followed by the repeated last value
        etree_new = append!(append!([nice[:emissions, :etree][1]], treemissions.value), repeat([treemissions.value[end]], 49))
        update_param!(variant_nice, :etree, etree_new)

        for prtp in prtps
            for whose_money in whose_moneys
                for eta in etas
                    scc = Main.MimiNICE.compute_scc(variant_nice, year=2025, last_year=2205, eta=eta, prtp=prtp, whose_money=whose_money, interpolate=true)
                    newrow = DataFrame(Model=model, Scenario=scenario, whose=whose_money, prtp=prtp, eta=eta, scc=scc)
                    results_table = vcat(results_table, newrow)
                end
            end
        end
    end
    return results_table
end

# Convert results to 2020 dollars
dollar_val_2020 = 1.35

results_interp = scenario_compare()
results_interp.scc = results_interp.scc * dollar_val_2020
CSV.write("./output/interpolated_scc_scenarios_v5_2025.csv", results_interp)