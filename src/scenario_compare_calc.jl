include("../src/MimiNICE.jl")
using OptiMimi
using Mimi
using CSVFiles, DataFrames

nice = Main.MimiNICE.create_nice()
prtp = 0.015
etas = [1.0]
whose_moneys = ["independent", "average", "poorest"]

nice = Main.MimiNICE.create_nice()
run(nice)
variant_nice = Main.MimiNICE.create_nice()
run(variant_nice)
base_emissions = nice[:emissions, :EIND]
scenarios_file = DataFrame(load(joinpath(@__DIR__, "data", "SSP_IAM_V2_201811.csv")))
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
        push!(emissions_ratios, [year, emissions_ratio - 1])
        treemission = (variant_treemissions[1, string(year-5)] + variant_treemissions[1, string(year+5)]) / 2 * 12/44 / 1000
        push!(treemissions, [year, treemission])
        year_ind += 1
    end
    # The code cannot deal with cases more extreme than the baseline, so we remove emissions values below 0
    emissions_ratios.value[emissions_ratios.value .< 0 ] .= 0
    new_miu = vcat(vcat(zeros(1, 12), repeat(emissions_ratios.value, 1, 12)), emissions_ratios.value[end] * ones(50, 12))
    update_param!(variant_nice, :MIU, new_miu)
    # Our new construction requires the original starting value, followed by the repeated last value
    etree_new = append!(append!([nice[:emissions, :etree][1]], treemissions.value), repeat([treemissions.value[end]], 50))
    update_param!(variant_nice, :etree, etree_new)

    for whose_money in whose_moneys
        for eta in etas
            """
            scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2595, eta=eta, prtp=prtp, whose_money=whose_money)
            print("\n eta = ", eta, " whose = ", whose_money, " scc = ", scc, "\n")
            """
            # Calculate for the optimised trajectory
        
            print("\n", maxx)
            print("\n", maxf)
            scc = Main.MimiNICE.compute_scc(variant_nice, year=2015, last_year=2595, eta=eta, prtp=prtp, whose_money=whose_money)
            print("\n optimised eta = ", eta, " whose = ", whose_money, " scc = ", scc, "\n")

        end
    end
end

