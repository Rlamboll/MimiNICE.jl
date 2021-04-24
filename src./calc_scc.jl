"""
compute_scc(m::Model=nothing, year::Int = nothing, last_year::Int = 2595), prtp::Float64 = 0.03)
Computes the social cost of CO2 for an emissions pulse in `year` for the provided MimiRICE2010 model. 
Constant discounting is used from the specified pure rate of time preference `prtp`. 
`whose_money` determines whether we calculate the SCC relative to the poorest group at time `year` (use the option 
`poorest`), the world average (`average`) or simply summing over each region (`independent`). 
"""
function compute_scc(
    m::Model = nothing; year::Union{Int, Nothing} = nothing, last_year::Int = m.mi.last, 
    prtp::Float64 = 0.015, eta::Float64 = 1.5, whose_money::String = "poorest"
)
    !(whose_money in ["poorest", "average", "independent"]) ? error(
        "whose_money must be one of `poorest`, `average` or `independent`."
    ) : nothing
    m === nothing ? error("Must specify a model. Try `compute_scc(create_nice(), year=2020)`.") : nothing    
    year === nothing ? error("Must specify an emission year. Try `compute_scc(m, year=2020)`.") : nothing
    !(last_year in model_years) ? error(
        "Invalid value of $last_year for last_year. last_year must be within the model's time index $model_years."
        ) : nothing
    !(year in model_years[1]:10:last_year) ? error(
        "Cannot compute the scc for year $year, year must be within the model's time index $(model_years[1]):10:$last_year."
        ) : nothing

    mm = get_marginal_model(m; year = year)

    return _compute_scc(mm, year=year, last_year=last_year, prtp=prtp, eta=eta, whose_money=whose_money)
end

"""
compute_scc_mm(m::Model=nothing; year::Int = nothing, last_year::Int = 2595, prtp::Float64 = 0.03)
Returns a NamedTuple (scc=scc, mm=mm) of the social cost of carbon and the MarginalModel used to compute it.
Computes the social cost of CO2 for an emissions pulse in `year` for the provided MimiRICE2010 marginal model. 
Constant discounting is used from the specified pure rate of time preference `prtp`.
`whose_money` determines whether we calculate the SCC relative to the poorest group at time `year` (use the option 
`poorest`), the world average (`average`) or simply summing over each region (`independent`). 
"""
function compute_scc_mm(
    m::Model=nothing; year::Union{Int, Nothing} = nothing, last_year::Int = m.mi.last, prtp::Float64 = 0.015, 
    eta::Float64 = 1.5, whose_money::String = "poorest"
)
    !(whose_money in ["poorest", "average", "independent"]) ? error(
        "whose_money must be one of `poorest`, `average` or `independent`."
    ) : nothing
    m === nothing ? error(
        "Must specify a model. Try `compute_scc_mm(create_nice(), year=2020)`."
        ) : nothing
    year === nothing ? error("Must specify an emission year. Try `compute_scc_mm(m, year=2015)`.") : nothing
    !(last_year in model_years) ? error(
        "Invlaid value of $last_year for last_year. last_year must be within the model's time index $model_years."
        ) : nothing
    !(year in model_years[1]:10:last_year) ? error(
        "Cannot compute the scc for year $year, year must be within the model's time index $(model_years[1]):10:$last_year."
        ) : nothing

    mm = get_marginal_model(m; year = year)
    scc = _compute_scc(mm; year=year, last_year=last_year, prtp=prtp, eta=eta, whose_money=whose_money)
    
    return (scc = scc, mm = mm)
end

# helper function for computing SCC from a MarginalModel
function _compute_scc(mm::MarginalModel; year::Int, last_year::Int, prtp::Float64, eta::Float64, whose_money::String)
    # Will run through the timestep of the specified last_year
    ntimesteps = findfirst(isequal(last_year), model_years)     
    run(mm, ntimesteps=ntimesteps)

    # convert from trillion $/ton C to $/ton CO2;
    marginal_damages = (mm.base[:nice_neteconomy, :quintile_c_post] 
        - mm.modified[:nice_neteconomy, :quintile_c_post]) * 10.0^12 * 12/44 
    # We determine the weights using the pre-damages consumption in $
    cpc = 1000 .* mm.base[:nice_neteconomy, :quintile_c_pre]

    year_index = findfirst(isequal(year), model_years)
    df = zeros(size(cpc))
    if whose_money == "poorest"
        null_cpc = minimum(cpc[year_index, :, 1], dims=1)
    end
    if whose_money == "average"
        null_cpc = sum(sum(cpc[year_index, :, :] .* mm.base[:nice_neteconomy, :l][year_index, :, :], dims=2), dims=1) /
         sum(sum(mm.base[:nice_neteconomy, :l][year_index, :, :], dims=2), dims=1)
    end
    if whose_money in ["poorest", "average"]
        for (i,t) in enumerate(model_years) 
            if year<=t<=last_year
                df[i, :, :] = (null_cpc./cpc[i, :, :]).^eta * 1/(1+prtp)^(t-year)
            end
        end
    end
    if whose_money == "independent"
        for (i,t) in enumerate(model_years) 
            if year<=t<=last_year
                df[i, :, :] = (cpc[year, :, :]./cpc[i, :, :])^eta * 1/(1+prtp)^(t-year)
            end
        end
    end
    # currently implemented as a 10year step function; so each timestep of discounted marginal damages is multiplied by 10
    # TODO: check and interpolate
    scc = sum(df .* marginal_damages * 10)
    return scc
end

"""
get_marginal_model(m::Model = get_model(); year::Int = nothing)
Creates a Mimi MarginalModel where the provided m is the base model, and the marginal model has additional emissions of CO2 in year `year`.
If no Model m is provided, the default model from MimiDICE2010.get_model() is used as the base model.
"""
function get_marginal_model(m::Model=nothing; year::Union{Int, Nothing} = nothing)
    m === nothing ? error(
        "Must specify a model. Try `get_marginal_model(create_nice(), year=2020)`."
        ) : nothing        
    year === nothing ? error("Must specify an emission year. Try `get_marginal_model(m, year=2015)`.") : nothing
    !(year in model_years) ? error("Cannot add marginal emissions in $year, year must be within the model's time index $(model_years[1]):10:$last_year.") : nothing

    mm = create_marginal_model(m, 1e10) # Pulse has a value of 1GtC per year for ten years
    add_marginal_emissions!(mm.modified, year)

    return mm
end

"""
Adds a marginal emission component to year m which adds 1Gt of additional C emissions per year for ten years starting in the specified `year`.
"""
function add_marginal_emissions!(m::Model, year::Int) 
    add_comp!(m, Mimi.adder, :marginalemission, before=:co2cycle)

    time = Mimi.dimension(m, :time)
    addem = zeros(length(time))
    addem[time[year]] = 1.0     # 1 GtC per year for ten years

    set_param!(m, :marginalemission, :add, addem)
    connect_param!(m, :marginalemission, :input, :emissions, :E)
    connect_param!(m, :co2cycle, :E, :marginalemission, :output)
end
