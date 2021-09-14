include("../src/MimiNICE.jl")

# First test that the replicator works
array = [1 1.1 1.2; 2 2.1 2.2; 3 3.1 3.2]
reparray = Main.MimiNICE.interp_t(array, 2)
@assert(all((reparray - [1 1.1 1.2; 1.5 1.6 1.7; 2 2.1 2.2; 2.5 2.6 2.7; 3 3.1 3.2]) .^2 .< 1e-12))

# Then test the different types of social cost calculator work
nice = Main.MimiNICE.create_nice()
prtp = 0.015
eta = 2.0
scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="independent")
sccpoor = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="poorest")
mm = Main.MimiNICE.get_marginal_model(nice, year=2015)
run(mm)
# Number of years * difference * weighting, converted to units of tCO2, per unit emitted
cost = (
    10 * (mm.base[:nice_neteconomy, :quintile_c_post][3, :, :] - mm.modified[:nice_neteconomy, :quintile_c_post][3, :, :]) .* 
    repeat(mm.base[:nice_neteconomy, :l][3, :], outer=[1, 5]) * 1000000 / 5 
)

weight = (
    mm.base[:nice_neteconomy, :quintile_c_post][2, :, :] ./ mm.base[:nice_neteconomy, :quintile_c_post][3, :, :]
).^ eta
computed_val = sum(weight .*  cost) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^10 
# Also test for the poorest weighting
poorweight = (
    minimum(mm.base[:nice_neteconomy, :quintile_c_post][2, :, 1]) ./ 
    mm.base[:nice_neteconomy, :quintile_c_post][3, :, :]
) .^ eta
computed_val_poor = sum(poorweight .*  cost) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^10 
@assert(abs((scc-computed_val)/computed_val) < 1e-12)
@assert(abs((sccpoor-computed_val_poor)/computed_val_poor) < 1e-12)

# Compare the effect of Interpolation
sccpoor_int = Main.MimiNICE.compute_scc(
    nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="poorest", interpolate=true
)
@assert(sccpoor_int < 2* sccpoor)
@assert(sccpoor_int > sccpoor)

# Compare the effect of average weighting vs poorest weighting
sccav = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="average")
conversion_factor = (minimum(mm.base[:nice_neteconomy, :quintile_c_post][2, :, :]) /
    sum(sum(mm.base[:nice_neteconomy, :quintile_c_post][2, :, :].*mm.base[:nice_neteconomy, :l][2, :, :]./
    sum(mm.base[:nice_neteconomy, :l][2, :, :])))).^eta
@assert(abs(sccpoor/sccav-conversion_factor) < 1e-15)

# Then calculate the 2035 contribution 
scc2035 = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2035, eta=eta, prtp=prtp, whose_money="independent")
poorscc2035 = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2035, eta=eta, prtp=prtp, whose_money="poorest")
cost = (
    10 * (mm.base[:nice_neteconomy, :quintile_c_post][4, :, :] - mm.modified[:nice_neteconomy, :quintile_c_post][4, :, :]) .* 
    repeat(mm.base[:nice_neteconomy, :l][4, :], outer=[1, 5]) * 1000000 / 5 
)

weight = (
    mm.base[:nice_neteconomy, :quintile_c_post][2, :, :] ./ mm.base[:nice_neteconomy, :quintile_c_post][4, :, :]
).^ eta
computed_val2035 = sum(weight .*  cost) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^20 + computed_val
@assert(abs((scc2035-computed_val2035)/computed_val2035) < 1e-12)

poorweight = (
    minimum(mm.base[:nice_neteconomy, :quintile_c_post][2, :, 1]) ./ 
    mm.base[:nice_neteconomy, :quintile_c_post][4, :, :]
) .^ eta
poorcomputed_val2035 = sum(poorweight .*  cost) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^20 + computed_val_poor
@assert(abs((poorscc2035-poorcomputed_val2035)/poorcomputed_val2035) < 1e-12)

print("\n Test passed")