using Mimi
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

# Test that the value returns 0 if all emissions are 0
nice = Main.MimiNICE.create_nice()
nice_2 = Main.MimiNICE.create_nice()
run(nice)
run(nice_2)
scc_orig = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="poorest")
scc_orig2105 = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2105, eta=eta, prtp=prtp, whose_money="poorest")
new_miu = zeros(60, 12)
new_miu[4:end, :] .= 0.5
update_param!(nice_2, :MIU, new_miu)
run(nice_2)
scc_upd = Main.MimiNICE.compute_scc(nice_2, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="poorest")
scc_upd2105 = Main.MimiNICE.compute_scc(nice_2, year=2015, last_year=2105, eta=eta, prtp=prtp, whose_money="poorest")
@assert(scc_orig == scc_upd)
@assert(scc_orig2105 > scc_upd2105)
@assert( nice_2[:damages, :TATM][1:3] == nice[:damages, :TATM][1:3])
@assert( nice_2[:damages, :TATM][4:end] < nice[:damages, :TATM][4:end])
etree_new = nice[:emissions, :etree]
etree_new[3:end] = etree_new[3:end] /2
update_param!(nice_2, :etree, etree_new)
run(nice_2)
scc_updtree = Main.MimiNICE.compute_scc(nice_2, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="poorest")
scc_upd2105tree = Main.MimiNICE.compute_scc(nice_2, year=2015, last_year=2105, eta=eta, prtp=prtp, whose_money="poorest")
@assert(nice_2[:damages, :TATM][1:2] == nice[:damages, :TATM][1:2])
@assert(nice_2[:damages, :TATM][3] < nice[:damages, :TATM][3])

mm = Main.MimiNICE.get_marginal_model(nice_2, year=2015)
run(mm)
# Number of years * difference * weighting, converted to units of tCO2, per unit emitted
cost = (
    10 * (mm.base[:nice_neteconomy, :quintile_c_post][3, :, :] - mm.modified[:nice_neteconomy, :quintile_c_post][3, :, :]) .* 
    repeat(mm.base[:nice_neteconomy, :l][3, :], outer=[1, 5]) * 1000000 / 5 
)

poorweight = (
    minimum(mm.base[:nice_neteconomy, :quintile_c_post][2, :, 1]) ./ 
    mm.base[:nice_neteconomy, :quintile_c_post][3, :, :]
) .^ eta
computed_val_poor = sum(poorweight .*  cost) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^10
@assert(scc_updtree == computed_val_poor)

print("\n Test passed")