include("../src/MimiNICE.jl")

# Try with Main? 
nice = Main.MimiNICE.create_nice()
prtp = 0.015
eta = 2.0
scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2025, eta=eta, prtp=prtp, whose_money="independent")
print("\n scc calc complete, SCC: \n")
print(scc)
mm = Main.MimiNICE.get_marginal_model(nice, year=2015)
run(mm)
# Number of years * difference * weighting, converted to units of tCO2, per unit emitted
computed_val = 10 * sum( (
    (mm.base[:nice_neteconomy, :quintile_c_post][3, :, :] - mm.modified[:nice_neteconomy, :quintile_c_post][3, :, :]) .* 
    repeat(mm.base[:nice_neteconomy, :l][3, :], outer=[1, 5]) * 1000000 / 5
) .* (
    mm.base[:nice_neteconomy, :quintile_c_post][2, :, :] ./ mm.base[:nice_neteconomy, :quintile_c_post][3, :, :]
).^ eta ) * 12 / 44 / mm.delta * 1000 / (1 + prtp).^10 
print("\n Computed value: \n")
print(computed_val)

@assert(abs((scc-computed_val)/computed_val) < 1e-12)