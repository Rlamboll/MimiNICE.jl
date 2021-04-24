print("test calc scc")

include("../src/MimiNICE.jl")

# Try with Main? 
nice = Main.MimiNICE.create_nice()
#marginal_nice = Main.MimiNICE.get_marginal_model(nice, year=2015)
#run(marginal_nice)
#print(marginal_nice[:nice_neteconomy, :l])
scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2095, eta=1.0)
print("\n scc calc complete, SCC: \n")
print(scc)