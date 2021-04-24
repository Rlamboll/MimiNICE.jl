using Mimi

include("./MimiNICE.jl")

nice = MimiNICE.create_nice()

println("created nice")
println("calculating scc")
scc = MimiNICE.compute_scc(nice, year=2015, last_year=2095)
println("calculated scc")
print(scc)
#explore(nice)
