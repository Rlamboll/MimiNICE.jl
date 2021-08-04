include("../src/MimiNICE.jl")
using OptiMimi
using Mimi

nice = Main.MimiNICE.create_nice()
prtp = 0.015
etas = [2.0]
whose_moneys = ["independent", "average", "poorest"]

opt_nice = Main.MimiNICE.create_nice()
function objective(model::Model)
    model[:nice_welfare, :welfare]
end


# TODO: correct for inflation
for whose_money in whose_moneys
    for eta in etas
        """
        scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2595, eta=eta, prtp=prtp, whose_money=whose_money)
        print("\n eta = ", eta, " whose = ", whose_money, " scc = ", scc, "\n")
        """
        # Calculate for the optimised trajectory
        optprob = problem(opt_nice, [:nice_welfare], [:welfare], [-1e6], [1e6], objective)
        (maxx, maxf) = solution(optprob, () -> [0. for i in 1:5])
        print("\n" maxx)
        print("\n" maxf)
        scc = Main.MimiNICE.compute_scc(opt_nice, year=2015, last_year=2595, eta=eta, prtp=prtp, whose_money=whose_money)
        print("\n optimised eta = ", eta, " whose = ", whose_money, " scc = ", scc, "\n")

    end
end


