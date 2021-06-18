include("../src/MimiNICE.jl")

nice = Main.MimiNICE.create_nice()
prtp = 0.015
etas = [2.0]
whose_moneys = ["independent", "average", "poorest"]
# TODO: correct for inflation
for whose_money in whose_moneys
    for eta in etas
        scc = Main.MimiNICE.compute_scc(nice, year=2015, last_year=2595, eta=eta, prtp=prtp, whose_money=whose_money)
        print("\n eta = ", eta, " whose = ", whose_money, " scc = ", scc, "\n")
    end
end