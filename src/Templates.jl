include("Species.jl")

function true_network()
    net = BioNetwork()

    add_species!(net, Gene, "gCbf1")
    add_species!(net, Gene, "gGal4")
    add_species!(net, Gene, "gSwi5")
    add_species!(net, Gene, "gGal80")
    add_species!(net, Gene, "gAsh1")

    add_species!(net, mRNA, "mCbf1",  transcriptor=net["gCbf1"])
    add_species!(net, mRNA, "mGal4",  transcriptor=net["gGal4"])
    add_species!(net, mRNA, "mSwi5",  transcriptor=net["gSwi5"])
    add_species!(net, mRNA, "mGal80", transcriptor=net["gGal80"])
    add_species!(net, mRNA, "mAsh1",  transcriptor=net["gAsh1"])

    add_species!(net, Protein, "Swi5", translator=net["mSwi5"])
    add_species!(net, Protein, "Cbf1", translator=net["mCbf1"])
    add_species!(net, Protein, "Gal4", translator=net["mGal4"])
    add_species!(net, Protein, "Gal80", translator=net["mGal80"])
    add_species!(net, Protein, "Ash1", translator=net["mAsh1"])

    add_species!(net, Complex, "gp(gGal80-Swi5)",   species1=net["gGal80"], species2=net["Swi5"]),
    add_species!(net, Complex, "gp(gASh1-Swi5)",    species1=net["gASh1"],  species2=net["Swi5"]),
    add_species!(net, Complex, "gp(gCbf1-Swi5)",    species1=net["gCbf1"],  species2=net["Swi5"]),
    add_species!(net, Complex, "gp(gCbf1-Ash1)",    species1=net["gCbf1"],  species2=net["Ash1"]),
    add_species!(net, Complex, "gp(gGal4-Cbf1)",    species1=net["gGal4"],  species2=net["Cbf1"]),
    add_species!(net, Complex, "gp(gSwi5-Gal4)",    species1=net["gSwi5"],  species2=net["Gal4"]),
    add_species!(net, Complex, "c(Gal80-Gal4)",     species1=net["Gal80"],  species2=net["Gal4"])

    return net
    # species::Vector{Species} = []
    # append!(species, [
    #     Gene("gCbf1", 1),  # 1
    #     Gene("gGal4", 2),  # 2
    #     Gene("gSwi5", 3),  # 3
    #     Gene("gGal80",4), # 4
    #     Gene("gAsh1", 5)   # 5
    # ])

    # append!(species, [
    #     mRNA("mCbf1",  6,  species[1]),  # 6
    #     mRNA("mGal4",  7,  species[2]),  # 7
    #     mRNA("mSwi5",  8,  species[3]),  # 8
    #     mRNA("mGal80", 9, species[4]),  # 9
    #     mRNA("mAsh1", 10, species[5])   # 10
    # ])

    # append!(species, [
    #     Protein("Cbf1", 11, species[6]),  # 11
    #     Protein("Gal4", 12, species[7]),  # 12
    #     Protein("Swi5", 13, species[8]),  # 13
    #     Protein("Gal80",14, species[9]),  # 14
    #     Protein("Ash1", 15, species[10])  # 15
    # ])

    # append!(species, [
    #     GeneProtein("gp(gGal80-Swi5)", 16, species[4], species[13]),
    #     GeneProtein("gp(gASh1-Swi5)",  17, species[5], species[13]),
    #     GeneProtein("gp(gCbf1-Swi5)",  18, species[1], species[13]),
    #     GeneProtein("gp(gCbf1-Ash1)",  19, species[1], species[15]),
    #     GeneProtein("gp(gGal4-Cbf1)",  20, species[2], species[11]),
    #     GeneProtein("gp(gSwi5-Gal4)",  21, species[3], species[12]),
    #     ProteinComplex("c(Gal80-Gal4)",22, species[14],species[12])
    # ])
    # return species
end