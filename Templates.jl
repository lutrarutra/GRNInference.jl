function true_network()
    species::Vector{Species} = []
    append!(species, [
        Gene("gCbf1"),  # 1
        Gene("gGal4"),  # 2
        Gene("gSwi5"),  # 3
        Gene("gGal80"), # 4
        Gene("gAsh1")   # 5
    ])

    append!(species, [
        mRNA("mCbf1",  species[1]),  # 6
        mRNA("mGal4",  species[2]),  # 7
        mRNA("mSwi5",  species[3]),  # 8
        mRNA("mGal80", species[4]),  # 9
        mRNA("mAsh1",  species[5])   # 10
    ])

    append!(species, [
        Protein("Cbf1", species[6]),  # 11
        Protein("Gal4", species[7]),  # 12
        Protein("Swi5", species[8]),  # 13
        Protein("Gal80",species[9]),  # 14
        Protein("Ash1", species[10])  # 15
    ])

    append!(species, [
        GeneProtein("gp(gGal80-Swi5)", species[4], species[13]),
        GeneProtein("gp(gASh1-Swi5)",  species[5], species[13]),
        GeneProtein("gp(gCbf1-Swi5)",  species[1], species[13]),
        GeneProtein("gp(gCbf1-Ash1)",  species[1], species[15]),
        GeneProtein("gp(gGal4-Cbf1)",  species[2], species[11]),
        GeneProtein("gp(gSwi5-Gal4)",  species[3], species[12]),
        ProteinComplex("c(Gal80-Gal4)", species[14],species[12])
    ])
    return species
end