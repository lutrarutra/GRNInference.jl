include("Species.jl")

function simplified_network()
    net = BioNetwork()
    net.desc = "simplified_network"

    add_species!(net, Gene, "gCbf1")
    add_species!(net, Gene, "gGal4")
    add_species!(net, Gene, "gSwi5")
    add_species!(net, Gene, "gGal80")
    add_species!(net, Gene, "gAsh1")

    add_species!(net, mRNA, "mCbf1")
    add_species!(net, mRNA, "mGal4")
    add_species!(net, mRNA, "mSwi5")
    add_species!(net, mRNA, "mGal80")
    add_species!(net, mRNA, "mAsh1")

    add_species!(net, Protein, "Swi5")
    add_species!(net, Protein, "Cbf1")
    add_species!(net, Protein, "Gal4")
    add_species!(net, Protein, "Gal80")
    add_species!(net, Protein, "Ash1")

    add_species!(net, Complex, "gp(gGal80-Swi5)")
    add_species!(net, Complex, "gp(gASh1-Swi5)")
    add_species!(net, Complex, "gp(gCbf1-Swi5)")
    add_species!(net, Complex, "gp(gCbf1-Ash1)")
    add_species!(net, Complex, "gp(gGal4-Cbf1)")
    add_species!(net, Complex, "gp(gSwi5-Gal4)")
    add_species!(net, Complex, "c(Gal80-Gal4)")

    add_transcriptor!(net["mCbf1"], net["gCbf1"])
    add_transcriptor!(net["mGal4"], net["gGal4"])
    add_transcriptor!(net["mSwi5"], net["gSwi5"])
    add_transcriptor!(net["mGal80"],net["gGal80"])
    add_transcriptor!(net["mAsh1"], net["gAsh1"])

    add_translator!(net["Swi5"], net["mSwi5"])
    add_translator!(net["Cbf1"], net["mCbf1"])
    add_translator!(net["Gal4"], net["mGal4"])
    add_translator!(net["Gal80"],net["mGal80"])
    add_translator!(net["Ash1"], net["mAsh1"])

    add_reactants!(net["gp(gGal80-Swi5)"], net["gGal80"], net["Swi5"])
    add_reactants!(net["gp(gASh1-Swi5)"], net["gAsh1"], net["Swi5"])
    add_reactants!(net["gp(gCbf1-Swi5)"], net["gCbf1"], net["Swi5"])
    add_reactants!(net["gp(gCbf1-Ash1)"], net["gCbf1"], net["Ash1"])
    add_reactants!(net["gp(gGal4-Cbf1)"], net["gGal4"], net["Cbf1"])
    add_reactants!(net["gp(gSwi5-Gal4)"], net["gSwi5"], net["Gal4"])
    add_reactants!(net["c(Gal80-Gal4)"], net["Gal80"], net["Gal4"])

    return net
end

function true_network()
    net = BioNetwork()
    net.desc = "true_network"

    add_species!(net, Gene, "gCbf1")
    add_species!(net, Gene, "gGal4")
    add_species!(net, Gene, "gSwi5")
    add_species!(net, Gene, "gGal80")
    add_species!(net, Gene, "gAsh1")

    add_species!(net, mRNA, "mCbf1")
    add_species!(net, mRNA, "mGal4")
    add_species!(net, mRNA, "mSwi5")
    add_species!(net, mRNA, "mGal80")
    add_species!(net, mRNA, "mAsh1")

    add_species!(net, Protein, "Swi5")
    add_species!(net, Protein, "Cbf1")
    add_species!(net, Protein, "Gal4")
    add_species!(net, Protein, "Gal80")
    add_species!(net, Protein, "Ash1")

    add_species!(net, Complex, "gp(gGal80-Swi5)")
    add_species!(net, Complex, "gp(gASh1-Swi5)")
    add_species!(net, Complex, "gp(gCbf1-Swi5)")
    add_species!(net, Complex, "gp(gCbf1-Ash1)")
    add_species!(net, Complex, "gp(gGal4-Cbf1)")
    add_species!(net, Complex, "gp(gSwi5-Gal4)")
    add_species!(net, Complex, "c(Gal80-Gal4)")

    add_transcriptor!(net["mCbf1"], net["gp(gCbf1-Swi5)"])
    add_transcriptor!(net["mGal4"], net["gGal4"])
    add_transcriptor!(net["mSwi5"], net["gp(gSwi5-Gal4)"])
    add_transcriptor!(net["mGal80"],net["gp(gGal80-Swi5)"])
    add_transcriptor!(net["mAsh1"], net["gp(gASh1-Swi5)"])

    add_translator!(net["Swi5"], net["mSwi5"])
    add_translator!(net["Cbf1"], net["mCbf1"])
    add_translator!(net["Gal4"], net["mGal4"])
    add_translator!(net["Gal80"],net["mGal80"])
    add_translator!(net["Ash1"], net["mAsh1"])

    add_reactants!(net["gp(gGal80-Swi5)"], net["gGal80"], net["Swi5"])
    add_reactants!(net["gp(gASh1-Swi5)"], net["gAsh1"], net["Swi5"])
    add_reactants!(net["gp(gCbf1-Swi5)"], net["gCbf1"], net["Swi5"])
    add_reactants!(net["gp(gCbf1-Ash1)"], net["gCbf1"], net["Ash1"])
    add_reactants!(net["gp(gGal4-Cbf1)"], net["gGal4"], net["Cbf1"])
    add_reactants!(net["gp(gSwi5-Gal4)"], net["gSwi5"], net["Gal4"])
    add_reactants!(net["c(Gal80-Gal4)"], net["Gal80"], net["Gal4"])

    return net
end

function random_network(num_interactions::Int64)
    net = BioNetwork()
    net.desc = "random_network_$num_interactions"

    add_species!(net, Gene, "gCbf1")
    add_species!(net, Gene, "gGal4")
    add_species!(net, Gene, "gSwi5")
    add_species!(net, Gene, "gGal80")
    add_species!(net, Gene, "gAsh1")

    add_species!(net, Protein, "Swi5")
    add_species!(net, Protein, "Cbf1")
    add_species!(net, Protein, "Gal4")
    add_species!(net, Protein, "Gal80")
    add_species!(net, Protein, "Ash1")

    used = []
    while length(used) < num_interactions
        s1_idx = rand(1:10)
        s2_idx = rand(6:10)
        if s1_idx == s2_idx || (s1_idx, s2_idx) ∈ used
            continue
        end
        push!(used, (s1_idx, s2_idx))
        s1 = net[s1_idx]
        s2 = net[s2_idx]
        name = "c($(s1.name)-$(s2.name))"
        add_species!(net, Complex, name)
        add_reactants!(net[name], s1, s2)
    end

    add_species!(net, mRNA, "mCbf1")
    add_species!(net, mRNA, "mGal4")
    add_species!(net, mRNA, "mSwi5")
    add_species!(net, mRNA, "mGal80")
    add_species!(net, mRNA, "mAsh1")

    add_transcriptor!(net["mCbf1"], net["gCbf1"])
    add_transcriptor!(net["mGal4"], net["gGal4"])
    add_transcriptor!(net["mSwi5"], net["gSwi5"])
    add_transcriptor!(net["mGal80"],net["gGal80"])
    add_transcriptor!(net["mAsh1"], net["gAsh1"])

    add_translator!(net["Swi5"], net["mSwi5"])
    add_translator!(net["Cbf1"], net["mCbf1"])
    add_translator!(net["Gal4"], net["mGal4"])
    add_translator!(net["Gal80"],net["mGal80"])
    add_translator!(net["Ash1"], net["mAsh1"])

    return net
end