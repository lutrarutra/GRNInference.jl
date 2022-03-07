abstract type Species end

toM(x) = mapreduce(permutedims, vcat, x)

mutable struct Gene <: Species
    name::String
    idx::Int64
    u₀::Float64
    num_params::Int64

    Gene(name, idx::Int64) = new(name, idx, 0.05, 0)
end

csv(g::Gene) = "Gene,$(g.idx),$(g.name),$(g.u₀)"

mutable struct mRNA <: Species
    name::String
    idx::Int64
    u₀::Float64
    num_params::Int64

    transcriptor::Species

    mRNA(name, idx::Int64) = new(name, idx, 0.05, 2)
end

function add_transcriptor!(m::mRNA, s::Species)
    m.transcriptor = s
end

csv(m::mRNA) = "mRNA,$(m.idx),$(m.name),$(m.u₀),$(m.transcriptor.name)"

mutable struct Protein <: Species
    name::String
    idx::Int64
    u₀::Float64
    num_params::Int64

    translator::mRNA

    Protein(name::String, idx::Int64) = new(name, idx, 0.0, 1)
end

function add_translator!(p::Protein, t::mRNA)
    p.translator = t
end

csv(p::Protein) = "Protein,$(p.idx),$(p.name),$(p.u₀),$(p.translator.name)"

mutable struct Complex <: Species
    name::String
    idx::Int64
    u₀::Float64
    num_params::Int64

    species1::Species
    species2::Species
    
    Complex(name, idx::Int64) = new(name, idx, 0.0, 2)
end

function add_reactants!(c::Complex, s1::Species, s2::Species)
    c.species1 = s1
    c.species2 = s2
end

csv(c::Complex) = "Complex,$(c.idx),$(c.name),$(c.u₀),$(c.species1.name),$(c.species2.name)"





