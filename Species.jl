abstract type Species end

mutable struct Gene <: Species
    name::String
    u₀::Float64
    u::Float64
    gene_proteins::Vector{Species}

    function Gene(name) 
        u::Float64 = rand()
        new(name, u, u, [])
    end
end

function d_dt(gene::Gene)::Float64
    (
        sum(map(gp -> gp.u * gp.release_rate, gene.gene_proteins), init=0.0)
        - gene.u * sum(map(gp -> gp.protein.u * gp.binding_rate, gene.gene_proteins), init=0.0)
    )
end

function θ!(gene::Gene, vals::Vector{Float64})
    gene.u₀ = vals[1]
end

function θ(gene::Gene)::Vector{Float64}
    [gene.u₀]
end


mutable struct mRNA <: Species
    name::String
    u₀::Float64
    u::Float64
    transcriptor::Species

    # params
    transcription_rate::Float64
    degradation_rate::Float64

    mRNA(name, transcriptor::Species) =
        new(name, 0.0, 0.0, transcriptor, rand(), rand())

end


function d_dt(mrna::mRNA)
    (
        mrna.transcriptor.u * mrna.transcription_rate
        - mrna.u * mrna.degradation_rate
    )
end

function θ!(m::mRNA, vals::Vector{Float64})
    m.transcription_rate = vals[1]
    m.degradation_rate = vals[2]
end

function θ(m::mRNA)::Vector{Float64}
    [m.transcription_rate, m.degradation_rate]
end

mutable struct Protein <: Species
    name::String
    u₀::Float64
    u::Float64
    targets::Vector{Species}
    translator::mRNA
    
    # params
    translation_rate::Float64

    Protein(name, translator::mRNA) =
        new(name, 0.0, 0.0, [], translator, rand())
end

function d_dt(protein::Protein)
    (
        protein.translator.u * protein.translation_rate
        - protein.u * sum(map(target -> target.u * target.binding_rate, protein.targets), init=0.0)
        + sum(map(c -> c.u * c.release_rate, protein.targets), init=0.0)
    )
end

function θ!(p::Protein, vals::Vector{Float64})
    p.translation_rate = vals[1]
end

function θ(p::Protein)::Vector{Float64}
    [p.translation_rate]
end

mutable struct GeneProtein <: Species
    name::String
    u₀::Float64
    u::Float64
    gene::Gene
    protein::Protein

    # params
    binding_rate::Float64
    release_rate::Float64

    function GeneProtein(name, gene::Gene, protein::Protein)
        gp = new(name, 0.0, 0.0, gene, protein, rand(), rand())
        push!(protein.targets, gp)
        push!(gene.gene_proteins, gp)
        gp
    end
end

function d_dt(gp::GeneProtein)
    (
        gp.protein.u * gp.gene.u * gp.binding_rate
        - gp.u * gp.release_rate
    )
end

function θ!(gp::GeneProtein, vals::Vector{Float64})
    gp.binding_rate = vals[1]
    gp.release_rate = vals[2]
end

function θ(gp::GeneProtein)::Vector{Float64}
    [gp.binding_rate, gp.release_rate]
end

mutable struct ProteinComplex <: Species
    name::String
    u₀::Float64
    u::Float64
    protein1::Protein
    protein2::Protein

    # params
    binding_rate::Float64
    release_rate::Float64

    function ProteinComplex(name, protein1::Protein, protein2::Protein)
        new(name, 0.0, 0.0, protein1, protein2, rand(), rand())
    end
end

function d_dt(c::ProteinComplex)
    (
        c.protein1.u * c.protein2.u * c.binding_rate
        - c.u * c.release_rate
    )
end

function θ!(c::ProteinComplex, vals::Vector{Float64})
    c.binding_rate = vals[1]
    c.release_rate = vals[2]
end

function θ(c::ProteinComplex)::Vector{Float64}
    [c.binding_rate, c.release_rate]
end


