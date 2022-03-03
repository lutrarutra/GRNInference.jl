abstract type Species end

mutable struct Gene <: Species
    name::String
    idx::Int64
    u₀::Float64
    gene_proteins::Vector{Species}

    num_params::Int64
    param_idx::Vector{Int64}

    Gene(name, idx::Int64) = new(name, idx, 0.05, [], 0, [])
end


function d_dt(gene::Gene)::Float64
    (
        sum(map(gp -> gp.u * gp.release_rate, gene.gene_proteins), init=0.0)
        - gene.u * sum(map(gp -> gp.protein.u * gp.binding_rate, gene.gene_proteins), init=0.0)
    )
end

function params!(gene::Gene, θ)
    # gene.u₀ = θ[1]
end

function params(gene::Gene)::Vector{Float64}
    []
end


mutable struct mRNA <: Species
    name::String
    idx::Int64
    u₀::Float64
    transcriptor::Species
    
    # params
    transcription_rate::Float64
    degradation_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}


    mRNA(name, idx::Int64, transcriptor::Species) =
        new(name, idx, 0.0, transcriptor, rand()*0.01, rand()*0.01, 2, [])

end


function d_dt(mrna::mRNA)
    (
        mrna.transcriptor.u * mrna.transcription_rate
        - mrna.u * mrna.degradation_rate
    )
end

function params!(m::mRNA, θ)
    m.transcription_rate = θ[1]
    m.degradation_rate = θ[2]
end

function params(m::mRNA)::Vector{Float64}
    [m.transcription_rate, m.degradation_rate]
end

mutable struct Protein <: Species
    name::String
    idx::Int64
    u₀::Float64
    targets::Vector{Species}
    translator::mRNA
    
    # params
    translation_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}

    Protein(name::String, idx::Int64, translator::mRNA) =
        new(name, idx, 0.0, [], translator, rand()*0.01, 1, [])
end

function d_dt(protein::Protein)
    (
        protein.translator.u * protein.translation_rate
        - protein.u * sum(map(target -> target.u * target.binding_rate, protein.targets), init=0.0)
        + sum(map(c -> c.u * c.release_rate, protein.targets), init=0.0)
    )
end

function params!(p::Protein, θ)
    p.translation_rate = θ[1]
end

function params(p::Protein)::Vector{Float64}
    [p.translation_rate]
end

mutable struct GeneProtein <: Species
    name::String
    idx::Int64
    u₀::Float64
    gene::Gene
    protein::Protein

    # params
    binding_rate::Float64
    release_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}
    

    function GeneProtein(name, idx::Int64, gene::Gene, protein::Protein)
        gp = new(name, idx, 0.0, gene, protein, rand()*0.01, rand()*0.01, 2, [])
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

function params!(gp::GeneProtein, θ)
    gp.binding_rate = θ[1]
    gp.release_rate = θ[2]
end

function params(gp::GeneProtein)::Vector{Float64}
    [gp.binding_rate, gp.release_rate]
end

mutable struct ProteinComplex <: Species
    name::String
    idx::Int64
    u₀::Float64
    protein1::Protein
    protein2::Protein

    # params
    binding_rate::Float64
    release_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}


    function ProteinComplex(name, idx::Int64, protein1::Protein, protein2::Protein)
        new(name, idx, 0.0, protein1, protein2, rand()*0.01, rand()*0.01, 2, [])
    end
end

function d_dt(c::ProteinComplex)
    (
        c.protein1.u * c.protein2.u * c.binding_rate
        - c.u * c.release_rate
    )
end

function params!(c::ProteinComplex, θ)
    c.binding_rate = θ[1]
    c.release_rate = θ[2]
end

function params(c::ProteinComplex)::Vector{Float64}
    [c.binding_rate, c.release_rate]
end


