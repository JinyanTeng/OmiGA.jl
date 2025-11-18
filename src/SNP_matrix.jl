export SNP_matrix
struct SNP_matrix <: AbstractMatrix{UInt8}
    data::Matrix{UInt8} 
    columncounts::Matrix{Int} 
    rowcounts::Matrix{Int} 
    m::Int 
    coding_rule::String         
end
Base.size(s::SNP_matrix) = s.m, size(s.data, 2) 
@inline function Base.getindex(f::SNP_matrix, i::Integer, j::Integer)
    @boundscheck checkbounds(f, i, j)
    ip3 = i + 3
    val = (f.data[ip3>>2, j] >> ((ip3 & 0x03) << 1)) & 0x03
    if f.coding_rule == "012"
        if val == 0x00
            0x02
        elseif val == 0x03
            0x00
        elseif val == 0x02
            0x01
        else  
            error("There are missing genotypes! Please provide the imputed genotype data.")
        end
    elseif f.coding_rule == "010"
        if val == 0x00 
            0x00
        elseif val == 0x03
            0x00
        elseif val == 0x02
            0x01
        else  
            error("There are missing genotypes! Please provide the imputed genotype data.")
        end
    end
end
function checkplinkfilename(filename::AbstractString, suffix::AbstractString)
    suffix ∈ ["bed", "fam", "bim"] || 
    throw(ArgumentError("suffix should be bed, fam or bim"))
    endswith(filename, "." * suffix) || 
    any(endswith.(filename, ".$suffix." .* ["gz", "zlib", "zz", "xz", "zst", "bz2"])) || 
    throw(ArgumentError("compressed format should be one of $ALLOWED_FORMAT"))
end
function SNP_matrix(
    bednm::AbstractString, 
    m::Integer,
    chooseID::Union{Nothing,AbstractArray{AbstractString}}=nothing, 
    args...;
    coding_rule::String="012", 
    kwargs...
    )
    checkplinkfilename(bednm, "bed") 
    data = open(bednm, args...; kwargs...) do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bednm"))
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bednm, is not in correct orientation"))
        Mmap.mmap(io) 
    end
    drows = (m + 3) >> 2
    n, r = divrem(length(data), drows)
    if isnothing(chooseID)
        SNP_matrix(reshape(data, (drows, n)), zeros(Int, (3, n)), zeros(Int, (3, m)), m, coding_rule)
    else
        @views IDname = CSV.read(string(splitext(bednm)[1], ".fam"), DataFrame, header=false)[:, 2]
        positionID = findall(in(chooseID), IDname)
        @views SNP_matrix(reshape(data, (drows, n)),zeros(Int, (4, n)), zeros(Int, (4, m)), m, coding_rule)[positionID, :]
    end
end
SNP_matrix(nm::AbstractString, args...; kwargs...) = SNP_matrix(nm, countlines(string(splitext(nm)[1], ".fam")), nothing, args...; coding_rule="012", kwargs...)
SNP_matrix(nm::AbstractString, m::Integer, args...; kwargs...) = SNP_matrix(nm, countlines(string(splitext(nm)[1], ".fam")), chooseID = nothing, args...;  coding_rule="012", kwargs...)
SNP_matrix(nm::AbstractString ,chooseID::AbstractArray{AbstractString}, args...;kwargs...) = SNP_matrix(nm, countlines(string(splitext(nm)[1], ".fam")),chooseID::AbstractArray{AbstractString}, args...;  coding_rule="012", kwargs...)
function _counts(s::AbstractMatrix{UInt8}, dims::Integer)
    if isone(dims)
        cc = (typeof(s) == SNP_matrix) ? s.columncounts : zeros(Int, (4, size(s, 2))) 
        if all(iszero, cc)
            m, n = size(s)
            @inbounds for j in 1:n
                for i in 1:m
                    cc[s[i, j] + 1, j] += 1  
                end
            end
        end
        return cc
    elseif dims == 2
        rc = (typeof(s) == SnpArray) ? s.rowcounts : zeros(Int, (4, size(s, 1)))
        if all(iszero, rc)
            m, n = size(s)
            @inbounds for j in 1:n
                for i in 1:m
                    rc[s[i, j] + 1, i] += 1
                end
            end
        end
        return rc
    else
        throw(ArgumentError("counts(s::SnpArray, dims=k) only defined for k = 1 or 2"))
    end
end        
_counts(s::AbstractMatrix{UInt8}, ::Colon) = sum(_counts(s, 1), dims=2) 
function maf!(out::AbstractVector{T}, s::AbstractMatrix{UInt8}) where T <: AbstractFloat
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2) 
        out[j] = (cc[2, j] + 2cc[3, j]) / 2(cc[1, j] + cc[2, j] + cc[3, j])
        (out[j] > 0.5) && (out[j] = 1 - out[j])
    end
    out
end
maf(s::AbstractMatrix{UInt8}) = maf!(Vector{Float64}(undef, size(s, 2)), s)
function af!(out::AbstractVector{T}, s::AbstractMatrix{UInt8}) where T <: AbstractFloat
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2) 
        out[j] = (cc[2, j] + 2cc[3, j]) / 2(cc[1, j] + cc[2, j] + cc[3, j])
        (out[j] < 0.5) && (out[j] = 1 - out[j])
    end
    out
end
af(s::AbstractMatrix{UInt8}) = af!(Vector{Float64}(undef, size(s, 2)), s) 
function het!(out::AbstractVector{T}, s::AbstractMatrix{UInt8}) where T <: AbstractFloat
    cc = _counts(s, 1)
    @inbounds for j in 1:size(s, 2) 
        out[j] = 2cc[2, j] / 2(cc[1, j] + cc[2, j] + cc[3, j])
    end
    out
end
het(s::AbstractMatrix{UInt8}) = het!(Vector{Float64}(undef, size(s, 2)), s) 
