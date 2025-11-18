function get_expvar_optimized!(M::AbstractMatrix{T1}, n_samples::Int, n_snps::Int, _2p::Matrix{T2}, _2pq::Matrix{T2}) where {T1, T2<:Float32}
    _p = zeros(Float32, 1, n_snps)
    _2nsamples = 1.0f0 / (2.0f0 * n_samples)
    sacle_vec = fill(_2nsamples, 1, n_samples)  
    mul!(_p, sacle_vec, M)
    mul!(_2p, _p, 2.0f0)  
    @. _p = 1.0f0 - _p
    @. _2pq = _2p * _p
end
function getG_fast(M::AbstractMatrix; alpha::Int=0, afs::Union{Nothing,AbstractArray{T}}=nothing, return_eltype::DataType=Float32) where {T<:AbstractFloat} 
    n_samples, n_snps = size(M)
    α = alpha
    _halfα = α / 2.0f0
    scalvar = 0.0f0  
    blocks = 1000  
    block_size = ceil(Int, n_snps / blocks)
    G = zeros(Float32, n_samples, n_samples)
    tempG = zeros(Float32, n_samples, n_samples)  
    X_block = zeros(Float32, n_samples, block_size)
    if α != 0.0f0
        _2pqα_block = zeros(Float32, 1, block_size)
    end
    _2p = zeros(Float32, 1, n_snps)
    _2pq = zeros(Float32, 1, n_snps)
    if isnothing(afs)
        get_expvar_optimized!(M, n_samples, n_snps, _2p, _2pq)
    else
        mul!(_2p, afs, 2.0f0)  
        @. _2pq[:] = 1.0f0 - afs
        @. _2pq *= _2p
    end
    if α == -1.0f0
        scalvar = n_snps
    elseif α == 0.0f0
        scalvar = sum(_2pq)
    else
        scalvar = sum(_2pq .^ (1.0f0 + α))
    end
    for b in 1:blocks
        i = (b - 1) * block_size + 1
        j = min(b * block_size, n_snps)  
        cblk_size = j - i + 1
        if cblk_size < block_size
            X_block = view(X_block, :, 1:cblk_size)
            if α != 0.0f0
                _2pqα_block = view(_2pqα_block, :, 1:cblk_size)
            end
        end
        M_block = view(M, :, i:j)
        _2p_block = view(_2p, :, i:j)
        _2pq_block = view(_2pq, :, i:j)
        broadcast!(-, X_block, M_block, _2p_block)
        if α != 0.0f0
            _2pqα_block .= _2pq_block .^ _halfα
            broadcast!(*, X_block, X_block, _2pqα_block)
        end
        mul!(tempG, X_block, X_block')
        tempG ./= scalvar
        G .+= tempG
    end
    G .= Symmetric(G)
    if return_eltype != Float32
        return return_eltype.(G)
    end
    return G
end
function getG_fast_threaded(M::Matrix{Int8}; alpha::Int=0, afs::Union{Nothing,AbstractArray{T}}=nothing, return_eltype::DataType=Float32) where {T<:AbstractFloat}
    n_samples, n_snps = size(M)
    α = Float32(alpha)
    _halfα = α / 2.0f0
    _2p = zeros(Float32, 1, n_snps)
    _2pq = zeros(Float32, 1, n_snps)
    if isnothing(afs)
        get_expvar_optimized!(M, n_samples, n_snps, _2p, _2pq)
    else
        mul!(_2p, afs, 2.0f0)
        @. _2pq[:] = 1.0f0 - afs
        @. _2pq *= _2p
    end
    scalvar = α == -1.0f0 ? n_snps :
              α == 0.0f0 ? sum(_2pq) :
              sum(_2pq .^ (1.0f0 + α))
    n_threads = Threads.nthreads()
    block_size = 1024
    blocks = ceil(Int, n_snps / block_size)
    buffer_pool = [
        (
            X_block = Matrix{Float32}(undef, n_samples, block_size),
            _2pqα_block = α != 0 ? Vector{Float32}(undef, block_size) : nothing,
            tempG = Matrix{Float32}(undef, n_samples, n_samples)
        ) for _ in 1:n_threads
    ]
    G = zeros(Float32, n_samples, n_samples)
    lock = Threads.SpinLock()
    Threads.@threads for b in 1:blocks
        tid = Threads.threadid()
        buf = buffer_pool[tid]
        i = (b-1)*block_size + 1
        j = min(b*block_size, n_snps)
        actual_size = j - i + 1
        X_block = view(buf.X_block, :, 1:actual_size)
        M_block = view(M, :, i:j)
        _2p_block = view(_2p, :, i:j)
        broadcast!(-, X_block, M_block, _2p_block)
        if α != 0
            _2pq_block = view(_2pq, :, i:j)
            _2pqα = view(buf._2pqα_block, 1:actual_size)
            @. _2pqα = _2pq_block ^ _halfα
            broadcast!(*, X_block, X_block, _2pqα')
        end
        mul!(buf.tempG, X_block, X_block')
        buf.tempG ./= scalvar
        Threads.lock(lock)
        G .+= buf.tempG
        Threads.unlock(lock)
    end
    G .= Symmetric(G)
    return_eltype != Float32 ? return_eltype.(G) : G
end
function getG_VanRaden(geno::Matrix; adjusted::Bool=false, afs::Union{Nothing,Vector{T}}=nothing, byrow::Bool=false, return_eltype::DataType=Float64) where {T <: Float64} 
    n_samples, n_snps = size(geno)[ifelse(byrow, [1, 2], [2, 1])]
    if isnothing(afs)
        if adjusted
            error("'afs' must be provided for pre-adjusted genotype data")
        end
        afs = vec(sum(geno, dims=ifelse(byrow, 1, 2))) / (2 * n_samples)
    end
    _2p = 2 .* afs
    sum_2pq = dot(_2p, (1 .- afs))
    G = zeros(T, n_samples, n_samples)
    if n_snps < 100000
        M = copy(geno)
        if byrow 
            if !adjusted
                M = M .- _2p'
            end
            mul!(G, M, M') 
        else 
            if !adjusted
                M = M .- _2p
            end
            mul!(G, M', M)  
        end
    else 
        block_size = 50000
        if byrow
            M = zeros(T, n_samples, block_size)
            iter_collects = collect(Iterators.partition(1:length(_2p), block_size))
            for i in 1:(length(iter_collects)-1)
                M .= geno[:, iter_collects[i]]
                M .-= _2p[iter_collects[i]]'
                G .+= mat_mul(M, M')
            end
            M = geno[:, iter_collects[end]] |> Matrix{T}
            M .-= _2p[iter_collects[end]]'
            G .+= mat_mul(M, M')
        else
            M = zeros(T, block_size, n_samples)
            iter_collects = collect(Iterators.partition(1:length(_2p), block_size))
            for i in 1:(length(iter_collects)-1)
                M .= geno[iter_collects[i], :]
                M .-= _2p[iter_collects[i]]
                G .+= mat_mul(M', M)
            end
            M = geno[iter_collects[end], :] |> Matrix{T}
            M .-= _2p[iter_collects[end]]
            G .+= mat_mul(M', M)
        end
    end
    G ./= sum_2pq
    G .= Symmetric(G)
    if return_eltype == Float32
        G = Float32.(G)
    end
    return G
end
function getG_VanRaden!(G::Matrix{T}, geno::Matrix{T}; adjusted::Bool=false, afs::Union{Nothing,Vector{T}}=nothing, byrow::Bool=false) where {T <: Float64} 
    n_samples, n_snps = size(geno)[ifelse(byrow, [1, 2], [2, 1])]
    if isnothing(afs)
        if adjusted
            error("'afs' must be provided for pre-adjusted genotype data")
        end
        afs = vec(sum(geno, dims=ifelse(byrow, 1, 2))) / (2 * n_samples)
    end
    _2p = 2 .* afs
    sum_2pq = dot(_2p, (1 .- afs))
    fill!(G, 0.0)
    if n_snps < 100000
        M = copy(geno)
        if byrow 
            if !adjusted
                M = M .- _2p'
            end
            mul!(G, M, M')
        else 
            if !adjusted
                M = M .- _2p
            end
            mul!(G, M', M)  
        end
    else 
        block_size = 50000
        if byrow
            M = zeros(T, n_samples, block_size)
            iter_collects = collect(Iterators.partition(1:length(_2p), block_size))
            for i in 1:(length(iter_collects)-1)
                M .= geno[:, iter_collects[i]]
                M .-= _2p[iter_collects[i]]'
                G .+= mat_mul(M, M')
            end
            M = geno[:, iter_collects[end]] |> Matrix{T}
            M .-= _2p[iter_collects[end]]'
            G .+= mat_mul(M, M')
        else
            M = zeros(T, block_size, n_samples)
            iter_collects = collect(Iterators.partition(1:length(_2p), block_size))
            for i in 1:(length(iter_collects)-1)
                M .= geno[iter_collects[i], :]
                M .-= _2p[iter_collects[i]]
                G .+= mat_mul(M', M)
            end
            M = geno[iter_collects[end], :] |> Matrix{T}
            M .-= _2p[iter_collects[end]]
            G .+= mat_mul(M', M)
        end
    end
    G ./= sum_2pq
    G .= Symmetric(G)
    return G
end
function getG_dominance(geno::Matrix; code_type::String, adjusted::Bool=false, afs::Union{Nothing,Vector{T}}=nothing, byrow::Bool=false, return_eltype::DataType=Float64) where {T <: Float64}
    if code_type == "additive"
        if adjusted
            error("genotype data must be coded as 0,1,2")
        end
        if isnothing(afs)
            afs = vec(sum(geno, dims=ifelse(byrow, 1, 2))) / (2 * n_samples)
        end
        S = replace(geno, 2 => 0)
    elseif code_type == "dominance"
        S = copy(geno)
        if isnothing(afs)
            error("'afs' must be provided for dominance code genotype data")
        end
    end
    _2pq = 2 .* afs .* (1 .- afs)
    scale_factor = dot(_2pq, 1 .- _2pq)
    n_samples, n_snps = size(geno)[ifelse(byrow, [1, 2], [2, 1])]
    G = zeros(T, n_samples, n_samples)
    if n_snps < 100000
        if byrow 
            if !adjusted
                S = S .- _2pq'
            end
            mul!(G, S, S')
        else 
            if !adjusted
                S = S .- _2pq
            end
            mul!(G, S', S)  
        end
    else 
        block_size = 50000
        if byrow
            S = zeros(T, n_samples, block_size)
            iter_collects = collect(Iterators.partition(1:length(_2pq), block_size))
            for i in 1:(length(iter_collects)-1)
                S .= geno[:, iter_collects[i]]
                S .-= _2pq[iter_collects[i]]'
                G .+= mat_mul(S, S')
            end
            S = geno[:, iter_collects[end]] |> Matrix{T}
            S .-= _2pq[iter_collects[end]]'
            G .+= mat_mul(S, S')
        else
            S = zeros(T, block_size, n_samples)
            iter_collects = collect(Iterators.partition(1:length(_2pq), block_size))
            for i in 1:(length(iter_collects)-1)
                S .= geno[iter_collects[i], :]
                S .-= _2pq[iter_collects[i]]
                G .+= mat_mul(S', S)
            end
            S = geno[iter_collects[end], :] |> Matrix{T}
            S .-= _2pq[iter_collects[end]]
            G .+= mat_mul(S', S)
        end
    end
    G ./= scale_factor
    G .= Symmetric(G)
    if return_eltype == Float32
        G = Float32.(G)
    end
    return G
end
function getG_dominance!(G::Matrix{T}, geno::Matrix{T}; code_type::String, adjusted::Bool=false, afs::Union{Nothing,Vector{T}}=nothing, byrow::Bool=false) where {T <: Float64}
    if code_type == "additive"
        if adjusted
            error("genotype data must be coded as 0,1,2")
        end
        if isnothing(afs)
            afs = vec(sum(geno, dims=ifelse(byrow, 1, 2))) / (2 * n_samples)
        end
        S = replace(geno, 2 => 0)
    elseif code_type == "dominance"
        S = copy(geno)
        if isnothing(afs)
            error("'afs' must be provided for dominance code genotype data")
        end
    end
    _2pq = 2 .* afs .* (1 .- afs)
    scale_factor = dot(_2pq, 1 .- _2pq)
    n_samples, n_snps = size(geno)[ifelse(byrow, [1, 2], [2, 1])]
    fill!(G, 0.0)
    if n_snps < 100000
        if byrow 
            if !adjusted
                S = S .- _2pq'
            end
            mul!(G, S, S')
        else 
            if !adjusted
                S = S .- _2pq
            end
            mul!(G, S', S)  
        end
    else 
        block_size = 50000
        if byrow
            S = zeros(T, n_samples, block_size)
            iter_collects = collect(Iterators.partition(1:length(_2pq), block_size))
            for i in 1:(length(iter_collects)-1)
                S .= geno[:, iter_collects[i]]
                S .-= _2pq[iter_collects[i]]'
                G .+= mat_mul(S, S')
            end
            S = geno[:, iter_collects[end]] |> Matrix{T}
            S .-= _2pq[iter_collects[end]]'
            G .+= mat_mul(S, S')
        else
            S = zeros(T, block_size, n_samples)
            iter_collects = collect(Iterators.partition(1:length(_2pq), block_size))
            for i in 1:(length(iter_collects)-1)
                S .= geno[iter_collects[i], :]
                S .-= _2pq[iter_collects[i]]
                G .+= mat_mul(S', S)
            end
            S = geno[iter_collects[end], :] |> Matrix{T}
            S .-= _2pq[iter_collects[end]]
            G .+= mat_mul(S', S)
        end
    end
    G ./= scale_factor
    G .= Symmetric(G)
    return G
end
