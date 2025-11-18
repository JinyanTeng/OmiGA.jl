function append_top_summary!(_df_full, _df_tops, gene_i, n_tests)
    if n_tests == 1
        index_top = first_nonnan_argmin(_df_full.pval_g1) 
        if index_top == -1
            return false
        end
        _df_tops.num_var[gene_i] = size(_df_full, 1)
        _df_tops.variant_id[gene_i] = _df_full.variant_id[index_top]
        _df_tops.start_distance[gene_i] = _df_full.start_distance[index_top]
        _df_tops.af[gene_i] = _df_full.af[index_top]
        _df_tops.beta_g1[gene_i] = _df_full.beta_g1[index_top]
        _df_tops.beta_se_g1[gene_i] = _df_full.beta_se_g1[index_top]
        _df_tops.pval_g1[gene_i] = _df_full.pval_g1[index_top]
    elseif n_tests == 2
        index_top = first_nonnan_argmin(_df_full.pval_g2)
        if index_top == -1
            return false
        end
        _df_tops.num_var[gene_i] = size(_df_full, 1)
        _df_tops.variant_id[gene_i] = _df_full.variant_id[index_top]
        _df_tops.start_distance[gene_i] = _df_full.start_distance[index_top]
        _df_tops.af[gene_i] = _df_full.af[index_top]
        _df_tops.beta_g1[gene_i] = _df_full.beta_g1[index_top]
        _df_tops.beta_se_g1[gene_i] = _df_full.beta_se_g1[index_top]
        _df_tops.pval_g1[gene_i] = _df_full.pval_g1[index_top]
        _df_tops.beta_g2[gene_i] = _df_full.beta_g2[index_top]
        _df_tops.beta_se_g2[gene_i] = _df_full.beta_se_g2[index_top]
        _df_tops.pval_g2[gene_i] = _df_full.pval_g2[index_top]
    elseif n_tests > 2
        index_top = first_nonnan_argmin(_df_full.pval_joint)
        if index_top == -1
            return false
        end
        _df_tops.num_var[gene_i] = size(_df_full, 1)
        _df_tops.variant_id[gene_i] = _df_full.variant_id[index_top]
        _df_tops.start_distance[gene_i] = _df_full.start_distance[index_top]
        _df_tops.af[gene_i] = _df_full.af[index_top]
        _df_tops[gene_i, string.("beta_g", 1:n_tests)] = _df_full[index_top, string.("beta_g", 1:n_tests)]
        _df_tops[gene_i, string.("beta_se_g", 1:n_tests)] = _df_full[index_top, string.("beta_se_g", 1:n_tests)]
        _df_tops[gene_i, string.("pval_g", 1:n_tests)] = _df_full[index_top, string.("pval_g", 1:n_tests)]
        _df_tops.pval_joint[gene_i] = _df_full.pval_joint[index_top]
    end
    return true
end
function ACATest(Pvals::AbstractVector{T1}; Weights::Union{Nothing,T2}=nothing, is_check::Bool=true) where {T1<:AbstractFloat,T2<:AbstractArray}
    if is_check
        if any(ismissing.(Pvals))
            error("Cannot have NAs in the p-values!")
        end
        if any(Pvals .< 0) || any(Pvals .> 1)
            error("P-values must be between 0 and 1!")
        end
        has_zero = any(Pvals .== 0)
        has_one = any(Pvals .== 1)
        if has_zero && has_one
            error("Cannot have both 0 and 1 p-values!")
        end
        if has_zero
            return 0.0
        end
        if has_one
            @warn "There are p-values that are exactly 1!"
            return 1.0
        end
    end
    n = length(Pvals)
    if isnothing(Weights)
        Weights = fill(1.0 / n, n)  
    else
        if length(Weights) != n
            error("The length of weights should be the same as that of the p-values")
        end
        if any(Weights .< 0)
            error("All the weights must be positive!")
        end
        Weights = Weights / sum(Weights)
    end
    is_small = Pvals .< 1e-16
    cct_stat = 0.0
    if sum(is_small) == 0
        cct_stat = sum(Weights .* tan.((0.5 .- Pvals) .* π))
    else
        cct_stat += sum((Weights[is_small] ./ Pvals[is_small]) ./ π)
        cct_stat += sum(Weights[.!is_small] .* tan.((0.5 .- Pvals[.!is_small]) .* π))
    end
    pval = if cct_stat > 1e15
        (1.0 / cct_stat) / π
    else
        1.0 - cdf(Cauchy(), cct_stat)
    end
    return pval
end
function get_GC_lambda(median_P::T) where {T}
    if isnan(median_P)
        lambda = NaN
    else
        lambda = quantile(Chisq(1), 1 - median_P) / 0.4549364
    end
    return lambda
end
function get_GC_lambda(median_P::T, _ni, _X_c) where {T}
    if isnan(median_P)
        lambda = NaN
    else
        d = FDist(1, _ni - _X_c - 1)
        lambda = quantile(d, 1 - median_P) / quantile(d, 0.5)
    end
    return lambda
end
function WaldTest(v1::Number, v2::Number)
    if _args_wald_test == "FDist"
        return FDist(v1, v2)
    else
        return Chisq(v1)
    end
end
function residualizer(X_MME::Matrix{T}, center=true) where {T}
    C_t = copy(X_MME) 
    if center
        C_t = C_t .- Statistics.mean(C_t, dims=1)
    end
    Q, _ = qr(C_t)
    Q_t = Q[:, 1:size(C_t, 2)]
    dof = size(X_MME, 1) - 2 - size(X_MME, 2)
    return Q_t, dof
end
function residualizer(X_MME::AbstractMatrix{T}, center=true) where {T}
    C_t = copy(X_MME) 
    if center
        C_t = C_t .- Statistics.mean(C_t, dims=1)
    end
    Q, _ = qr(C_t)
    @runif USE_GPU Q = Q |> CuArray
    Q_t = Q[:, 1:size(C_t, 2)]
    dof = size(X_MME, 1) - 2 - size(X_MME, 2)
    return Q_t, dof
end
function residualizer_transform(M_t::Vector, Q::AbstractMatrix; center::Bool=true)
    M0_t = M_t .- Statistics.mean(M_t)
    if center
        M0_t .-= Q * Q' * M0_t
    else
        M0_t .= M_t .- Q * Q' * M0_t
    end
    return M0_t
end
function residualizer_transform(M_t::AbstractMatrix, Q::AbstractMatrix; center=true)
    M0_t = M_t .- Statistics.mean(M_t, dims=1)
    if center
        M0_t .-= Q * Q' * M0_t 
    else
        M0_t = M_t .- Q * Q' * M0_t
    end
    return M0_t
end
function residualizer_transform!(M_t::AbstractMatrix, Q::AbstractMatrix; center=true)
    if center
        broadcast!(-, M_t, M_t, Statistics.mean(M_t, dims=1))
    end
    M_t .-= Q * Q' * M_t 
end
function center_normalize(M_t::AbstractMatrix{T}, dim::Int=1) where {T}
    """Center and normalize M"""
    N_t = M_t - Statistics.mean(M_t, dims=dim)
    return N_t ./ sqrt.(sum(abs2.(N_t), dims=dim))
end
function center_normalize!(M_t::AbstractMatrix{T}, dim::Int=1; center::Bool=true) where {T}
    """Center and normalize M"""
    if center
        M_t .-= Statistics.mean(M_t, dims=dim)
    end
    M_t ./= sqrt.(sum(abs2.(M_t), dims=dim))
end
function center_normalize!(M_t::AbstractVector{T}, dim::Int=1; center::Bool=true) where {T}
    """Center and normalize M"""
    if center
        M_t .-= Statistics.mean(M_t)
    end
    M_t ./= sqrt(sum(abs2.(M_t)))
end
function get_matrix_resid(exppheno::AbstractVector{T}, Q::Union{Nothing,AbstractMatrix{T}}; return_std::Bool=false) where {T}
    if !isnothing(Q)
        exppheno_res = residualizer_transform(exppheno, Q; center=true)
    else
        exppheno_res = copy(exppheno)
    end
    if return_std
        exppheno_std = std(exppheno_res)
    end
    phenotype_res_t = exppheno_res
    center_normalize!(phenotype_res_t)
    if return_std
        return phenotype_res_t, exppheno_std
    else
        return phenotype_res_t
    end
end
function get_matrix_resid(exppheno::AbstractMatrix{T}, Q::Union{Nothing,AbstractMatrix{T}}; return_std::Bool=false) where {T}
    if !isnothing(Q)
        exppheno_res = residualizer_transform(exppheno, Q; center=true)
    else
        exppheno_res = copy(exppheno)
    end
    if return_std
        exppheno_std = vec(std(exppheno_res, dims=1))
    end
    phenotype_res_t = exppheno_res
    center_normalize!(phenotype_res_t)
    if return_std
        return phenotype_res_t, exppheno_std
    else
        return phenotype_res_t
    end
end
function get_matrix_resid!(exppheno::AbstractMatrix{T}, Q::Union{Nothing,AbstractMatrix{T}}; return_std::Bool=false, center::Bool=true) where {T}
    @timeit to "residualizer_transform!" if !isnothing(Q)
        residualizer_transform!(exppheno, Q; center=center)
    end
    @timeit to "std" if return_std
        exppheno_std = vec(std(exppheno, dims=1))
    end
    @timeit to "center_normalize!" center_normalize!(exppheno;center=center)
    if return_std
        return exppheno_std
    end
end
function get_phenotype_resid(exppheno::Vector{T}, I_minus_Q_Qt::Matrix{T}; byrow::Bool=false) where {T}
    y_minus_mean = exppheno .- Statistics.mean(exppheno)
    phenotype_res_t = mat_mul(I_minus_Q_Qt, y_minus_mean)
    phenotype_res_t ./= norm(phenotype_res_t)
    return phenotype_res_t
end
function get_phenotype_resid!(y_perms_res::Matrix{T}, y_perms::Matrix{T}, I_minus_Q_Qt::Matrix{T}; byrow::Bool=false) where {T}
    mat_mul!(y_perms_res, I_minus_Q_Qt, y_perms)
    y_perms_res ./= sqrt.(sum(abs2.(y_perms_res), dims=1))
end
function get_genotype_residual!(genotype_resid::Matrix{T}, genotype::Matrix{Int8}, afs::Vector{T}, I_minus_Q_Qt::Matrix{T}, dominance::Bool=false) where {T} 
    if !dominance
        genotype_resid .= genotype .- 2 .* afs
    else
        genotype_resid .= genotype .- 2 .* afs .* (1 .- afs)
    end
    genotype_resid .= genotype_resid * I_minus_Q_Qt
    sqrt_sum_abs2_eachrow = similar(genotype_resid, size(genotype_resid, 1))
    @inbounds @views @threads for i in 1:size(genotype_resid, 1)
        sqrt_sum_abs2_eachrow[i] = norm(genotype_resid[i, :])
    end
    genotype_resid ./= sqrt_sum_abs2_eachrow
end
function get_genotype_residual!(genotype_resid::Matrix{T}, genotype::Matrix{T}, I_minus_Q_Qt::Matrix{T}; byrow::Bool=true) where {T}
    if byrow
        genotype_resid .= I_minus_Q_Qt * genotype
        sqrt_sum_abs2_eachcol = similar(genotype_resid, 1, size(genotype_resid, 2))
        @inbounds @views @threads for i in 1:size(genotype_resid, 2)
            sqrt_sum_abs2_eachcol[i] = norm(genotype_resid[:, i])
        end
        genotype_resid ./= sqrt_sum_abs2_eachcol
    else
        genotype_resid .= genotype * I_minus_Q_Qt
        sqrt_sum_abs2_eachrow = similar(genotype_resid, size(genotype_resid, 1))
        @inbounds @views @threads for i in 1:size(genotype_resid, 1)
            sqrt_sum_abs2_eachrow[i] = norm(genotype_resid[i, :])
        end
        genotype_resid ./= sqrt_sum_abs2_eachrow
    end
end
function get_genotype_residual(genotype::Matrix{T}, I_minus_Q_Qt::Matrix{T}; byrow::Bool=false) where {T}
    if byrow
        genotype_resid = genotype' * I_minus_Q_Qt
    else
        genotype_resid = genotype * I_minus_Q_Qt
    end
    copyto!(genotype_resid, genotype_resid ./ sqrt.(sum(abs2, genotype_resid, dims=2)))
    if byrow
        genotype_resid = genotype_resid' |> Matrix
    end
    return genotype_resid
end
function getVinv(EA::Eigen{T,T,Matrix{T},Vector{T}}, Σ_i::Vector{T}) where {T}
    Dinv = spdiagm(1.0 ./ (Σ_i[1] .* EA.values .+ Σ_i[2]))
    Vi = EA.vectors * Dinv * EA.vectors'
    return Vi
end
function getVinv!(Vi::Matrix{T}, Di::Diagonal{T,Vector{T}}, EA::Eigen{T,T,Matrix{T},Vector{T}}, Σ_i::Vector{T}, EAvec_Di::Matrix{T}) where {T}
    Di[diagind(Di)] .= 1.0 ./ (Σ_i[1] .* EA.values .+ Σ_i[2])
    mul!(EAvec_Di, Di, EA.vectors')
    mat_mul!(Vi, EA.vectors, EAvec_Di)
end
function getVinv!(Vi::Matrix{T1}, Σ_i::Vector{T1}, partialV::Union{Nothing,PartialV{T2}}, Di::Diagonal{T1,Vector{T1}}, EAvec_Di::Matrix{T1}, EAvec::Matrix{T1}) where {T1<:AbstractFloat,T2<:AbstractFloat}
    ratio = maximum([Σ_i[1] / Σ_i[2], Σ_i[2] / Σ_i[1]])
    println("Ratio between vg1 and vg2: ", ratio)
    _n_cmp = length(Σ_i)
    ratio_index = argmin(abs.(ratio .- partialV.vec_ratio)) * 2
    if Σ_i[1] >= Σ_i[2]
        ratio_index -= 1
    end
    j = ratio_index
    n_samples = size(Vi, 1)
    range_inds = n_samples*(j-1)+1:n_samples*j
    vals_inds = CartesianIndices((range_inds,))
    vecs_inds = CartesianIndices((range_inds, 1:n_samples))
    copyto!(EAvec, CartesianIndices((1:n_samples, 1:n_samples)), partialV.eigvecs, vecs_inds)
    begin
        Di[diagind(Di)] .= partialV.eigvals[vals_inds]
        Di[diagind(Di)] .*= sum(Σ_i[1:2])
        Di[diagind(Di)] .+= Σ_i[3]
        Di[diagind(Di)] .\= 1.0
    end
    mul!(EAvec_Di, Di, EAvec')
    mat_mul!(Vi, EAvec, EAvec_Di)
end
function getVinv!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, Σ_i::Vector{T}) where {T}
    _n_cmp = length(Σ_i)
    fill!(Vi, 0.0)
    Vi[diagind(Vi)] .= Σ_i[_n_cmp]
    for i in 1:(_n_cmp-1)
        Vi .+= VCM_Arr[i] .* Σ_i[i]
    end
    if isposdef(Vi)
        begin
            cVi = cholesky!(Vi)
        end
        Vi .= inv(cVi)
    else
        Vi .= inv(Vi)
    end
end
function get_cis_snp_info(_gene_annot::DataFrame, _snp_annot::DataFrame, gene::AbstractString, window::Int; nsnp_only::Bool=false)
    geneinfo = _gene_annot[_gene_annot.pheno_id.==gene, :]
    tss_pos = geneinfo.start[1]
    end_pos = geneinfo.end[1]
    window_start_pos = max(0, tss_pos - window)
    window_end_pos = end_pos + window - 1
    if !nsnp_only
        cissnps_annot = _snp_annot[(_snp_annot.position.>=window_start_pos).&(_snp_annot.position.<=window_end_pos), :]
        cissnps_annot.pheno_id .= gene
        cissnps_annot.start_distance .= cissnps_annot.position .- tss_pos
        cissnps_annot.end_distance .= cissnps_annot.position .- end_pos .+ 1
        n_snps = size(cissnps_annot, 1)
        return cissnps_annot, n_snps
    else
        return sum((_snp_annot.position .>= window_start_pos) .& (_snp_annot.position .<= window_end_pos))
    end
end
function calculate_test_stats(Vi::Matrix{T}, xs0::Matrix{T}, yadj::Vector{T}; byrow::Bool=false) where {T}
    n_snps = size(xs0, ifelse(byrow, 2, 1))
    Vi_yc = Vi * yadj
    df_test = zeros(Float64, n_snps, 3)
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        if byrow
            xs = xs0[:, xi]
            St_Vi_S = dot(xs, mat_mul(xs, Vi)) 
        else
            xs = xs0[xi, :]
            St_Vi_S = dot(xs, mat_mul(Vi, xs)) 
        end
        var_s = 1 / St_Vi_S
        beta_s = dot(xs, Vi_yc) * var_s
        se_beta = sqrt(abs(var_s))
        zs = beta_s / se_beta
        pv_s = ccdf(Chisq(1), Float64(zs^2))
        df_test[xi, :] = [beta_s, se_beta, pv_s]
    end
    return df_test
end
function calculate_test_stats!(df_stats::Matrix{T}, Vi::Matrix{T}, xs0::Matrix{T}, yadj::Vector{T}) where {T}
    n_snps = size(xs0, 1)
    Vi_yc = mat_mul(Vi, yadj)
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = xs0[xi, :]
        St_Vi_S = dot(xs, mat_mul(Vi, xs)) 
        var_s = 1 / St_Vi_S
        beta_s = dot(xs, Vi_yc) * var_s
        se_beta = sqrt(abs(var_s))
        zs = beta_s / se_beta
        pv_s = ccdf(Chisq(1), Float64(zs^2))
        df_stats[xi, :] .= [beta_s, se_beta, pv_s]
    end
end
function calcu_one_assoc_test!(df_stats::Matrix{T1}, Vi::AbstractMatrix{T}, xs0::Matrix{T}, yadj::Vector{T}) where {T1<:AbstractFloat,T<:AbstractFloat}
    Vi_yc = mat_mul(Vi, yadj)
    Vi_s = similar(xs0)
    mat_mul!(Vi_s, xs0, Vi)
    Vi_s .*= xs0
    σ = 1.0 ./ sum(Vi_s, dims=2)
    df_stats[:, 1] .= mat_mul(xs0, Vi_yc) .* σ
    df_stats[:, 2] .= sqrt.(abs.(σ))
    df_stats[:, 3] .= ccdf(Chisq(1), abs2.(df_stats[:, 1] ./ df_stats[:, 2]))
end
function calcu_one_assoc_test!(df_stats::Matrix{T1}, Vi::AbstractMatrix{T}, xs0::Matrix{Int8}, norm_factors::Vector{T}, yadj::Vector{T}, block_size::Int64=10000) where {T1<:AbstractFloat,T<:AbstractFloat}
    n_snps, n_samples = size(xs0)
    Vi_yc = mat_mul(Vi, yadj)
    iter_collects = collect(Iterators.partition(1:n_snps, block_size))
    S0 = zeros(Float64, block_size, n_samples)
    Vi_S0 = similar(S0)
    σ0 = zeros(block_size)
    S0_Vi_yc = zeros(block_size)
    block_size = length(iter_collects[end])
    S1 = zeros(Float64, block_size, n_samples)
    Vi_S1 = similar(S1)
    σ1 = zeros(block_size)
    S1_Vi_yc = zeros(block_size)
    @inbounds @views for i in 1:(length(iter_collects)-1)
        S0 .= xs0[iter_collects[i], :]
        S0 .-= norm_factors[iter_collects[i]]
        mat_mul!(Vi_S0, S0, Vi)
        Vi_S0 .*= S0
        σ0 .= sum(Vi_S0, dims=2)
        σ0 .= 1.0 ./ σ0
        mat_mul!(S0_Vi_yc, S0, Vi_yc)
        df_stats[iter_collects[i], 1] .= S0_Vi_yc .* σ0
        df_stats[iter_collects[i], 2] .= sqrt.(abs.(σ0))
    end
    @inbounds @views begin
        S1 .= xs0[iter_collects[end], :]
        S1 .-= norm_factors[iter_collects[end]]
        mat_mul!(Vi_S1, S1, Vi)
        Vi_S1 .*= S1
        σ1 .= sum(Vi_S1, dims=2)
        σ1 .= 1.0 ./ σ1
        mat_mul!(S1_Vi_yc, S1, Vi_yc)
        df_stats[iter_collects[end], 1] .= S1_Vi_yc .* σ1
        df_stats[iter_collects[end], 2] .= sqrt.(abs.(σ1))
        df_stats[:, 3] .= ccdf(Chisq(1), abs2.(df_stats[:, 1] ./ df_stats[:, 2]))
    end
end
function calcu_one_assoc_test!(df_stats::Matrix{T}, Vi::AbstractMatrix{T}, xs0::Matrix{Int8}, norm_factors::AbstractVecOrMat{T}, yadj::Vector{T}, Vi_yc::AbstractVecOrMat{T}, iter_collects::Vector{UnitRange{Int64}}, S0::Matrix{T}, Vi_S0::Matrix{T}, σ0::Vector{T}, S0_Vi_yc::AbstractVecOrMat{T}, S1::Matrix{T}, Vi_S1::Matrix{T}, σ1::Vector{T}, S1_Vi_yc::AbstractVecOrMat{T}, byrow::Bool=false) where {T}
    if byrow
        mat_mul!(Vi_yc, yadj', Vi)
        @inbounds @views for i in 1:(length(iter_collects)-1)
            S0 .= xs0[:, iter_collects[i]]
            S0 .-= norm_factors[:, iter_collects[i]]
            mat_mul!(Vi_S0, Vi, S0)
            Vi_S0 .*= S0
            σ0 .= sum(Vi_S0, dims=1)[:]
            σ0 .\= 1
            mat_mul!(S0_Vi_yc, Vi_yc, S0)
            df_stats[iter_collects[i], 1] .= S0_Vi_yc' .* σ0
            df_stats[iter_collects[i], 2] .= sqrt.(abs.(σ0))
        end
        @inbounds @views begin
            S1 .= xs0[:, iter_collects[end]]
            S1 .-= norm_factors[:, iter_collects[end]]
            mat_mul!(Vi_S1, Vi, S1)
            Vi_S1 .*= S1
            σ1 .= sum(Vi_S1, dims=1)[:]
            σ1 .\= 1
            mat_mul!(S1_Vi_yc, Vi_yc, S1)
            df_stats[iter_collects[end], 1] .= S1_Vi_yc' .* σ1
            df_stats[iter_collects[end], 2] .= sqrt.(abs.(σ1))
            df_stats[:, 3] .= ccdf(Chisq(1), abs2.(df_stats[:, 1] ./ df_stats[:, 2]))
        end
    else
        mat_mul!(Vi_yc, Vi, yadj)
        @inbounds @views for i in 1:(length(iter_collects)-1)
            S0 .= xs0[iter_collects[i], :]
            S0 .-= norm_factors[iter_collects[i]]
            mat_mul!(Vi_S0, S0, Vi)
            Vi_S0 .*= S0
            σ0 .= sum(Vi_S0, dims=2)
            σ0 .\= 1
            mat_mul!(S0_Vi_yc, S0, Vi_yc)
            df_stats[iter_collects[i], 1] .= S0_Vi_yc .* σ0
            df_stats[iter_collects[i], 2] .= sqrt.(abs.(σ0))
        end
        @inbounds @views begin
            S1 .= xs0[iter_collects[end], :]
            S1 .-= norm_factors[iter_collects[end]]
            mat_mul!(Vi_S1, S1, Vi)
            Vi_S1 .*= S1
            σ1 .= sum(Vi_S1, dims=2)
            σ1 .\= 1
            mat_mul!(S1_Vi_yc, S1, Vi_yc)
            df_stats[iter_collects[end], 1] .= S1_Vi_yc .* σ1
            df_stats[iter_collects[end], 2] .= sqrt.(abs.(σ1))
            df_stats[:, 3] .= ccdf(Chisq(1), abs2.(df_stats[:, 1] ./ df_stats[:, 2]))
        end
    end
end
function calcu_one_assoc_test!(df_stats::Matrix{T}, Vi::AbstractMatrix{T}, xs0::Matrix{T}, yadj::Vector{T}, Vi_s::Matrix{T}) where {T}
    Vi_yc = mat_mul(Vi, yadj)
    mat_mul!(Vi_s, xs0, Vi)
    Vi_s .*= xs0
    σ = 1.0 ./ sum(Vi_s, dims=2)
    df_stats[:, 1] .= mat_mul(xs0, Vi_yc) .* σ
    df_stats[:, 2] .= sqrt.(abs.(σ))
    df_stats[:, 3] .= ccdf(Chisq(1), abs2.(df_stats[:, 1] ./ df_stats[:, 2]))
end
function calcu_two_assoc_test!(df_stats::Matrix{T1}, Vi::Matrix{T}, S0::Matrix{T}, S1::Matrix{T}, yadj::Vector{T}, no_perm::Bool) where {T1<:AbstractFloat,T<:AbstractFloat}
    n_snps, n_samples = size(S0)
    Vi_yc = Vi * yadj
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = zeros(Float64, n_samples, 2)
        xs[:, 1] .= S0[xi, :]
        xs[:, 2] .= S1[xi, :]
        St_Vi_S = xs' * mat_mul(Vi, xs) 
        var_s = diag(inv(St_Vi_S))
        β = xs' * Vi_yc .* var_s
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        pv_each = ccdf(Chisq(1), abs2.(zs))
        WaldStat = β' * St_Vi_S * β
        pv_comb = ccdf.(Chisq(2), WaldStat)
        df_stats[xi, :] .= [β[1], β_se[1], pv_each[1], β[2], β_se[2], pv_each[2], pv_comb]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, Vi::Matrix{T}, S0::Matrix{T}, S1::Matrix{T}, yadj::Vector{T}) where {T}
    n_snps, n_samples = size(S0)
    Vi_yc = Vi * yadj
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = zeros(Float64, n_samples, 2)
        xs[:, 1] .= S0[xi, :]
        xs[:, 2] .= S1[xi, :]
        St_Vi_S = mat_mul(mat_mul(xs', Vi), xs) 
        var_s = diag(inv(St_Vi_S))
        β = xs' * Vi_yc .* var_s
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        pv_each = ccdf(Chisq(1), abs2.(zs))
        WaldStat = β' * St_Vi_S * β
        pv_comb = ccdf.(Chisq(2), WaldStat)
        df_stats[xi, :] .= [β[1], β_se[1], pv_each[1], β[2], β_se[2], pv_each[2], pv_comb]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, Vi::Matrix{T}, S0::Matrix{T}, iterm::Vector{T}, yadj::Vector{T}; byrow::Bool=false) where {T}
    n_snps = size(S0, ifelse(byrow, 2, 1))
    Vi_yc = Vi * yadj
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = hcat(S0[xi, :], S0[xi, :] .* iterm)
        St_Vi_S = xs' * mat_mul(Vi, xs) 
        var_s = diag(inv(St_Vi_S))
        β = xs' * Vi_yc .* var_s
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        pv_each = ccdf(Chisq(1), abs2.(zs))
        WaldStat = β' * St_Vi_S * β
        pv_comb = ccdf.(Chisq(2), WaldStat)
        df_stats[xi, :] .= [β[1], β_se[1], pv_each[1], β[2], β_se[2], pv_each[2], pv_comb]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, vp::Float64, S0::Matrix{T}, S1::Matrix{T}, yadj::Vector{T}; byrow::Bool=false) where {T}
    n_snps, n_samples = size(S0)
    Vi_yc = vp * yadj
    @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = zeros(Float64, n_samples, 2)
        @views xs[:, 1] .= S0[xi, :]
        @views xs[:, 2] .= S1[xi, :]
        St_Vi_S = xs' * vp * xs 
        var_s = diag(inv(St_Vi_S))
        β = xs' * Vi_yc .* var_s
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        pv_each = ccdf(Chisq(1), abs2.(zs))
        WaldStat = β' * St_Vi_S * β
        pv_comb = ccdf.(Chisq(2), WaldStat)
        df_stats[xi, :] .= [β[1], β_se[1], pv_each[1], β[2], β_se[2], pv_each[2], pv_comb]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, vp::T, S0::Matrix{T}, S1::Matrix{T}, yadj::Vector{T}, xs_p_1_1::Vector{T}, xs_p_2_2::Vector{T}, xs_p_2_1::Vector{T}; byrow::Bool=false) where {T}
    n_snps, n_samples = size(S0)[ifelse(byrow, [2, 1], [1, 2])]
    Vi_yc = vp * yadj
    S0_Vi_yc = zeros(n_snps)
    S1_Vi_yc = zeros(n_snps)
    if byrow
        mul!(S0_Vi_yc, S0', Vi_yc) 
        mul!(S1_Vi_yc, S1', Vi_yc) 
    else
        mul!(S0_Vi_yc, S0, Vi_yc) 
        mul!(S1_Vi_yc, S1, Vi_yc) 
    end
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        St_Vi_S = SMatrix{2,2}(xs_p_1_1[xi] * vp, xs_p_2_1[xi] * vp, xs_p_2_1[xi] * vp, xs_p_2_2[xi] * vp) 
        var_s = diag(inv(St_Vi_S))
        β = MVector{2}(0.0, 0.0)
        β[1] = S0_Vi_yc[xi] * var_s[1]
        β[2] = S1_Vi_yc[xi] * var_s[2]
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        df_stats[xi, :] .= [zs[1], zs[2], NaN]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, vp::T, S0::Matrix{T}, S1::Vector{T}, yadj::Vector{T}, xs_p_1_1::Vector{T}, xs_p_2_2::Vector{T}, xs_p_2_1::Vector{T}; byrow::Bool=false) where {T}
    n_snps, n_samples = size(S0)[ifelse(byrow, [2, 1], [1, 2])]
    Vi_yc = vp * yadj
    S0_Vi_yc = zeros(n_snps)
    S1_Vi_yc = zeros(n_snps)
    if byrow
        mul!(S0_Vi_yc, S0', Vi_yc) 
        mul!(S1_Vi_yc, S0', Vi_yc .* S1) 
    else
        mul!(S0_Vi_yc, S0, Vi_yc) 
        mul!(S1_Vi_yc, S0', Vi_yc .* S1) 
    end
    @views @inbounds Threads.@threads for xi in range(1, n_snps)
        St_Vi_S = SMatrix{2,2}(xs_p_1_1[xi] * vp, xs_p_2_1[xi] * vp, xs_p_2_1[xi] * vp, xs_p_2_2[xi] * vp) 
        var_s = diag(inv(St_Vi_S))
        β = MVector{2}(0.0, 0.0)
        β[1] = S0_Vi_yc[xi] * var_s[1]
        β[2] = S1_Vi_yc[xi] * var_s[2]
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        df_stats[xi, :] .= [zs[1], zs[2], NaN]
    end
    return df_stats
end
function calcu_two_assoc_test!(df_stats::Matrix{T}, vp::Float64, S0::Matrix{T}, iterm::Vector{T}, yadj::Vector{T}; byrow::Bool=false) where {T}
    n_snps = size(S0, ifelse(byrow, 2, 1))
    Vi_yc = vp * yadj
    @inbounds Threads.@threads for xi in range(1, n_snps)
        xs = hcat(S0[xi, :], S0[xi, :] .* iterm)
        St_Vi_S = xs' * vp * xs 
        var_s = diag(inv(St_Vi_S))
        β = xs' * Vi_yc .* var_s
        β_se = sqrt.(abs.(var_s))
        zs = β ./ β_se
        pv_each = ccdf(Chisq(1), abs2.(zs))
        WaldStat = β' * St_Vi_S * β
        pv_comb = ccdf.(Chisq(2), WaldStat)
        df_stats[xi, :] .= [β[1], β_se[1], pv_each[1], β[2], β_se[2], pv_each[2], pv_comb]
    end
    return df_stats
end
function calculate_perm_stats(yadj::Vector{T}, xs0T_DIVIDE_xs0T_xs0::Matrix{T}, xs0T_xs0_inv::Vector{T}) where {T}
    n_snps = length(xs0T_xs0_inv)
    varY = var(yadj)
    beta_s_xs0 = xs0T_DIVIDE_xs0T_xs0 * yadj
    var_s_xs0 = xs0T_xs0_inv * varY
    se_beta_xs0 = sqrt.(abs.(var_s_xs0))
    zs_xs0 = beta_s_xs0 ./ se_beta_xs0
    pv_s_xs0_min = ccdf(Chisq(1), maximum(abs2, zs_xs0))
    return pv_s_xs0_min
end
function fast_linear_permutation_test(genotype_resid::Matrix{T}, y_perms_res::Matrix{T}, dof::Int64) where {T}
    r_nominal = zeros(size(genotype_resid, 1), size(y_perms_res, 2))
    mat_mul!(r_nominal, genotype_resid, y_perms_res)
    r2_nominal_max = maximum(abs2, r_nominal, dims=1)
    tstat2 = dof .* r2_nominal_max ./ (1 .- r2_nominal_max)
    df_pval_perm = vec(ccdf(Chisq(1), tstat2))
    return df_pval_perm
end
function fast_linear_permutation_test!(genotype_resid::Matrix{T}, y_perms_res::Matrix{T}, max_absr_perms::Vector{T}; byrow::Bool=false) where {T}
    if byrow
        r_perms = mat_mul(genotype_resid', y_perms_res)
        max_absr_perms .= vec(maximum(abs, r_perms, dims=1))
    else
        r_perms = mat_mul(genotype_resid, y_perms_res)
        max_absr_perms .= vec(maximum(abs, r_perms, dims=1))
    end
end
function fast_linear_permutation_test!(genotype_resid::Matrix{T}, y_perms_res::Matrix{T}, max_absr_perms::Vector{T}, variant_subset_index::Union{Nothing,AbstractVector{Int}}) where {T}
    r_perms = mat_mul(y_perms_res', genotype_resid[:,variant_subset_index])
    return r_perms
end
function fast_linear_permutation_test(genotype_resid::Matrix{T}, y_perms_res::Matrix{T}) where {T}
    r_perms = mat_mul(genotype_resid', y_perms_res)
    return findmax(abs, r_perms, dims=1)
end
function fast_linear_permutation_test!(genotype_resid::Matrix{T}, y_perms_res::Matrix{T}, max_absr_perms::Vector{T}, nbatchs::Int; byrow::Bool=false) where {T}
    _ns = size(y_perms_res, 2)
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        if byrow
            r_perms = zeros(size(genotype_resid, 2), size(y_perms_res, 2))
            mat_mul!(r_perms, genotype_resid', y_perms_res)
            max_absr_perms .= vec(maximum(abs, r_perms, dims=1))
        else
            r_perms = zeros(size(genotype_resid, 1), size(y_perms_res, 2))
            mat_mul!(r_perms, genotype_resid, y_perms_res)
            max_absr_perms .= vec(maximum(abs, r_perms, dims=1))
        end
    end
end
function fast_permutation_clipper(genotype_resid::Matrix{T}, phenotype_res_t::Vector{T}, n_perms::Int; byrow::Bool=false) where {T}
    df_absr_perm_max = zeros(n_perms)
    @inbounds Threads.@threads for perm_i in range(1, n_perms)
        Random.seed!(_args_seed + perm_i)
        y_perm = shuffle(phenotype_res_t)
        if byrow
            r_perm = genotype_resid' * y_perm
        else
            r_perm = genotype_resid * y_perm
        end
        df_absr_perm_max[perm_i] = maximum(abs, r_perm)
    end
    return df_absr_perm_max
end
function fast_permutation_clipper!(max_absr_perms::Vector{T}, genotype_resid::Matrix{T}, phenotype_res_t::Vector{T}, n_perms::Int; byrow::Bool=false) where {T}
    @inbounds Threads.@threads for perm_i in range(1, n_perms)
        Random.seed!(_args_seed + perm_i)
        y_perm = shuffle(phenotype_res_t)
        if byrow
            r_perm = genotype_resid' * y_perm
        else
            r_perm = genotype_resid * y_perm
        end
        max_absr_perms[perm_i] = maximum(abs, r_perm)
    end
end
function fast_permutation_clipper!(r_perm_mat::AbstractArray{T}, genotype_resid::Matrix{T}, phenotype_res_t::Vector{T}, n_perms::Int, variant_subset_index::Union{Nothing,AbstractVector{Int}}) where {T}
    @inbounds Threads.@threads for perm_i in range(1, n_perms)
        Random.seed!(_args_seed + perm_i)
        y_perm = shuffle(phenotype_res_t)
        r_perm_mat[perm_i, :] .= vec(y_perm' * genotype_resid[:, variant_subset_index])
    end
end
function calculate_perm(perm_i, exppheno, X_MME, XtX, xs0T_DIVIDE_xs0T_xs0, xs0T_xs0_inv)
    Random.seed!(_args_seed + perm_i)
    y_perm = shuffle(exppheno)
    yadj = get_lm_residuals(X_MME, y_perm, XtX)
    pv_perm_min = calculate_perm_stats(yadj, xs0T_DIVIDE_xs0T_xs0, xs0T_xs0_inv)
    return pv_perm_min
end
function get_approx_slope_se_from_r!(df_test::AbstractArray{T1}, assoc_r::AbstractVector{T}, dof::Int, std_g::AbstractVector{T}, std_y::T; is_calcu_pv::Bool=true) where {T1<:AbstractFloat,T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    slope = assoc_r .* std_y ./ std_g
    slope_se = abs.(slope) ./ tstat
    if is_calcu_pv
        pv_s = cdf(TDist(dof), -tstat) * 2
        df_test .= [slope slope_se pv_s]
    else
        df_test .= [slope slope_se tstat]
    end
end
function get_approx_slope_se_from_r(assoc_r::Vector{T}, dof::Int, std_g::Vector{T}, std_y::Union{T,Vector{T}}; is_calcu_pv::Bool=true) where {T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    slope = assoc_r .* std_y ./ std_g
    slope_se = abs.(slope) ./ tstat
    if is_calcu_pv
        pv_s = cdf(TDist(dof), -tstat) * 2
        df_test = [slope slope_se pv_s]
    else
        df_test = [slope slope_se tstat]
    end
    return df_test
end
function get_approx_p_from_r(assoc_r::Vector{T}, dof::Int) where {T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    clamp!(assoc_r2, 0, 1)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    pv_s = cdf(TDist(dof), -tstat) * 2
    return pv_s
end
function get_approx_p_from_r(assoc_r::Vector{T}, dof::T) where {T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    clamp!(assoc_r2, 0, 1)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    pv_s = cdf(TDist(dof), -tstat) * 2
    return pv_s
end
function get_approx_p_from_r(assoc_r::Matrix{T}, dof::Int) where {T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    clamp!(assoc_r2, 0, 1)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    pv_s = cdf(TDist(dof), -tstat) * 2
    return pv_s
end
function get_approx_p_from_r(assoc_r::Matrix{T}, dof::T) where {T<:AbstractFloat}
    assoc_r2 = abs2.(assoc_r)
    clamp!(assoc_r2, 0, 1)
    tstat = abs.(assoc_r .* sqrt.(dof ./ (1 .- assoc_r2)))
    pv_s = cdf(TDist(dof), -tstat) * 2
    return pv_s
end
function get_approx_absr(slope::T1, slope_se::T1, dof::T2) where {T1<:AbstractFloat,T2<:Real}
    tstat2 = abs2(slope / slope_se)
    r2 = 1 / (dof / tstat2 + 1)
    absr = sqrt(r2)
    return absr
end
function get_approx_absr(zstat::T1, dof::T2; square::Bool=true) where {T1<:AbstractFloat,T2<:Real}
    if square
        tstat2 = abs2(zstat)
    else
        tstat2 = zstat
    end
    r2 = 1 / (dof / tstat2 + 1)
    absr = sqrt(r2)
    return absr
end
function get_lm_residuals(X_MME::Matrix{T}, y_perm::AbstractVecOrMat{T}) where {T}
    beta0 = (X_MME' * X_MME) \ (X_MME' * y_perm)
    yresid = y_perm - X_MME * beta0
    return yresid
end
function get_lm_residuals(X_MME::Matrix{T}, y_perm::AbstractVecOrMat{T}, XtX::Matrix{T}) where {T}
    beta0 = XtX \ (X_MME' * y_perm)
    yresid = y_perm - X_MME * beta0
    return yresid
end
function get_lm_residuals!(X_MME::Matrix{T}, y_perm::AbstractVecOrMat{T}, XtX::Matrix{T}) where {T}
    beta0 = XtX \ (X_MME' * y_perm)
    y_perm .-= X_MME * beta0
end
function get_lm_residuals!(SQ_adj::Matrix{T}, X_MME::Matrix{T}, SQ::Matrix{T}, SQ2::Matrix{T}) where {T}
    _n, _X_c = size(X_MME)
    X_s2 = zeros(_n, _X_c + 1)
    X_s2[:, 1:_X_c] .= X_MME
    X_s2_t_X_s2 = zeros(_X_c + 1, _X_c + 1)
    beta0 = zeros(_X_c + 1)
    @inbounds @views for xi in 1:size(SQ2, 2)
        X_s2[:, end] .= SQ2[:, xi]
        mat_mul!(X_s2_t_X_s2, X_s2', X_s2)
        mat_mul!(beta0, X_s2', SQ[:, xi])
        beta0 .= X_s2_t_X_s2 \ beta0
        SQ_adj[:, xi] .= SQ[:, xi] - X_s2 * beta0
    end
end
function get_t_pval(t, df)
    """
    Get p-value corresponding to t statistic and degrees of freedom (df). t and/or df can be arrays.
    If log=True, returns -log10(P).
    """
    return ccdf(TDist(df), abs.(t)) * 2
end
function pval_from_corr(r2::Union{Vector{T},T}, dof) where {T<:AbstractFloat}
    tstat2 = dof * r2 ./ (1 .- r2)
    return get_t_pval(sqrt.(tstat2), dof)
end
function df_cost(r2::Union{Vector{T},T}, dof) where {T<:AbstractFloat}
    """minimize abs(1-alpha) as a function of M_eff"""
    pval = pval_from_corr(r2, dof)
    m = Statistics.mean(pval)
    v = var(pval)
    return m * (m * (1 - m) / v - 1) - 1
end
function beta_log_likelihood(x::Vector{T}, shape1::T, shape2::T) where {T<:AbstractFloat}
    """negative log-likelihood of beta distribution"""
    logbeta = loggamma(shape1) + loggamma(shape2) - loggamma(shape1 + shape2)
    return (1 - shape1) * sum(log.(x)) + (1 - shape2) * sum(log.(1 .- x)) + length(x) * logbeta
end
function fit_beta_parameters(r2_perm::Vector{T}, dof_init, tol=1e-4, return_minp::Bool=false) where {T<:AbstractFloat}
    """
    r2_perm:    array of max. r2 values from permutations
    dof_init:   degrees of freedom
    """
    true_dof = dof_init
    Optim_converged = false
    try
        res = optimize(x -> abs(df_cost(r2_perm, x[1])), [Float64(dof_init)], Newton(), Optim.Options(x_tol=tol, iterations=50)) 
        Optim_converged = Optim.converged(res)
        true_dof = Optim.minimizer(res)[1]  
    catch e
        println(e)
        println("WARNING: Newton's method failed to converge (running Nelder-Mead)")
        res = optimize(x -> abs(df_cost(r2_perm, x[1])), [FloatT(dof_init)], NelderMead(), Optim.Options(x_tol=tol))
        true_dof = Optim.minimizer(res)[1]  
        Optim_converged = true
    end
    if !Optim_converged
        println("WARNING: Newton's method failed to converge (running Nelder-Mead)")
        res = optimize(x -> abs(df_cost(r2_perm, x[1])), [FloatT(dof_init)], NelderMead(), Optim.Options(x_tol=tol))
        true_dof = Optim.minimizer(res)[1]
    end
    pval = pval_from_corr(r2_perm, true_dof)
    if true
        mean_pval = Statistics.mean(pval)
        var_pval = var(pval)
        beta_shape1 = mean_pval * (mean_pval * (1 - mean_pval) / var_pval - 1)
        beta_shape2 = beta_shape1 * (1 / mean_pval - 1)
        res = optimize(s -> beta_log_likelihood(pval, s[1], s[2]), [beta_shape1, beta_shape2], NelderMead(), Optim.Options(x_tol=tol))
        beta_shape1, beta_shape2 = res.minimizer
    else
        fit_beta = fit(Beta, pval)
        beta_shape1, beta_shape2 = fit_beta.α, fit_beta.β
    end
    if return_minp
        return beta_shape1, beta_shape2, true_dof, pval
    else
        return beta_shape1, beta_shape2, true_dof
    end
end
function calculate_beta_approx_pval(r2_perm, r2_nominal, dof_init, tol=1e-4)
    """
      r2_nominal: nominal max. r2 (scalar or array)
      r2_perm:    array of max. r2 values from permutations
      dof_init:   degrees of freedom
    """
    beta_shape1, beta_shape2, true_dof = fit_beta_parameters(r2_perm, dof_init, tol)
    pval_true_dof = pval_from_corr(r2_nominal, true_dof)
    pval_beta = cdf(Beta(beta_shape1, beta_shape2), pval_true_dof)
    return [pval_beta, beta_shape1, beta_shape2, true_dof, pval_true_dof]
end
