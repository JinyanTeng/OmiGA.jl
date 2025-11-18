function cis_permutation_one_assoc!(permutation_method::String, gene_i::Int, exppheno::Vector{T}, Q::Union{Nothing,Matrix{T}}, cis_genotype::Matrix{T}, absr_perm::Vector{T}, _df_perm::DataFrame, _df_tops::DataFrame, _df_full::DataFrame, dof::Int; y_perms::Union{Nothing,Matrix{T}}=nothing, y_perms_res::Union{Nothing,Matrix{T}}=nothing, is_sample_byrow::Bool=true, fdr::T=T(0.05), variant_subset_index::Union{Nothing,AbstractVector{Int}}=nothing) where {T<:AbstractFloat}
    n_cis_snps = size(cis_genotype, 2)
    n_perms = length(absr_perm)
    @timeit to "Perm 1.1" if permutation_method == "standard_fast"
        y_perms .= exppheno
        Random.seed!(_args_seed)
        @inbounds @views for perm_i in 1:n_perms
            shuffle!(y_perms[:, perm_i])
        end
        get_matrix_resid!(y_perms, Q, return_std=false)
        if !isnothing(variant_subset_index)
           return fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm, variant_subset_index) 
        end
        fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm; byrow=is_sample_byrow) 
    elseif permutation_method == "clipper"
        if !isnothing(variant_subset_index)
            return fast_permutation_clipper!(absr_perm, cis_genotype, exppheno, n_perms, variant_subset_index)
         end
        fast_permutation_clipper!(absr_perm, cis_genotype, exppheno, n_perms; byrow=is_sample_byrow)
    end
    pval_g1_threshold = sort(get_approx_p_from_r(absr_perm, dof))[ceil(Int, n_perms * fdr)]
    @timeit to "Perm 1.2" begin
        index_top = first_nonnan_argmin(_df_full.pval_g1) 
        if index_top == -1
            return false
        else
            absr_exper = get_approx_absr(_df_full.beta_g1[index_top], _df_full.beta_se_g1[index_top], dof) 
            _df_perm[gene_i, 3:end] .= [absr_exper; absr_perm]
            _df_tops.num_var[gene_i] = n_cis_snps
            _df_tops.variant_id[gene_i] = _df_full.variant_id[index_top]
            _df_tops.start_distance[gene_i] = _df_full.start_distance[index_top]
            _df_tops.af[gene_i] = _df_full.af[index_top]
            _df_tops.beta_g1[gene_i] = _df_full.beta_g1[index_top]
            _df_tops.beta_se_g1[gene_i] = _df_full.beta_se_g1[index_top]
            _df_tops.pval_g1[gene_i] = _df_full.pval_g1[index_top]
            _df_tops.pval_g1_threshold[gene_i] = pval_g1_threshold
        end
    end
    return true
end
function cis_permutation_two_assoc!(permutation_method::String, gene_i::Int, exppheno::Vector{T}, cis_genotype::Matrix{T}, cis_genotype_iterm::Matrix{T},
    _df_perm::DataFrame, _df_tops::DataFrame, _df_full::DataFrame, y_perms::Union{Nothing,Matrix{T}}, cissnps_annot::DataFrame, stat_perm::Vector{T};
    X_MME::Union{Nothing,Matrix{T}}=nothing, Xt_X::Union{Nothing,Matrix{T}}=nothing, glo_EA::Union{Nothing,Eigen{T,T,Matrix{T},Vector{T}}}=nothing,
    chrom_xtx_SQt_xQ::Union{Nothing,Matrix{T}}=nothing, chrom_xtx_SQ2t_xQ::Union{Nothing,Matrix{T}}=nothing, xW0::Union{Nothing,Matrix{T}}=nothing, xty_init::Union{Nothing,Matrix{T}}=nothing, xtx_init::Union{Nothing,Matrix{T}}=nothing, fdr::T=T(0.05), variant_subset_index::Union{Nothing,AbstractVector{Int}}=nothing) where {T<:AbstractFloat}
    n_cis_snps = size(cissnps_annot, 1)
    n_perms = size(y_perms, 2)
    dof = size(X_MME, 1) - 2 - size(X_MME, 2)
    if permutation_method == "standard"
        _c2 = size(X_MME, 2) + 2
        @timeit to "Perm 1.1" begin
            y_perms .= exppheno 
            Random.seed!(_args_seed)
            @views @inbounds for perm_i in 1:n_perms
                shuffle!(y_perms[:, perm_i])
            end
            if _args_preadj_covar
                get_lm_residuals!(X_MME, y_perms, Xt_X)
            end
            y_perms .= glo_EA.vectors' * y_perms
            @timeit to "i0.0" begin 
                xtx_SQt_xQ = zeros(T, _c2, n_cis_snps)
                xtx_SQ2t_xQ = similar(xtx_SQt_xQ)
                Ainds = CartesianIndices(xtx_SQt_xQ)
                Binds = CartesianIndices((1:_c2, minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                copyto!(xtx_SQt_xQ, Ainds, chrom_xtx_SQt_xQ, Binds)
                copyto!(xtx_SQ2t_xQ, Ainds, chrom_xtx_SQ2t_xQ, Binds)
                mat_mul!(xty_init, xW0', y_perms)
            end
            @timeit to "i0.1" df_test_perm = zeros(T, n_perms, n_cis_snps * 2)
        end
        @timeit to "Perm 1.2" begin
            if _args_run_mode == "cis_interaction"
                idul_two_assoc_test_perm!(df_test_perm, cis_genotype, y_perms, cis_genotype_iterm, xW0, xtx_init, xty_init, xtx_SQt_xQ, xtx_SQ2t_xQ)
                max_absz_g1 = maximum(abs, df_test_perm[:, 1:2:size(df_test_perm, 2)], dims=2)
                max_absz_g2 = maximum(abs, df_test_perm[:, 2:2:size(df_test_perm, 2)], dims=2)
                stat_perm[1:3:length(stat_perm)] .= get_approx_absr.(max_absz_g1, dof)
                stat_perm[2:3:length(stat_perm)] .= get_approx_absr.(max_absz_g2, dof)
            elseif _args_run_mode == "cis"
                error("!")
            end
        end
    elseif permutation_method == "standard_fast"  
        @timeit to "Perm 1.1" begin
            y_perms .= exppheno
            Random.seed!(_args_seed)
            @views @inbounds for perm_i in 1:n_perms
                shuffle!(y_perms[:, perm_i])
            end
            get_lm_residuals!(X_MME, y_perms, Xt_X) 
            y_perms ./= sqrt.(sum(abs2, y_perms, dims=1))
            y_perms .-= Statistics.mean(y_perms, dims=1)
        end
        @timeit to "Perm 1.2" begin
            if _args_run_mode == "cis_interaction"
                if true
                    if !isnothing(variant_subset_index) 
                        perm_r_mat_g2 = y_perms' * cis_genotype_iterm[:, variant_subset_index]
                        return perm_r_mat_g2
                    end
                    perm_r_mat_g2 = y_perms' * cis_genotype_iterm
                    max_absr_g2 = vec(maximum(abs, perm_r_mat_g2, dims=2))
                    stat_perm[2:3:length(stat_perm)] .= max_absr_g2
                end
            elseif _args_run_mode == "cis"
                if !isnothing(variant_subset_index) 
                    perm_r_mat_g2 = y_perms' * cis_genotype_iterm[:, variant_subset_index]
                    return perm_r_mat_g2
                end
                perm_r_mat_g2 = y_perms' * cis_genotype_iterm
                max_absr_g2 = maximum(abs, perm_r_mat_g2, dims=2)
                stat_perm[2:3:length(stat_perm)] .= max_absr_g2
            end
        end
    end
    pval_g2_threshold = sort(get_approx_p_from_r(max_absr_g2, dof - 1))[ceil(Int, n_perms * fdr)] 
    @timeit to "Perm 1.2" begin
        index_top_g1 = first_nonnan_argmin(_df_full.pval_g1) 
        index_top_g2 = first_nonnan_argmin(_df_full.pval_g2) 
        sl1 = _df_full.beta_g1[index_top_g1]
        sl1_se = _df_full.beta_se_g1[index_top_g1]
        sl2 = _df_full.beta_g2[index_top_g2]
        sl2_se = _df_full.beta_se_g2[index_top_g2]
        r_g1 = get_approx_absr(sl1, sl1_se, dof)
        r_g2 = get_approx_absr(sl2, sl2_se, dof)
        r_g12 = T(NaN)
        _df_perm[gene_i, 3:end] .= [r_g1; r_g2; r_g12; stat_perm]
        if _args_run_mode == "cis_interaction"
            index_top = index_top_g2
        else
            index_top = index_top_g2
        end
        _df_tops.num_var[gene_i] = n_cis_snps
        _df_tops.variant_id[gene_i] = _df_full.variant_id[index_top]
        _df_tops.start_distance[gene_i] = _df_full.start_distance[index_top]
        _df_tops.af[gene_i] = _df_full.af[index_top]
        _df_tops.beta_g1[gene_i] = _df_full.beta_g1[index_top]
        _df_tops.beta_se_g1[gene_i] = _df_full.beta_se_g1[index_top]
        _df_tops.pval_g1[gene_i] = _df_full.pval_g1[index_top]
        _df_tops.beta_g2[gene_i] = _df_full.beta_g2[index_top]
        _df_tops.beta_se_g2[gene_i] = _df_full.beta_se_g2[index_top]
        _df_tops.pval_g2[gene_i] = _df_full.pval_g2[index_top]
        _df_tops.pval_g2_threshold[gene_i] = pval_g2_threshold
    end
end
