function runOmiGA_independent(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM::Union{Nothing,Dominance}=nothing, _struct_ITERM::Union{Nothing,InteractionTerm}=nothing)
    phenotype = _struct_PHENO.phenotype
    pheno_annotation = _struct_PHENO.annotation
    X_MME = _struct_COVAR.X_MME
    if !isnothing(_struct_PRIOR_HerB)
        prior_h2 = _struct_PRIOR_HerB.h2
        prior_B = _struct_PRIOR_HerB.B
    else
        prior_h2 = prior_B = nothing
    end
    snp_annotation = _struct_GENO.annotation
    if length(_args_cis_output) == 1
        df_tops_mt = CSV.File(_args_cis_output[1], header=true, buffer_in_memory=true) |> DataFrame
        if !any(.!isnothing.(match.(r"qval", names(df_tops_mt))))
            error("qval column must be provided in cis-file!")
        end
    else
        error("Please provide correct cis-file!")
    end
    non_pcols = isnothing.(match.(r"pval|qval", names(df_tops_mt)))
    df_tops_mt[!, names(df_to_float32!(df_tops_mt[:, non_pcols]))] = df_to_float32!(df_tops_mt[:, non_pcols])
    df_tops_mt = df_tops_mt[findall(df_tops_mt.pheno_id .∈ (pheno_annotation.pheno_id,)), :]
    with_group = !isnothing(_args_pheno_group_file)
    qtl_map_algo = _args_qtl_map_algo
    qtl_map_model = nothing
    if !isfile(replace(_args_cis_output[1], r"txt.gz$" => "info"))
        req_colnames = ["map_model"]
        if !issubset(req_colnames, names(df_tops_mt))
            error(req_colnames, " are required for `--mode cis_independent`")
        end
        qtl_map_model = df_tops_mt.map_model[1]
    else
        df_info = CSV.File(replace(_args_cis_output[1], r"txt.gz$" => "info")) |> DataFrame
        qtl_map_model = df_info.map_model[1]
    end
    _n_samples = _struct_GENO.n_samples
    _n, _X_c = size(X_MME)
    if with_group & _args_est_ind_h2
        est_ind_h2 = false
        @info "--est-ind-h2 do not supported for group cis-QTL!"
    else
        est_ind_h2 = _args_est_ind_h2
    end
    if isnothing(_args_h2_algo)
        h2_algo = "idul"
    else
        h2_algo = _args_h2_algo
    end
    n_tests = 1
    n_grms = 1
    if qtl_map_model in ["a" "a+ai" "a+ai+A"]
        @error "The current version does not support the conditional independent analysis for results from linear model or interaction cis-QTL."
    end
    if qtl_map_model == "a+A"
        glo_EA = _struct_KIN.EA
        genotype = _struct_GENO.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
            glo_KS = _struct_KIN.GRM
        end
        glo_h2_model = "Ag"
    end
    if qtl_map_model == "d+D"
        glo_EA = _struct_KIN.domEA
        genotype = _struct_DOM.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
            glo_KS = _struct_KIN.domGRM
        end
        glo_h2_model = "Dg"
    end
    if qtl_map_model == "d+A"
        glo_EA = _struct_KIN.EA
        genotype = _struct_DOM.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
            glo_KS = _struct_KIN.GRM
        end
        glo_h2_model = "Ag"
    end
    if qtl_map_model == "a+d+A+D"
        @error "The current version does not support the conditional independent analysis for '--qtl-map-model " * qtl_map_model * "'."
        genotype = _struct_GENO.genotype
        domGenotype = _struct_DOM.genotype
        n_tests = 2
        n_grms = 2
        glo_h2_model = "Ag+Dg"
    end
    if qtl_map_model == "d+A+D"
        @error "The current version does not support the conditional independent analysis for '--qtl-map-model " * qtl_map_model * "'."
        genotype = _struct_DOM.genotype
        n_grms = 2
        glo_h2_model = "Ag+Dg"
    end
    calcu_variant_threshold = _args_calcu_variant_threshold
    if !isnothing(_args_multiple_testing)
        println_to_file(string("Note: the option '--multiple-testing' is disabled for `--mode cis_independent`."), log_file)
    end
    multiple_testing_method = nothing
    n_perms = _args_n_perms
    nominal_only = _args_nominal_only
    if multiple_testing_method == "acat"
        nominal_only = true
    else
        calcu_variant_threshold = true
    end
    if calcu_variant_threshold
        nominal_only = false
        if multiple_testing_method == "acat"
            permutation_method = "standard_fast"
        elseif multiple_testing_method == "beta_approx"
            permutation_method = "standard_fast"
        else
            permutation_method = nothing
        end
    else
        permutation_method = multiple_testing_method
        println_to_file(string(" * The variant-level nominal P-value threshold for each phenotype will not be computed. If you want to compute it, please specify the flag '--calcu-variant-threshold'."), log_file)
    end
    if !nominal_only
        if isnothing(permutation_method)
            if n_tests == 1
                if _n_samples < 450
                    permutation_method = "standard_fast"
                else
                    permutation_method = "clipper"
                end
            elseif n_tests == 2
                permutation_method = "standard_fast"
            else
                permutation_method = "standard"
            end
        end
        @runif isnothing(n_perms) if permutation_method in ["standard", "standard_fast"]
            n_perms = 1000
        elseif permutation_method == "clipper"
            n_perms = 100
        end
        if (n_tests == 2) | _args_preadj_covar
            Xt_X = X_MME' * X_MME
        end
        if permutation_method in ["standard", "standard_fast", "clipper"]
            if permutation_method == "standard_fast"
                y_perms = fill(NAN, _n_samples, n_perms)
                absr_perm = fill(NAN, n_perms)
                perm_phen_type = "raw"
                perm_geno_type = "residualized"
            end
            if permutation_method == "clipper"
                absr_perm = fill(NAN, n_perms)
                perm_phen_type = "residualized"
                perm_geno_type = "residualized"
            end
            if (permutation_method == "standard") & (qtl_map_algo == "idul")
                y_perms = fill(NAN, _n_samples, n_perms)
                absr_perm = fill(NAN, n_perms)
                perm_phen_type = "raw"
                perm_geno_type = "raw"
            end
        end
        if (n_tests == 2) & (permutation_method == "standard") & (qtl_map_algo == "idul")
            _xQ_c = size(xQ, 2)
            _c2 = _xQ_c + 2
            xW0 = [ones(FloatT, _n, 2) xQ]
            xtx_init = mat_mul(xW0', xW0)
            xty_init = zeros(FloatT, _c2, n_perms)
        end
        calcu_variant_threshold = true
    else
        n_perms = 0
    end
    _args_PreV = n_grms == 2
    @timeit to "Pre-calcu V" if _args_PreV
        path_partialv = _args_path_partialv
        if !isnothing(path_partialv)
            _struct_PartialV = load(path_partialv)["data"]
        else
            error("--path-pre-eigen must be specified for QTL mapping models with two GRMs!")
        end
        if isnothing(prior_h2)
            error("--her-output must be specified for QTL mapping models with two GRMs!")
        end
    end
    if ((isnothing(prior_h2) & (h2_algo == "idul")) | (qtl_map_algo == "idul")) & (n_grms == 1)
        if _args_preadj_covar
            xQ = glo_EA.vectors' * X_MME[:, 1:1]
        else
            xQ = glo_EA.vectors' * X_MME 
        end
    end
    if (qtl_map_algo == "idul") & (qtl_map_model in ["d+A+D"])
        xQ = similar(X_MME)
        yQ = zeros(FloatT, _n_samples)
        cGRM = zeros(FloatT, _n_samples, _n_samples)
        GRM64 = zeros(Float64, _n_samples, _n_samples)
        glo_EA = Eigen{FloatT,FloatT,Matrix{FloatT},Vector{FloatT}}(zeros(FloatT, _n_samples), zeros(FloatT, _n_samples, _n_samples))
        Inds_AD_vecs = CartesianIndices(glo_EA.vectors)
    end
    if qtl_map_algo == "standard"
        Vi = zeros(FloatT, _n, _n)
        Di = Diagonal(zeros(FloatT, _n, _n))
        EAvec_Di = zeros(FloatT, _n, _n)
        if qtl_map_model in ["a+A+D", "d+A+D", "a+d+A+D"]
            EAvec = zeros(FloatT, _n, _n)
        end
    end
    if (qtl_map_algo == "standard") & (h2_algo in ["minmax", "em_aireml"])
        In = diagm(ones(FloatT, _n)) 
        if n_grms == 1
            glo_KSs = [glo_KS, In]
        end
        if (n_grms == 2) & (glo_h2_model == "Ag+Dg")
            glo_KSs = [_struct_KIN.GRM, _struct_KIN.domGRM, In]
        end
        glo_n_cmp = length(glo_KSs)
    end
    if est_ind_h2 | (qtl_map_algo == "standard")
        _Py = zeros(FloatT, _n)
        _P = zeros(FloatT, _n, _n)
        _PA = zeros(FloatT, _n, _n)
        _Vi = zeros(FloatT, _n, _n) 
        _cVi = zeros(FloatT, _n, _n) 
        _Di = Diagonal(zeros(FloatT, _n, _n))
        _EAvec_Di = zeros(FloatT, _n, _n)
        _EAvec = zeros(FloatT, _n, _n)
        @runif (qtl_map_algo == "standard") & (h2_algo == "em_aireml") glo_Hi = zeros(FloatT, glo_n_cmp, glo_n_cmp)
        cGRM = zeros(FloatT, _n, _n)
        GRM64 = zeros(Float64, _n, _n)
    end
    df_inds = DataFrame()
    sig_df_tops = df_tops_mt[df_tops_mt.qval_g1.<_args_fdr, :]
    if !with_group
        sig_genes_index = pheno_annotation.pheno_id .∈ (sig_df_tops.pheno_id,)
        pheno_annotation.index = 1:nrow(pheno_annotation)
        sig_gene_annot = pheno_annotation[sig_genes_index, :]
        sig_phenotype = phenotype[:, sig_genes_index]
    else 
        sig_genes_index = pheno_annotation.group_id .∈ (sig_df_tops.group_id,)
        pheno_annotation.index = 1:nrow(pheno_annotation)
        sig_gene_annot = pheno_annotation[sig_genes_index, :]
        sig_phenotype = phenotype[:, sig_genes_index]
    end
    println_to_file(string("* ", size(sig_df_tops, 1), " phenotypes (or groups) with qval_g1<", _args_fdr, " will be used for conditional independent analysis."), log_file)
    chroms = intersect(string.(unique(pheno_annotation.chrom)), string.(unique(snp_annotation.chromosome)))
    _task_write_full_pairs = @task @info "Using Tasks"
    for chrom in chroms
        snp_index_chrom = snp_annotation.chromosome .== chrom
        _gene_annot = sig_gene_annot[sig_gene_annot.chrom.==chrom, :]
        _snp_annot = snp_annotation[snp_annotation.chromosome.==chrom, :]
        if with_group 
            _n_genes = sum(_gene_annot.group_id_first) 
        else
            _n_genes = length(_gene_annot.pheno_id) 
        end
        _df_inds = DataFrame()
        df_full = DataFrame()
        _2p = 2 .* snp_annotation.af[snp_index_chrom]
        _2pq = _2p .* (1 .- snp_annotation.af[snp_index_chrom])
        if is_sample_byrow
            _2p = _2p' |> Matrix
            _2pq = _2pq' |> Matrix
            chrom_genotype = similar(genotype[:, snp_index_chrom], FloatT)
            chrom_genotype .= genotype[:, snp_index_chrom]
            if qtl_map_model in ["d+D", "d+A", "d+A+D"]
                chrom_genotype .-= _2pq
            else
                chrom_genotype .-= _2p
            end
            if (_args_run_mode in ["cis_independent"]) & (qtl_map_algo == "idul") & (n_grms == 1)
                chrom_SQ = glo_EA.vectors' * chrom_genotype
                chrom_yQ = glo_EA.vectors' * phenotype[:, _gene_annot.index] 
            end
        else
            chrom_genotype = similar(genotype[snp_index_chrom, :], FloatT)
            chrom_genotype .= genotype[snp_index_chrom, :]
            if qtl_map_model in ["d+D", "d+A", "d+A+D"]
                chrom_genotype .-= _2pq
            else
                chrom_genotype .-= _2p
            end
        end
        _snp_annot.reindex = 1:size(_snp_annot, 1)
        if with_group
            cis_genotype = nothing
            cis_genotype_iterm = nothing
            cis_genotype_2nd = nothing
        end
        GC.gc()
        time_start = now()
        println_to_file(string("Chromosome: ", chrom, ", start at: ", time_start), log_file)
        @runif with_group for i in 1:_n_genes
            gene = unique(_gene_annot.group_id)[i]
            println_to_file(string("    PHENO: ", i, "/", _n_genes, " <", gene,">"), log_file)
            group_pheno_ids = _gene_annot.pheno_id[_gene_annot.group_id.==gene]
            group_size = length(group_pheno_ids)
            index_tested_gene = findfirst(sig_df_tops.group_id .== gene)
            @timeit to "Keep SNPs within cis-region" cissnps_annot, n_cis_snps = get_cis_snp_info(_gene_annot, _snp_annot, group_pheno_ids[1], _args_cis_window)
            @timeit to "Pull cis-SNP genotypes" begin
                if is_sample_byrow
                    if qtl_map_algo == "standard"
                        cis_genotype = chrom_genotype[:, cissnps_annot.reindex]
                    elseif qtl_map_algo == "idul"
                        if n_grms == 1
                            cis_genotype = chrom_SQ[:, cissnps_annot.reindex]
                            @runif _args_run_mode == "cis_interaction" cis_genotype_iterm = chrom_SIQ[:, cissnps_annot.reindex]
                        elseif n_grms == 2
                            cis_genotype = chrom_genotype[:, cissnps_annot.reindex]
                        end
                    end
                else
                    cis_genotype = chrom_genotype[cissnps_annot.reindex, :]
                end
                if (qtl_map_model in ["a+d+A+D"])
                    cis_genotype_2nd = chrom_domGenotype[cissnps_annot.reindex, :]
                end
            end
            @timeit to "Temporary DF to store associations" begin
                df_test = zeros(FloatT, n_cis_snps, 3)
                _df_full = DataFrame([
                    "pheno_id" => "",
                    "variant_id" => cissnps_annot.variant,
                    "rank" => 1,
                    "corr" => NAN,
                ])
                if n_tests == 2
                    _df_full = hcat(_df_full, DataFrame(df_test[:, 1:end-1], Symbol.(["beta_g1", "beta_se_g1", "pval_g1", "beta_g2", "beta_se_g2", "pval_g2"])))
                    if USE_Float32
                        _df_full.pval_g1 = Float64.(_df_full.pval_g1)
                        _df_full.pval_g2 = Float64.(_df_full.pval_g2)
                    end
                end
                if n_tests == 1
                    _df_full = hcat(_df_full, DataFrame(df_test, Symbol.(["beta_g1", "beta_se_g1", "pval_g1"])))
                    if USE_Float32
                        _df_full.pval_g1 = Float64.(_df_full.pval_g1)
                    end
                end
            end
            @timeit to "Prepare for IDUL" if (qtl_map_algo == "idul") & (qtl_map_model in ["d+A+D"])
                if !isnothing(prior_h2)
                    index_prior_h2 = findfirst(prior_h2.pheno_id .== gene)
                    Σ_i = [prior_h2.glo_vg1[index_prior_h2], prior_h2.glo_vg2[index_prior_h2], prior_h2.glo_ve[index_prior_h2]]
                    glo_vc_B = vec(prior_B[index_prior_h2, 3:end])
                end
                ratio = maximum([Σ_i[1] / Σ_i[2], Σ_i[2] / Σ_i[1]])
                @runif _args_verbose println("Ratio between vg1 and vg2: ", ratio) 
                if ratio > maximum(_struct_PartialV.vec_ratio)
                    if ratio > 1000
                        if Σ_i[2] > Σ_i[1] 
                            glo_EA.vectors .= _struct_KIN.domEA.vectors
                            glo_EA.values .= _struct_KIN.domEA.values * Σ_i[2]
                        else 
                            glo_EA.vectors .= _struct_KIN.EA.vectors
                            glo_EA.values .= _struct_KIN.EA.values * Σ_i[1]
                        end
                    else
                        cGRM .= Σ_i[1] * _struct_KIN.GRM + Σ_i[2] * _struct_KIN.domGRM
                        glo_EA = eigen!(cGRM)
                    end
                else
                    ratio_index = argmin(abs.(ratio .- _struct_PartialV.vec_ratio)) * 2
                    if Σ_i[1] >= Σ_i[2]
                        ratio_index -= 1
                    end
                    j = ratio_index
                    range_inds = _n*(j-1)+1:_n*j
                    vals_inds = CartesianIndices((range_inds,))
                    vecs_inds = CartesianIndices((range_inds, 1:_n))
                    @timeit to "p1" copyto!(glo_EA.vectors, Inds_AD_vecs, _struct_PartialV.eigvecs, vecs_inds)
                    @timeit to "p2" glo_EA.values .= _struct_PartialV.eigvals[vals_inds] * sum(Σ_i[1:2])
                end
                @timeit to "p3" cis_genotype .= glo_EA.vectors' * cis_genotype
                @timeit to "p5" mul!(xQ, glo_EA.vectors', X_MME)
            end
            lead_variant = sig_df_tops.variant_id[index_tested_gene]
            index_lead_variant = findfirst(cissnps_annot.variant .== lead_variant)
            index_lead_variant_list = [index_lead_variant]
            if qtl_map_algo == "standard"
                lead_variant_and_X_MME = hcat(cis_genotype[index_lead_variant, :], X_MME)
            elseif qtl_map_algo == "idul"
                lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], xQ)
            end
            ranki = 1
            index_lead_variant_group = index_lead_variant
            best_nominal_pheno_i = -1
            index_lead_pheno_list = [findfirst(group_pheno_ids .== sig_df_tops.pheno_id[sig_df_tops.group_id .== gene])] 
            summary_test_minp = [NaN, NaN, NaN]
            @runif permutation_method != "acat" max_absr_perm = fill(NAN, n_perms)
            @runif permutation_method == "acat" best_p_group_for_acat = fill(NaN, n_cis_snps)
            snp_h2 = snp_h2_se = acc_h2 = acc_h2_se = NAN
            perm_r_mat = fill(NAN, n_perms, n_cis_snps)
            while true
                @runif _args_verbose println("Forward rank: ", ranki)
                pval_g1_threshold = NAN
                if !_args_nominal_only
                    Q, _dof = residualizer(glo_EA.vectors * lead_variant_and_X_MME)
                    QQt = Q * Q'
                end
                _pval_group_min_indices = fill(1, n_cis_snps)
                _df_test_group_max = fill(FloatT(0), n_cis_snps)
                for group_pheno_i in 1:group_size
                    pheno_id = group_pheno_ids[group_pheno_i]
                    chrom_index_tested_pheno = findfirst(_gene_annot.pheno_id .== pheno_id)
                    yQ = chrom_yQ[:, chrom_index_tested_pheno]
                    if qtl_map_algo == "idul"
                        Random.seed!(_args_seed + n_cis_snps)
                        idul_prior_index = sample(1:n_cis_snps, 10)
                        df_subs = fill(NAN, 10, 3 * n_tests + 2)
                        @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, cis_genotype[:, idul_prior_index], yQ, lead_variant_and_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), return_detail=true, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                        init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2]) 
                        @runif isnan(init_eta) begin
                            init_eta = FloatT(1e-5)
                            println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                            @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                        end
                        if !_args_exact_map
                            @timeit to "Idul test" idul_assoc_test_approx!(df_test, cis_genotype, yQ, lead_variant_and_X_MME, glo_EA.values; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                        else
                            @timeit to "Idul test" idul_assoc_test_plus!(df_test, cis_genotype, yQ, lead_variant_and_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                        end
                    end
                    greater_indices = findall(df_test[:,3] .> _df_test_group_max)
                    _pval_group_min_indices[greater_indices] .= group_pheno_i
                    _df_test_group_max[greater_indices] .= df_test[greater_indices,3]
                    index_lead_variant = first_nonnan_argmax(df_test[:, 3])
                    @runif index_lead_variant == -1 continue
                    if group_pheno_i == 1
                        summary_test_minp .= df_test[index_lead_variant, :]
                        summary_test_minp[3] = ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), summary_test_minp[3])
                        index_lead_variant_group = index_lead_variant
                        best_nominal_pheno_i = 1
                    else
                        current_group_lead_pval = ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), df_test[index_lead_variant, 3])
                        if current_group_lead_pval < summary_test_minp[3]
                            summary_test_minp .= df_test[index_lead_variant, :]
                            summary_test_minp[3] = current_group_lead_pval
                            index_lead_variant_group = index_lead_variant
                            best_nominal_pheno_i = group_pheno_i
                        end
                    end
                    if n_grms < 2
                        cis_genotype .= chrom_SQ[:, cissnps_annot.reindex]
                    else
                        cis_genotype .= glo_EA.vectors' * chrom_genotype[:, cissnps_annot.reindex]
                    end
                end
                perm_r_mat .= NAN
                for group_pheno_i in 1:group_size
                    pheno_id = group_pheno_ids[group_pheno_i]
                    global_index_tested_pheno = _gene_annot.index[_gene_annot.pheno_id.==pheno_id]
                    if !_args_nominal_only
                        if permutation_method != "acat"
                            cis_genotype .= chrom_genotype[:, cissnps_annot.reindex]
                            get_matrix_resid!(cis_genotype, QQt, return_std=false)
                            variant_subset_indices = findall(_pval_group_min_indices .== group_pheno_i)
                            if permutation_method == "standard_fast"
                                y_perms .= phenotype[:, global_index_tested_pheno] .- Statistics.mean(phenotype[:, global_index_tested_pheno])
                                Random.seed!(_args_seed)
                                @inbounds @views for perm_i in 1:n_perms
                                    shuffle!(y_perms[:, perm_i])
                                end
                                get_matrix_resid!(y_perms, QQt, return_std=false)
                                perm_r_mat[:, variant_subset_indices] = fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm, variant_subset_indices)
                            elseif permutation_method == "clipper"
                                pheno_perm = vec(get_matrix_resid(phenotype[:, global_index_tested_pheno], QQt, return_std=false)) 
                                @views fast_permutation_clipper!(perm_r_mat[:, variant_subset_indices], cis_genotype, pheno_perm, n_perms, variant_subset_indices)
                            end
                            if group_pheno_i == group_size
                                max_absr_perm .= vec(maximum(abs, perm_r_mat, dims=2))
                            end
                        end
                        if permutation_method == "acat"
                            if group_pheno_i == 1
                                best_p_group_for_acat .= df_test[:, 3]
                            else
                                best_p_group_for_acat[best_p_group_for_acat.<df_test[:, 3]] .= df_test[:, 3][best_p_group_for_acat.<df_test[:, 3]]
                            end
                        end
                    end
                end
                if n_grms < 2
                    cis_genotype .= chrom_SQ[:, cissnps_annot.reindex]
                else
                    cis_genotype .= glo_EA.vectors' * chrom_genotype[:, cissnps_annot.reindex]
                end
                if !isnothing(_args_ind_r2)
                    snp_r2 = abs2.(vec(cor(lead_variant_and_X_MME[:, 1:ranki], cis_genotype[:, index_lead_variant_group])))
                    if any(snp_r2 .> _args_ind_r2)
                        break
                    end
                end
                if permutation_method == "acat"
                    best_p_group_for_acat .= ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), best_p_group_for_acat)
                    pval_acat = ACATest(best_p_group_for_acat, is_check=false)
                    if (pval_acat <= _args_fdr) & !(index_lead_variant_group in index_lead_variant_list)
                        ranki = ranki + 1
                        index_lead_variant = index_lead_variant_group
                        push!(index_lead_variant_list, index_lead_variant)
                        push!(index_lead_pheno_list, best_nominal_pheno_i)
                        lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], lead_variant_and_X_MME)
                    else
                        break
                    end
                else
                    pval_g1_threshold = sort(get_approx_p_from_r(max_absr_perm, _dof))[ceil(Int, n_perms * _args_fdr)]
                    if (summary_test_minp[3] <= pval_g1_threshold) & !(index_lead_variant_group in index_lead_variant_list)
                        ranki = ranki + 1
                        index_lead_variant = index_lead_variant_group
                        push!(index_lead_variant_list, index_lead_variant)
                        push!(index_lead_pheno_list, best_nominal_pheno_i)
                        lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], lead_variant_and_X_MME)
                    else
                        break
                    end
                end
            end
            n_forward_variants = length(index_lead_variant_list)
            if n_forward_variants > 1
                backward_X_MME = fill(NAN, _n, size(lead_variant_and_X_MME, 2) - 1)
            else
                backward_X_MME = xQ
            end
            index_lead_variant_list_pass = []
            ranki = 1
            _df_test = copy(df_test)
            for bi in 1:n_forward_variants
                @runif _args_verbose println("Backward rank: ", bi)
                rm_col = n_forward_variants - bi + 1
                @runif n_forward_variants > 1 backward_X_MME .= lead_variant_and_X_MME[:, setdiff(1:size(lead_variant_and_X_MME, 2), rm_col)]
                index_lead_variant = index_lead_variant_list[bi]
                best_indsnp_pheno_index = 0
                if !_args_nominal_only
                    Q, _dof = residualizer(glo_EA.vectors * backward_X_MME)
                    QQt = Q * Q'
                end
                group_pheno_i = index_lead_pheno_list[bi]
                pheno_id = group_pheno_ids[group_pheno_i]
                chrom_index_tested_pheno = findfirst(_gene_annot.pheno_id .== pheno_id)
                global_index_tested_pheno = _gene_annot.index[_gene_annot.pheno_id.==pheno_id]
                yQ = chrom_yQ[:, chrom_index_tested_pheno]
                if qtl_map_algo == "idul"
                    Random.seed!(_args_seed + n_cis_snps)
                    idul_prior_index = sample(1:n_cis_snps, 10)
                    df_subs = fill(FloatT(NAN), 10, 3 * n_tests + 2)
                    @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, cis_genotype[:, idul_prior_index], yQ, backward_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), return_detail=true, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                    init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                    @runif isnan(init_eta) begin
                        init_eta = FloatT(1e-5)
                        println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                        @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                    end
                    if !_args_exact_map
                        @timeit to "Idul test" idul_assoc_test_approx!(df_test, cis_genotype, yQ, backward_X_MME, glo_EA.values; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                    else
                        @timeit to "Idul test" idul_assoc_test_plus!(df_test, cis_genotype, yQ, backward_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                    end
                end
                    summary_test_minp .= df_test[index_lead_variant, :]
                    summary_test_minp[3] = ccdf(WaldTest(1, _n - size(backward_X_MME, 2) - 1), summary_test_minp[3])
                    _df_test .= df_test
                    best_indsnp_pheno_index = chrom_index_tested_pheno
                pval_g1_threshold = NAN
                if !_args_nominal_only
                    if permutation_method != "acat"
                        cis_genotype .= chrom_genotype[:, cissnps_annot.reindex]
                        get_matrix_resid!(cis_genotype, QQt, return_std=false)
                        if permutation_method == "standard_fast"
                            y_perms .= phenotype[:, global_index_tested_pheno] .- Statistics.mean(phenotype[:, global_index_tested_pheno])
                            Random.seed!(_args_seed)
                            @inbounds @views for perm_i in 1:n_perms
                                shuffle!(y_perms[:, perm_i])
                            end
                            get_matrix_resid!(y_perms, QQt, return_std=false)
                            fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm, byrow=is_sample_byrow) 
                        elseif permutation_method == "clipper"
                            pheno_perm = vec(get_matrix_resid(phenotype[:, global_index_tested_pheno], QQt, return_std=false)) 
                            fast_permutation_clipper!(absr_perm, cis_genotype, pheno_perm, n_perms; byrow=is_sample_byrow)
                        end
                            max_absr_perm .= absr_perm
                    end
                    if permutation_method == "acat"
                            best_p_group_for_acat .= df_test[:, 3]
                    end
                end
                if n_grms < 2
                    cis_genotype .= chrom_SQ[:, cissnps_annot.reindex]
                else
                    cis_genotype .= glo_EA.vectors' * chrom_genotype[:, cissnps_annot.reindex]
                end
                if permutation_method == "acat"
                    best_p_group_for_acat .= ccdf(WaldTest(1, _n - size(backward_X_MME, 2) - 1), best_p_group_for_acat)
                    pval_acat = ACATest(best_p_group_for_acat, is_check=false)
                    go_append = (pval_acat <= _args_fdr) | (n_forward_variants == 1)
                else
                    pval_g1_threshold = sort(get_approx_p_from_r(max_absr_perm, _dof))[ceil(Int, n_perms * _args_fdr)]
                    go_append = (summary_test_minp[3] <= pval_g1_threshold) | (n_forward_variants == 1)
                end
                if go_append
                    snp_r = cor(lead_variant_and_X_MME[:, rm_col], lead_variant_and_X_MME[:, n_forward_variants])
                    append!(index_lead_variant_list_pass, index_lead_variant)
                    variant_id = cissnps_annot.variant[index_lead_variant]
                    start_distance = cissnps_annot.start_distance[index_lead_variant]
                    _af = cissnps_annot.af[index_lead_variant]
                    het_rate = cissnps_annot.het_rate[index_lead_variant]
                    beta_g1, beta_se_g1, pval_g1 = summary_test_minp
                    _df_backward = DataFrame([
                        "chrom" => chrom,
                        "group_id" => gene,
                        "pheno_id" => _gene_annot.pheno_id[best_indsnp_pheno_index],
                        "variant_id" => variant_id,
                        "start_distance" => start_distance,
                        "af" => _af,
                        "het_rate" => het_rate,
                        "beta_g1" => beta_g1,
                        "beta_se_g1" => beta_se_g1,
                        "pval_g1" => pval_g1,
                        "pval_g1_threshold" => pval_g1_threshold,
                        "rank" => ranki,
                        "corr" => snp_r,
                    ])
                    @runif !_args_nominal_only begin
                        cols_float = [eltype(_df_backward[:, x]) == FloatT for x in range(1, ncol(_df_backward))]
                        _df_backward[:, cols_float] = round.(_df_backward[:, cols_float], sigdigits=6)
                        is_append = (size(_df_inds, 1) > 0) | (size(df_inds, 1) > 0) ? true : false
                        if length(_args_chrom) > 0
                            CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_independent.", chrom, ".txt")), _df_backward, delim="\t", append=is_append)
                        else
                            CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_independent.txt.gz")), _df_backward, delim="\t", compress=true, append=is_append)
                        end
                    end
                    append!(_df_inds, _df_backward)
                    snp_rs = vec(cor(cis_genotype, cis_genotype[:, index_lead_variant]))
                    @timeit to "append the full summary" begin
                        _df_full.pheno_id .= _gene_annot.pheno_id[best_indsnp_pheno_index]
                        _df_full.rank .= ranki
                        _df_full.corr .= round.(snp_rs, sigdigits=3)
                        _df_full.beta_g1 .= _df_test[:, 1]
                        _df_full.beta_se_g1 .= _df_test[:, 2]
                        _df_full.pval_g1 .= ccdf(WaldTest(1, _n - size(backward_X_MME, 2) - 1), df_test[:, 3])
                        append!(df_full, _df_full) 
                    end
                    ranki += 1
                end
            end
            @runif _args_verbose if i % 10 == 0
                println_to_file(string("*** Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
                show(to)
                println()
            end
        end
        @runif !with_group for i in 1:_n_genes
            gene = _gene_annot.pheno_id[i]
            println_to_file(string("    PHENO: ", i, "/", _n_genes, " <", gene,">"), log_file)
            index_tested_gene = findfirst(sig_df_tops.pheno_id .== gene)
            global_index_tested_gene = _gene_annot.index[findfirst(_gene_annot.pheno_id .== gene)]
            @timeit to "Keep SNPs within cis-region" cissnps_annot, n_cis_snps = get_cis_snp_info(_gene_annot, _snp_annot, gene, _args_cis_window)
            @timeit to "Pull cis-SNP genotypes" begin
                if is_sample_byrow
                    if qtl_map_algo == "standard"
                        cis_genotype = chrom_genotype[:, cissnps_annot.reindex]
                    elseif qtl_map_algo == "idul"
                        if n_grms == 1
                            cis_genotype = chrom_SQ[:, cissnps_annot.reindex]
                            @runif _args_run_mode == "cis_interaction" cis_genotype_iterm = chrom_SIQ[:, cissnps_annot.reindex]
                        elseif n_grms == 2
                            cis_genotype = chrom_genotype[:, cissnps_annot.reindex]
                        end
                    end
                else
                    cis_genotype = chrom_genotype[cissnps_annot.reindex, :]
                end
                if (qtl_map_model in ["a+d+A+D"])
                    cis_genotype_2nd = chrom_domGenotype[cissnps_annot.reindex, :]
                end
            end
            exppheno = vec(sig_phenotype[:, sig_gene_annot.pheno_id.==gene]) 
            @timeit to "Temporary DF to store associations" begin
                df_test = zeros(FloatT, n_cis_snps, 3)
                _df_full = DataFrame([
                    "pheno_id" => gene,
                    "variant_id" => cissnps_annot.variant,
                    "rank" => 1,
                    "corr" => NAN,
                ])
                if n_tests == 2
                    _df_full = hcat(_df_full, DataFrame(df_test[:, 1:end-1], Symbol.(["beta_g1", "beta_se_g1", "pval_g1", "beta_g2", "beta_se_g2", "pval_g2"])))
                    if USE_Float32
                        _df_full.pval_g1 = Float64.(_df_full.pval_g1)
                        _df_full.pval_g2 = Float64.(_df_full.pval_g2)
                    end
                end
                if n_tests == 1
                    _df_full = hcat(_df_full, DataFrame(df_test, Symbol.(["beta_g1", "beta_se_g1", "pval_g1"])))
                    if USE_Float32
                        _df_full.pval_g1 = Float64.(_df_full.pval_g1)
                    end
                end
            end
            @timeit to "Prepare for IDUL" if (qtl_map_algo == "idul") & (qtl_map_model in ["d+A+D"])
                if !isnothing(prior_h2)
                    index_prior_h2 = findfirst(prior_h2.pheno_id .== gene)
                    Σ_i = [prior_h2.glo_vg1[index_prior_h2], prior_h2.glo_vg2[index_prior_h2], prior_h2.glo_ve[index_prior_h2]]
                    glo_vc_B = vec(prior_B[index_prior_h2, 3:end])
                end
                ratio = maximum([Σ_i[1] / Σ_i[2], Σ_i[2] / Σ_i[1]])
                @runif _args_verbose println("Ratio between vg1 and vg2: ", ratio) 
                if ratio > maximum(_struct_PartialV.vec_ratio)
                    if ratio > 1000
                        if Σ_i[2] > Σ_i[1] 
                            glo_EA.vectors .= _struct_KIN.domEA.vectors
                            glo_EA.values .= _struct_KIN.domEA.values * Σ_i[2]
                        else 
                            glo_EA.vectors .= _struct_KIN.EA.vectors
                            glo_EA.values .= _struct_KIN.EA.values * Σ_i[1]
                        end
                    else
                        cGRM .= Σ_i[1] * _struct_KIN.GRM + Σ_i[2] * _struct_KIN.domGRM
                        glo_EA = eigen!(cGRM)
                    end
                else
                    ratio_index = argmin(abs.(ratio .- _struct_PartialV.vec_ratio)) * 2
                    if Σ_i[1] >= Σ_i[2]
                        ratio_index -= 1
                    end
                    j = ratio_index
                    range_inds = _n*(j-1)+1:_n*j
                    vals_inds = CartesianIndices((range_inds,))
                    vecs_inds = CartesianIndices((range_inds, 1:_n))
                    @timeit to "p1" copyto!(glo_EA.vectors, Inds_AD_vecs, _struct_PartialV.eigvecs, vecs_inds)
                    @timeit to "p2" glo_EA.values .= _struct_PartialV.eigvals[vals_inds] * sum(Σ_i[1:2])
                end
                @timeit to "p3" cis_genotype .= glo_EA.vectors' * cis_genotype
                @timeit to "p4" mul!(yQ, glo_EA.vectors', exppheno)
                @timeit to "p5" mul!(xQ, glo_EA.vectors', X_MME)
            end
            lead_variant = sig_df_tops.variant_id[index_tested_gene]
            index_lead_variant = findfirst(cissnps_annot.variant .== lead_variant)
            index_lead_variant_list = [index_lead_variant]
            if qtl_map_algo == "standard"
                lead_variant_and_X_MME = hcat(cis_genotype[index_lead_variant, :], X_MME)
            elseif qtl_map_algo == "idul"
                if n_grms == 1
                    yQ = chrom_yQ[:, i]
                end
                lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], xQ)
            end
            ranki = 1
            summary_test_minp = [NaN, NaN, NaN]
            snp_h2 = snp_h2_se = acc_h2 = acc_h2_se = NAN
            while true
                @runif _args_verbose println("Forward rank: ", ranki)
                if qtl_map_algo == "idul"
                    Random.seed!(_args_seed + n_cis_snps)
                    idul_prior_index = sample(1:n_cis_snps, 10)
                    df_subs = fill(FloatT(NAN), 10, 3 * n_tests + 2)
                    @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, cis_genotype[:, idul_prior_index], yQ, lead_variant_and_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), return_detail=true, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                    init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                    @runif isnan(init_eta) begin
                        init_eta = FloatT(1e-5)
                        println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                        @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                    end
                    if !_args_exact_map
                        @timeit to "Idul test" idul_assoc_test_approx!(df_test, cis_genotype, yQ, lead_variant_and_X_MME, glo_EA.values; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                    else
                        @timeit to "Idul test" idul_assoc_test_plus!(df_test, cis_genotype, yQ, lead_variant_and_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                    end
                end
                index_lead_variant = first_nonnan_argmax(df_test[:, 3])
                @runif index_lead_variant == -1 break
                summary_test_minp .= df_test[index_lead_variant, :]
                summary_test_minp[3] = ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), summary_test_minp[3])
                pval_g1_threshold = NAN
                if permutation_method == "acat"
                    pvals_for_acat = ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), df_test[:, 3])
                    pval_acat = ACATest(pvals_for_acat; is_check=false)
                else
                    Q, _dof = residualizer(glo_EA.vectors * lead_variant_and_X_MME) 
                    QQt = Q * Q'
                    cis_genotype .= chrom_genotype[:, cissnps_annot.reindex]
                    get_matrix_resid!(cis_genotype, QQt, return_std=false)
                    if permutation_method == "standard_fast"
                        y_perms .= phenotype[:, global_index_tested_gene] .- Statistics.mean(phenotype[:, global_index_tested_gene])
                        Random.seed!(_args_seed)
                        @inbounds @views for perm_i in 1:n_perms
                            shuffle!(y_perms[:, perm_i])
                        end
                        get_matrix_resid!(y_perms, QQt, return_std=false)
                        fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm, byrow=is_sample_byrow) 
                    elseif permutation_method == "clipper"
                        pheno_perm = vec(get_matrix_resid(phenotype[:, global_index_tested_gene], QQt, return_std=false)) 
                        fast_permutation_clipper!(absr_perm, cis_genotype, pheno_perm, n_perms; byrow=is_sample_byrow)
                    end
                    pval_g1_threshold = sort(get_approx_p_from_r(absr_perm, _dof))[ceil(Int, n_perms * _args_fdr)]
                end
                if n_grms < 2
                    cis_genotype .= chrom_SQ[:, cissnps_annot.reindex]
                else
                    cis_genotype .= glo_EA.vectors' * chrom_genotype[:, cissnps_annot.reindex]
                end
                if !isnothing(_args_ind_r2)
                    snp_r2 = abs2.(vec(cor(lead_variant_and_X_MME[:, 1:ranki], cis_genotype[:, index_lead_variant])))
                    if any(snp_r2 .> _args_ind_r2)
                        break
                    end
                end
                if permutation_method == "acat"
                    if (pval_acat <= _args_fdr) & !(index_lead_variant in index_lead_variant_list)
                        ranki = ranki + 1
                        push!(index_lead_variant_list, index_lead_variant)
                        lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], lead_variant_and_X_MME)
                    else
                        break
                    end
                else
                    if (summary_test_minp[3] <= pval_g1_threshold) & !(index_lead_variant in index_lead_variant_list)
                        ranki = ranki + 1
                        push!(index_lead_variant_list, index_lead_variant)
                        lead_variant_and_X_MME = hcat(cis_genotype[:, index_lead_variant], lead_variant_and_X_MME)
                    else
                        break
                    end
                end
            end
            n_forward_variants = length(index_lead_variant_list)
            if n_forward_variants > 1
                backward_X_MME = fill(NAN, _n, size(lead_variant_and_X_MME, 2) - 1)
            else
                backward_X_MME = xQ
            end
            index_lead_variant_list_pass = []
            ranki = 1
            for bi in 1:n_forward_variants
                @runif _args_verbose println("Backward rank: ", bi)
                rm_col = n_forward_variants - bi + 1
                @runif n_forward_variants > 1 backward_X_MME .= lead_variant_and_X_MME[:, setdiff(1:size(lead_variant_and_X_MME, 2), rm_col)]
                index_lead_variant = index_lead_variant_list[bi]
                if qtl_map_algo == "standard"
                    @timeit to "Adjusted phenotype" yadj = exppheno - lead_variant_and_X_MME * glo_vc[:Fix_eff]
                    @timeit to "Vi" if qtl_map_model in ["a+A", "d+D", "d+A"]
                        getVinv!(_Vi, _Di, glo_EA, glo_Σ_i, _EAvec_Di)
                    elseif qtl_map_model in ["a+A+D", "d+A+D", "a+d+A+D"]
                        getVinv!(_Vi, glo_Σ_i, _struct_PartialV, _Di, _EAvec_Di, _EAvec)
                    end
                    @timeit to "Get associations" calcu_one_assoc_test!(df_test, _Vi, cis_genotype, yadj)
                elseif qtl_map_algo == "idul"
                    Random.seed!(_args_seed + n_cis_snps)
                    idul_prior_index = sample(1:n_cis_snps, 10)
                    df_subs = fill(FloatT(NAN), 10, 3 * n_tests + 2)
                    @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, cis_genotype[:, idul_prior_index], yQ, backward_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), return_detail=true, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                    init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                    @runif isnan(init_eta) begin
                        init_eta = FloatT(1e-5)
                        println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                        @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                    end
                    if !_args_exact_map
                        @timeit to "Idul test" idul_assoc_test_approx!(df_test, cis_genotype, yQ, backward_X_MME, glo_EA.values; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                    else
                        @timeit to "Idul test" idul_assoc_test_plus!(df_test, cis_genotype, yQ, backward_X_MME, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                    end
                end
                summary_test_minp .= df_test[index_lead_variant, :]
                beta_g1, beta_se_g1, pval_g1 = summary_test_minp
                pval_g1 = ccdf(WaldTest(1, _n - size(backward_X_MME, 2) - 1), pval_g1)
                pval_g1_threshold = NAN
                if permutation_method == "acat"
                    pvals_for_acat = ccdf(WaldTest(1, _n - size(lead_variant_and_X_MME, 2) - 1), df_test[:, 3])
                    pval_acat = ACATest(pvals_for_acat; is_check=false)
                else
                    Q, _dof = residualizer(glo_EA.vectors * backward_X_MME)
                    QQt = Q * Q'
                    cis_genotype .= chrom_genotype[:, cissnps_annot.reindex]
                    get_matrix_resid!(cis_genotype, QQt, return_std=false)
                    if permutation_method == "standard_fast"
                        y_perms .= phenotype[:, global_index_tested_gene] .- Statistics.mean(phenotype[:, global_index_tested_gene])
                        Random.seed!(_args_seed)
                        @inbounds @views for perm_i in 1:n_perms
                            shuffle!(y_perms[:, perm_i])
                        end
                        get_matrix_resid!(y_perms, QQt, return_std=false)
                        fast_linear_permutation_test!(cis_genotype, y_perms, absr_perm, byrow=is_sample_byrow) 
                    elseif permutation_method == "clipper"
                        pheno_perm = vec(get_matrix_resid(phenotype[:, global_index_tested_gene], QQt, return_std=false)) 
                        fast_permutation_clipper!(absr_perm, cis_genotype, pheno_perm, n_perms; byrow=is_sample_byrow)
                    end
                    pval_g1_threshold = sort(get_approx_p_from_r(absr_perm, _dof))[ceil(Int, n_perms * _args_fdr)]
                end
                if n_grms < 2
                    cis_genotype .= chrom_SQ[:, cissnps_annot.reindex]
                else
                    cis_genotype .= glo_EA.vectors' * chrom_genotype[:, cissnps_annot.reindex]
                end
                snp_r = cor(lead_variant_and_X_MME[:, rm_col], lead_variant_and_X_MME[:, n_forward_variants])
                skip_current_variant = false
                if permutation_method == "acat"
                    if (pval_acat <= _args_fdr) | (n_forward_variants == 1)
                        append!(index_lead_variant_list_pass, index_lead_variant)
                    else
                        skip_current_variant = true
                    end
                else
                    if (pval_g1 <= pval_g1_threshold) | (n_forward_variants == 1)
                        append!(index_lead_variant_list_pass, index_lead_variant)
                    else
                        skip_current_variant = true
                    end
                end
                if skip_current_variant
                    continue
                end
                if est_ind_h2 
                    @timeit to "SNP heritability" begin
                        if qtl_map_algo == "idul"
                            geno_lead_variant = glo_EA.vectors * cis_genotype[:, index_lead_variant:index_lead_variant] 
                        end
                        if qtl_map_model in ["d+D", "d+A+D"]
                            getG_dominance!(GRM64, Float64.(geno_lead_variant), code_type="dominance", adjusted=true, afs=Float64.([cissnps_annot.af[index_lead_variant]]), byrow=true)
                        elseif qtl_map_model in ["a+A", "a+A+D", "d+A"]
                            getG_VanRaden!(GRM64, Float64.(geno_lead_variant), adjusted=true, afs=Float64.([cissnps_annot.af[index_lead_variant]]), byrow=true)
                        end
                        cGRM .= GRM64
                        _Vi .= cGRM
                        snp_GRM_EA = eigen!(_Vi)
                        if h2_algo == "minmax"
                            snp_vc = getMRVCModel([cGRM, In], exppheno; X=X_MME)
                        elseif (h2_algo == "idul") | true
                            snp_xQ = snp_GRM_EA.vectors' * X_MME
                            snp_vc = get_IDUL_VarianceComponent(snp_GRM_EA, exppheno, X_MME, xQ=snp_xQ)
                        end
                        snp_h2 = snp_vc[:h2]
                        snp_h2_se = snp_vc[:h2_se]
                    end
                    @timeit to "Acc SNP heritability" begin
                        if ranki > 1 
                            if qtl_map_algo == "idul"
                                geno_lead_variants = glo_EA.vectors * cis_genotype[:, index_lead_variant_list_pass] 
                            end
                            if qtl_map_model in ["d+D", "d+A+D"]
                                getG_dominance!(GRM64, Float64.(geno_lead_variants), code_type="dominance", adjusted=true, afs=Float64.(cissnps_annot.af[index_lead_variant_list]), byrow=true)
                            elseif qtl_map_model in ["a+A", "a+A+D", "d+A"]
                                getG_VanRaden!(GRM64, Float64.(geno_lead_variants), adjusted=true, afs=Float64.(cissnps_annot.af[index_lead_variant_list]), byrow=true)
                            end
                            cGRM .= GRM64
                            _Vi .= cGRM
                            acc_GRM_EA = eigen!(_Vi)
                            if h2_algo == "minmax"
                                acc_vc = getMRVCModel([cGRM, In], exppheno; X=X_MME)
                            elseif (h2_algo == "idul") | true
                                acc_xQ = acc_GRM_EA.vectors' * X_MME
                                acc_vc = get_IDUL_VarianceComponent(acc_GRM_EA, exppheno, X_MME, xQ=acc_xQ)
                            end
                            acc_h2 = acc_vc[:h2]
                            acc_h2_se = acc_vc[:h2_se]
                        else
                            acc_h2 = snp_h2
                            acc_h2_se = snp_h2_se
                        end
                    end
                end
                variant_id = cissnps_annot.variant[index_lead_variant]
                start_distance = cissnps_annot.start_distance[index_lead_variant]
                _af = cissnps_annot.af[index_lead_variant]
                het_rate = cissnps_annot.het_rate[index_lead_variant]
                if est_ind_h2
                    _df_backward = DataFrame([
                        "chrom" => chrom,
                        "pheno_id" => gene,
                        "variant_id" => variant_id,
                        "start_distance" => start_distance,
                        "af" => _af,
                        "het_rate" => het_rate,
                        "beta_g1" => beta_g1,
                        "beta_se_g1" => beta_se_g1,
                        "pval_g1" => pval_g1,
                        "pval_g1_threshold" => pval_g1_threshold,
                        "rank" => ranki,
                        "corr" => snp_r,
                        "snp_h2" => snp_h2,
                        "snp_h2_se" => snp_h2_se,
                        "acc_h2" => acc_h2,
                        "acc_h2_se" => acc_h2_se,
                    ])
                else
                    _df_backward = DataFrame([
                        "chrom" => chrom,
                        "pheno_id" => gene,
                        "variant_id" => variant_id,
                        "start_distance" => start_distance,
                        "af" => _af,
                        "het_rate" => het_rate,
                        "beta_g1" => beta_g1,
                        "beta_se_g1" => beta_se_g1,
                        "pval_g1" => pval_g1,
                        "pval_g1_threshold" => pval_g1_threshold,
                        "rank" => ranki,
                        "corr" => snp_r,
                    ])
                end
                @runif !_args_nominal_only begin
                    cols_float = [eltype(_df_backward[:, x]) == FloatT for x in range(1, ncol(_df_backward))]
                    _df_backward[:, cols_float] = round.(_df_backward[:, cols_float], sigdigits=6)
                    is_append = (size(_df_inds, 1) > 0) | (size(df_inds, 1) > 0) ? true : false
                    if length(_args_chrom) > 0
                        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_independent.", chrom, ".txt")), _df_backward, delim="\t", append=is_append)
                    else
                        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_independent.txt.gz")), _df_backward, delim="\t", compress=true, append=is_append)
                    end
                end
                append!(_df_inds, _df_backward)
                snp_rs = vec(cor(cis_genotype, cis_genotype[:, index_lead_variant]))
                @timeit to "append the full summary" begin
                    _df_full.rank .= ranki
                    _df_full.corr .= round.(snp_rs, sigdigits=3)
                    _df_full.beta_g1 .= df_test[:, 1]
                    _df_full.beta_se_g1 .= df_test[:, 2]
                    _df_full.pval_g1 .= ccdf(WaldTest(1, _n - size(backward_X_MME, 2) - 1), df_test[:, 3]) 
                    append!(df_full, _df_full) 
                end
                ranki += 1
            end
            @runif _args_verbose if i % 10 == 0
                println_to_file(string("*** Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
                show(to)
                println()
            end
        end
        cols_float = [df_full[1, x] isa AbstractFloat for x in range(1, ncol(df_full))]
        df_full[:, cols_float] = round.(df_full[:, cols_float], sigdigits=6)
        @timeit to "Write the full summary to disk" if !_args_write_top
            if chrom != chroms[1]
                wait(_task_write_full_pairs) 
            end
            task_df_full = copy(df_full)
            task_chrom = chrom
            _task_write_full_pairs = @spawn begin
                chrom_full_out_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_independent_pairs.", task_chrom, ".txt.gz"))
                @runif _args_use_gzip chrom_full_out_file = replace(chrom_full_out_file, r".gz$" => "")
                CSV.write(chrom_full_out_file, task_df_full, delim="\t", compress=!_args_use_gzip)
                @runif _args_use_gzip run(`gzip -f $chrom_full_out_file`; wait=chrom == chroms[end])
            end
            if chrom == chroms[end]
                wait(_task_write_full_pairs)
            end
        end
        println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
        append!(df_inds, _df_inds)
        if _args_verbose
            println_to_file(string(to), log_file)
            println()
        end
    end
    if _args_verbose
        show(to)
    end
end
