function runOmiGA_cis(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM::Union{Nothing,Dominance}=nothing, _struct_ITERM::Union{Nothing,InteractionTerm}=nothing)
    phenotype = _struct_PHENO.phenotype
    pheno_annotation = _struct_PHENO.annotation 
    X_MME = _struct_COVAR.X_MME
    @runif USE_GPU & _args_linear_model X_MME = X_MME |> CuArray
    if !isnothing(_struct_PRIOR_HerB)
        prior_h2 = _struct_PRIOR_HerB.h2
        prior_B = _struct_PRIOR_HerB.B
    else
        prior_h2 = prior_B = nothing
    end
    snp_annotation = _struct_GENO.annotation
    qtl_map_model = _args_qtl_map_model
    qtl_map_algo = _args_qtl_map_algo
    _n_samples = _struct_GENO.n_samples
    _n, _X_c = size(X_MME)
    with_strand = issubset(["strand"], names(pheno_annotation))
    test_phen_type = "raw" 
    test_geno_type = "raw" 
    perm_phen_type = "raw" 
    perm_geno_type = "raw" 
    n_tests = 1 
    if _args_linear_model
        n_grms = 0
        qtl_map_algo = "standard_fast" 
        if qtl_map_model == "a"
            genotype = _struct_GENO.genotype
            test_phen_type = "residualized"
            test_geno_type = "residualized"
        end
        if qtl_map_model == "d"
            genotype = _struct_DOM.genotype
            test_phen_type = "residualized"
            test_geno_type = "residualized"
        end
        if qtl_map_model == "a+d" 
            error("'a+d' is not yet supported!")
            genotype = _struct_GENO.genotype
            domGenotype = _struct_DOM.genotype
            n_tests = 2
        end
        if _args_run_mode == "cis_interaction" 
            qtl_map_algo = "idul"
            glo_EA = eigen(diagm(ones(FloatT, _n_samples))) 
            test_phen_type = "rotate"
            test_geno_type = "rotate"
            xQ = X_MME
            n_tests = _struct_ITERM.n_iterms + ifelse(gxe_omit_main_eff, 0, 1)
        end
    else
        n_grms = 1
        if qtl_map_model == "a+A"
            glo_EA = _struct_KIN.EA
            genotype = _struct_GENO.genotype
            if isnothing(prior_h2) & (_args_h2_algo == "minmax")
                glo_KS = _struct_KIN.GRM
            end
        end
        if qtl_map_model == "d+D"
            glo_EA = _struct_KIN.domEA
            @runif USE_GPU glo_EA = Eigen{FloatT,FloatT,CuMatrix{FloatT},CuVector{FloatT}}(FloatT.(glo_EA.values), FloatT.(glo_EA.vectors))
            genotype = _struct_DOM.genotype
            if isnothing(prior_h2) & (_args_h2_algo == "minmax")
                glo_KS = _struct_KIN.domGRM
            end
        end
        if qtl_map_model == "d+A"
            glo_EA = _struct_KIN.EA
            @runif USE_GPU glo_EA = Eigen{FloatT,FloatT,CuMatrix{FloatT},CuVector{FloatT}}(FloatT.(glo_EA.values), FloatT.(glo_EA.vectors))
            genotype = _struct_DOM.genotype
            if isnothing(prior_h2) & (_args_h2_algo == "minmax")
                glo_KS = _struct_KIN.GRM
            end
        end
        if qtl_map_model == "a+d+A+D"
            genotype = _struct_GENO.genotype
            domGenotype = _struct_DOM.genotype
            n_tests = 2
            n_grms = 2
        end
        if qtl_map_model == "d+A+D"
            genotype = _struct_DOM.genotype
            n_grms = 2
        end
        if (qtl_map_algo == "idul") & (n_grms < 2)
            test_phen_type = "rotate"
            test_geno_type = "rotate"
        end
        if _args_run_mode in ["cis_interaction"]
            n_tests = _struct_ITERM.n_iterms + ifelse(gxe_omit_main_eff, 0, 1)
            test_phen_type = "rotate"
            test_geno_type = "rotate"
        end
        _args_PreV = n_grms == 2
        @timeit to "Pre-calcu V" if _args_PreV
            path_partialv = _args_path_partialv
            if !isnothing(path_partialv)
                _struct_PartialV = load(path_partialv)["data"]
            else
                error("--paed-file must be specified for QTL mapping models with two GRMs!")
            end
            if isnothing(prior_h2)
                error("--her-file must be specified for QTL mapping models with two GRMs!")
            end
        end
        if isnothing(prior_h2) & (_args_h2_algo == "idul") & (n_grms < 2)
            test_phen_type = "rotate"
            test_geno_type = "rotate"
        end
        if test_phen_type == "rotate"
            if _args_preadj_covar
                xQ = glo_EA.vectors' * X_MME[:, 1:1]
            else
                xQ = glo_EA.vectors' * X_MME 
            end
        end
        if (qtl_map_algo == "idul") & (n_grms == 2)
            xQ = similar(X_MME)
            yQ = zeros(FloatT, _n_samples)
            cGRM = zeros(FloatT, _n_samples, _n_samples)
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
    end
    if _X_c == 1
        QQt = nothing
        Q = nothing
        DOF = size(X_MME, 1) - 2
    elseif _X_c > 1
        Q, DOF = residualizer(X_MME[:, 2:end])
        QQt = Q * Q'
        I_minus_Q_Qt = I - QQt
    end
    if qtl_map_algo == "idul"
        df_subs = fill(NAN, 5, 3 * n_tests + 2)
    end
    with_group = !isnothing(_args_pheno_group_file)
    @runif with_group println_to_file(string(" * A group annotation file was specified for the cis-QTL analysis."), log_file)
    calcu_variant_threshold = true
    if !with_group
        multiple_testing_method = isnothing(_args_multiple_testing) ? "acat" : _args_multiple_testing 
    else
        multiple_testing_method = isnothing(_args_multiple_testing) ? "acat2" : _args_multiple_testing 
    end
    if multiple_testing_method == "acat2"
        multiple_testing_method = "acat"
        use_acat2 = true
    else
        use_acat2 = false
    end
    if multiple_testing_method == "acat"
        use_acat2 = true
    end
    println_to_file(string(" * MultipleTesting method: ", multiple_testing_method), log_file)
    if multiple_testing_method == "acat"
        @runif with_group && !use_acat2 println_to_file(string(" * The ACAT test for each phenotype group is performed using the minimum p-values of variants across phenotypes within that group."), log_file)
        @runif with_group && !use_acat2 println_to_file(string(" * Recommend using the nested ACAT test for each phenotype group by specifying '--multiple-testing acat2' or keeping default option."), log_file)
        @runif with_group && use_acat2 println_to_file(string(" * Use nested ACAT test for each phenotype group."), log_file)
    end
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
        if multiple_testing_method == "acat" && calcu_variant_threshold
            permutation_method = nothing
            nominal_only = true
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
    if n_tests > 2
        _c2 = size(xQ, 2) + n_tests
        thread_local_storage_for_idul_multi = [Dict("X" => zeros(FloatT, _n_samples, _c2), "xtx" => zeros(FloatT2, _c2, _c2), "xty" => zeros(FloatT, _c2), "beta1" => zeros(FloatT, _c2), "r2" => zeros(FloatT, _n_samples)) for _ in 1:nthreads()] 
    end
    df_perm = DataFrame()
    df_tops = DataFrame()
    df_tops_detail = DataFrame()
    df_test_gene_annot = DataFrame()
    list_full_summary_files = String[]
    pheno_annotation.index = 1:nrow(pheno_annotation)
    chroms = intersect(string.(unique(pheno_annotation.chrom)), string.(unique(snp_annotation.chromosome)))
    if length(chroms) == 0
        error("Non-overlapping chromosomes between genotype and phenotype data!")
    end
    chunk_map = DataFrame()
    @timeit to "Prepare chunk map" if _args_low_mem
        snp_annotation.kept .= false
        pheno_annotation.kept .= false
        chunk_size = isnothing(_args_chunk_size) ? size(snp_annotation, 1) : _args_chunk_size
        for ch in chroms
            ch_snp_annotation = snp_annotation[snp_annotation.chromosome.==ch, :]
            ch_pheno_annotation = pheno_annotation[pheno_annotation.chrom.==ch, :]
            ch_pheno_annotation.start .= ch_pheno_annotation.start .- 1_000_000
            ch_pheno_annotation.end .= ch_pheno_annotation.end .+ 1_000_000 .- 1
            ranges = [ch_pheno_annotation.start[x]:ch_pheno_annotation.end[x] for x in 1:size(ch_pheno_annotation, 1)]
            result = range_vector_intersection_set(ranges, ch_snp_annotation.position)
            pheno_annotation.kept[pheno_annotation.chrom.==ch] .= result
            result = vec_in_ranges(ch_snp_annotation.position, ranges)
            snp_annotation.kept[snp_annotation.chromosome.==ch] .= result
        end
        @views genotype = genotype[:, snp_annotation.kept] 
        if length(_struct_GENO.kept_variants) != sum(snp_annotation.kept)
            _struct_GENO.kept_variants = _struct_GENO.kept_variants[snp_annotation.kept]
        end
        snp_annotation = snp_annotation[snp_annotation.kept, :]
        snp_annotation.index = 1:size(snp_annotation, 1)
        phenotype = phenotype[:, pheno_annotation.kept] 
        pheno_annotation = pheno_annotation[pheno_annotation.kept, :]
        pheno_annotation.index = 1:size(pheno_annotation, 1)
        _chid = 1
        for ch in chroms
            ch_snp_annotation = snp_annotation[(snp_annotation.chromosome.==ch).&snp_annotation.kept, :]
            ch_pheno_annotation = pheno_annotation[(pheno_annotation.chrom.==ch).&pheno_annotation.kept, :]
            ch_snp_annotation.reindex .= 1:size(ch_snp_annotation, 1)
            ch_pheno_annotation.reindex .= 1:size(ch_pheno_annotation, 1)
            ch_pheno_annotation.start .= ch_pheno_annotation.start .- 1_000_000
            ch_pheno_annotation.start[(ch_pheno_annotation.start.<ch_snp_annotation.position[1])] .= ch_snp_annotation.position[1]
            ch_pheno_annotation.end .= ch_pheno_annotation.end .+ 1_000_000 .- 1
            ch_pheno_annotation.end[(ch_pheno_annotation.end.>ch_snp_annotation.position[end])] .= ch_snp_annotation.position[end]
            batch_var_start = 1
            batch_var_end = batch_var_start
            batch_phe_start = 1
            batch_phe_end = batch_phe_start
            while (batch_phe_start <= ch_pheno_annotation.reindex[end]) & (batch_var_start <= ch_snp_annotation.reindex[end])
                logi_0 = ch_snp_annotation.position .>= ch_pheno_annotation.start[batch_phe_start]
                batch_var_start = findfirst(logi_0)
                batch_var_end = batch_var_start + chunk_size - 1
                next_var = batch_var_end + 1
                if batch_var_end > ch_snp_annotation.reindex[end]
                    batch_var_end = ch_snp_annotation.reindex[end]
                    next_var = batch_var_end
                end
                logi_1 = ch_pheno_annotation.start .<= ch_snp_annotation.position[batch_var_end]
                if next_var != batch_var_end
                    logi_2 = ch_pheno_annotation.end .< ch_snp_annotation.position[next_var]
                else
                    logi_2 = ch_pheno_annotation.end .<= ch_snp_annotation.position[next_var]
                end
                batch_phe_end = ch_pheno_annotation.reindex[findlast(logi_1 .& logi_2)] 
                logi_3 = ch_snp_annotation.position .<= ch_pheno_annotation.end[batch_phe_end]
                batch_var_end = ch_snp_annotation.reindex[findlast(logi_3)] 
                append!(chunk_map, DataFrame(chrom=ch, chunkindex=_chid, chunk_variant=[ch_snp_annotation.index[batch_var_start]:ch_snp_annotation.index[batch_var_end]], chunk_phenotype=[ch_pheno_annotation.index[batch_phe_start]:ch_pheno_annotation.index[batch_phe_end]]))
                batch_phe_start = batch_phe_end + 1
                batch_var_start = batch_var_end + 1
                _chid += 1
            end
        end
    else
        _chid = 1
        for ch in chroms
            ch_snp_annotation = snp_annotation[snp_annotation.chromosome.==ch, :]
            ch_pheno_annotation = pheno_annotation[pheno_annotation.chrom.==ch, :]
            append!(chunk_map, DataFrame(chrom=ch, chunkindex=_chid, chunk_variant=[ch_snp_annotation.index[1]:ch_snp_annotation.index[end]], chunk_phenotype=[ch_pheno_annotation.index[1]:ch_pheno_annotation.index[end]]))
            _chid += 1
        end
    end
    _task_write_full_pairs = @task @info "Using Tasks"
    chrom_mode = length(chroms) == size(chunk_map, 1)
    total_chunks = chunk_map.chunkindex[end]
    df_info = DataFrame(num_pheno=size(pheno_annotation, 1),
        num_var=size(snp_annotation, 1),
        num_perm=n_perms,
        mt_method=multiple_testing_method,
        dof=DOF,
        map_model=qtl_map_model,
    )
    if _args_run_mode == "cis_interaction"
        df_info.map_model .= replace(qtl_map_model, "a" => "a+ai")
    end
    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.info")), df_info, delim="\t")
    for chid in chunk_map.chunkindex
        chrom = chunk_map.chrom[chid]
        time_start = now()
        println_to_file(string("Chunk: ", chid, "/", total_chunks, ", Chromosome: ", chrom, ", start at: ", time_start), log_file)
        snp_index_chunk = chunk_map.chunk_variant[chid]
        _gene_annot = pheno_annotation[chunk_map.chunk_phenotype[chid], :]
        _snp_annot = snp_annotation[snp_index_chunk, :]
        _snp_annot.reindex = USE_Float32 ? Int32.(1:size(_snp_annot, 1)) : 1:size(_snp_annot, 1)
        _gene_snps = [get_cis_snp_info(_gene_annot, _snp_annot, gene, _args_cis_window; nsnp_only=true, window_type=_args_window_type) for gene in _gene_annot.pheno_id]
        _gene_annot = _gene_annot[_gene_snps.>0, :]
        _gene_snps = _gene_snps[_gene_snps.>0]
        _n_phenos = length(_gene_annot.pheno_id) 
        _gene_end_index = accumulate(+, _gene_snps)
        _gene_start_index = _gene_end_index[1:end-1] .+ 1
        prepend!(_gene_start_index, 1)
        _gene_annot.summary_start=_gene_start_index
        _gene_annot.summary_end=_gene_end_index
        append!(df_test_gene_annot, _gene_annot)
        _group_annot = nothing
        if issubset(["group_id_first"], names(_gene_annot))
            _group_annot = _gene_annot[_gene_annot.group_id_first, :]
            _group_annot.group_size .= DataFrames.combine(groupby(_gene_annot, :group_id), nrow).nrow
            select!(_group_annot, Not("pheno_id", "index"))
        end
        if with_group
            cis_genotype = nothing
            cis_genotype_iterm = nothing
            cis_genotype_2nd = nothing
            cis_g1_std = nothing
            cis_g2_std = nothing
            _group_start_index = nothing
            _group_end_index = nothing
            top_absr_exper_group = nothing
        end
        @timeit to "Prepare for nominal" begin
            nrow_df_full = sum(_gene_snps)
            df_full = DataFrame([
                "pheno_id" => repeat([""], nrow_df_full),
                "variant_id" => "",
                "start_distance" => 0,
                "af" => NAN,
                "beta_g1" => NAN,
                "beta_se_g1" => NAN,
                "pval_g1" => NaN,
            ])
            @runif with_strand insertcols!(df_full, "start_distance", "end_distance" => 0, after=true)
            _df_tops = DataFrame(chrom=_gene_annot.chrom,
                pheno_id=_gene_annot.pheno_id,
                num_var=_gene_snps,
                variant_id="",
                start_distance=0,
                af=NAN,
                beta_g1=NAN, 
                beta_se_g1=NAN,
                pval_g1=NaN,
            )
            @runif with_strand insertcols!(_df_tops, "start_distance", "end_distance" => 0, after=true)
            if with_group
                insertcols!(_df_tops, "chrom", "group_id" => _gene_annot.group_id, after=true)
                _df_tops_group = DataFrame(
                    chrom=chrom,
                    group_id=unique(_gene_annot.group_id),
                    group_size=DataFrames.combine(groupby(_gene_annot, :group_id), nrow).nrow,
                    pheno_id="",
                    num_var=0,
                    variant_id="",
                    start_distance=0,
                    af=NAN,
                    beta_g1=NAN,
                    beta_se_g1=NAN,
                    pval_g1=NaN,
                )
                @runif with_strand insertcols!(_df_tops, "start_distance", "end_distance" => 0, after=true)
            end
            chunk_genotype = Matrix{FloatT}(undef, (_n_samples, length(snp_index_chunk)))
            if !USE_GPU
                _2p = Matrix{FloatT}(undef, (1, length(_snp_annot.af)))
            else
                _snp_annot.af = _snp_annot.af |> CuArray
                _2p = CuMatrix{FloatT}(undef, (1, length(_snp_annot.af)))
            end
            _2p .= 2 .* _snp_annot.af'
            _2pq = similar(_2p)
            _2pq .= _2p .- _2p .* _snp_annot.af'
        end
        @timeit to "Prepare for 1 assoc" if n_tests == 1
            @runif calcu_variant_threshold _df_tops.pval_g1_threshold .= NaN
            @runif multiple_testing_method == "acat" _df_tops.pval_g1_pheno .= NaN
            @runif !nominal_only _df_perm = hcat(DataFrame(pheno_id=_gene_annot.pheno_id, chrom=chrom), DataFrame(Dict(Symbol(lpad(i, 5, '0')) => repeat([NAN], _n_phenos) for i in 1:(n_tests*(n_perms+1)))))
            @runif !nominal_only rename!(_df_perm, ["pheno_id", "chrom", string.("X", 1:(n_tests*(n_perms+1)))...])
            @timeit to "Prepare chunk_g" if !USE_GPU
                chunk_genotype .= genotype[:, snp_index_chunk]
                if qtl_map_model in ["d+D", "d+A", "d+A+D", "d"]
                    broadcast!(-, chunk_genotype, chunk_genotype, _2pq)
                else
                    broadcast!(-, chunk_genotype, chunk_genotype, _2p)
                end
            else
                chunk_genotype_res = genotype[:, snp_index_chunk] |> CuArray{FloatT}
                if qtl_map_model in ["d+D", "d+A", "d+A+D", "d"]
                    broadcast!(-, chunk_genotype_res, chunk_genotype_res, _2pq)
                else
                    broadcast!(-, chunk_genotype_res, chunk_genotype_res, _2p)
                end
                chunk_genotype .= chunk_genotype_res |> Array
            end
            @timeit to "Prepare chunk_g_res" if (test_geno_type == "residualized") | (perm_geno_type == "residualized")
                @timeit to "Prepare 1" if !USE_GPU
                    chunk_genotype_res = copy(chunk_genotype)
                end
                @timeit to "Prepare 2" if qtl_map_algo == "standard_fast"
                    chunk_g1_std = get_matrix_resid!(chunk_genotype_res, Q, return_std=true, center=false) 
                else
                    get_matrix_resid!(chunk_genotype_res, Q, return_std=false, center=false)
                end
            end
            @timeit to "Prepare chunk_p_res" if (test_phen_type == "residualized") | (perm_phen_type == "residualized")
                if _args_linear_model
                    chunk_phenotype_res, chunk_phenotype_std = get_matrix_resid(phenotype[:, chunk_map.chunk_phenotype[chid]], isnothing(Q) ? nothing : Array(Q), return_std=true) 
                else
                    chunk_phenotype_res = get_matrix_resid(phenotype[:, chunk_map.chunk_phenotype[chid]], Q, return_std=false) 
                end
            end
            @timeit to "Prepare rotate_g" if test_geno_type == "rotate"
                if !USE_GPU
                    chunk_SQ = glo_EA.vectors' * chunk_genotype
                else
                    chunk_SQ = cu(glo_EA.vectors') * chunk_genotype_res 
                    CUDA.unsafe_free!(chunk_genotype_res)
                end
            end
            @timeit to "Prepare rotate_p" if test_phen_type == "rotate"
                if _args_preadj_covar
                    chunk_yQ = glo_EA.vectors' * get_lm_residuals(X_MME, phenotype[:, _gene_annot.index], Xt_X)
                else
                    chunk_yQ = glo_EA.vectors' * phenotype[:, _gene_annot.index] 
                end
            end
        end
        @timeit to "Prepare for 2 assoc" if n_tests == 2
            df_full.beta_g2 .= NAN
            df_full.beta_se_g2 .= NAN
            df_full.pval_g2 .= NaN
            _df_tops.beta_g2 .= NAN
            _df_tops.beta_se_g2 .= NAN
            _df_tops.pval_g2 .= NaN
            @runif calcu_variant_threshold _df_tops.pval_g2_threshold .= NaN
            @runif multiple_testing_method == "acat" _df_tops.pval_g2_pheno .= NaN
            @runif !nominal_only _df_perm = hcat(DataFrame(pheno_id=_gene_annot.pheno_id, chrom=chrom), DataFrame(Dict(Symbol(lpad(i, 5, '0')) => repeat([NAN], _n_phenos) for i in 1:((n_tests+1)*(n_perms+1)))))
            @runif !nominal_only rename!(_df_perm, ["pheno_id", "chrom", string.("X", 1:((n_tests+1)*(n_perms+1)))...])
            if with_group
                if n_tests == 2
                    _df_tops_group.beta_g2 .= NAN
                    _df_tops_group.beta_se_g2 .= NAN
                    _df_tops_group.pval_g2 .= NaN
                end
            end
            chunk_genotype .= genotype[:, snp_index_chunk]
            chunk_genotype .-= _2p
            if _args_run_mode == "cis" 
                chunk_genotype_2nd = similar(domGenotype[:, snp_index_chunk], FloatT)
                chunk_genotype_2nd .= domGenotype[:, snp_index_chunk]
                chunk_genotype_2nd .-= _2pq
            elseif _args_run_mode == "cis_interaction"
                if qtl_map_algo != "idul"
                    chunk_genotype_2nd = copy(chunk_genotype)
                    chunk_genotype_2nd .*= _struct_ITERM.ITERM
                else
                    chunk_SQ = glo_EA.vectors' * chunk_genotype 
                    chunk_SIQ = glo_EA.vectors' .* _struct_ITERM.ITERM' * chunk_genotype 
                    if _args_preadj_covar
                        chunk_yQ = glo_EA.vectors' * get_lm_residuals(X_MME, phenotype[:, _gene_annot.index], Xt_X)
                    else
                        chunk_yQ = glo_EA.vectors' * phenotype[:, _gene_annot.index] 
                    end
                end
            end
        end
        @timeit to "Prepare for >2 assoc" if n_tests > 2
            for gi in 2:n_tests
                df_full[!,string("beta_g", gi)] .= NAN
                df_full[!,string("beta_se_g", gi)] .= NAN
                df_full[!,string("pval_g", gi)] .= NaN
                _df_tops[!,string("beta_g", gi)] .= NAN
                _df_tops[!,string("beta_se_g", gi)] .= NAN
                _df_tops[!,string("pval_g", gi)] .= NaN
                if with_group
                    if n_tests == 2
                        _df_tops_group[!,string("beta_g", gi)] .= NAN
                        _df_tops_group[!,string("beta_se_g", gi)] .= NAN
                        _df_tops_group[!,string("pval_g", gi)] .= NaN
                    end
                end
            end
            df_full.pval_joint .= NaN
            _df_tops.pval_joint .= NaN
            @runif multiple_testing_method == "acat" _df_tops.pval_joint_pheno .= NaN
            chunk_genotype .= genotype[:, snp_index_chunk]
            chunk_genotype .-= _2p
            if _args_run_mode == "cis" 
                chunk_genotype_2nd = similar(domGenotype[:, snp_index_chunk], FloatT)
                chunk_genotype_2nd .= domGenotype[:, snp_index_chunk]
                chunk_genotype_2nd .-= _2pq
            elseif _args_run_mode == "cis_interaction"
                if qtl_map_algo != "idul"
                    chunk_genotype_2nd = copy(chunk_genotype)
                    chunk_genotype_2nd .*= _struct_ITERM.ITERM
                else
                    @runif !gxe_omit_main_eff chunk_SQ = glo_EA.vectors' * chunk_genotype 
                    chunk_SIQ_list = [similar(chunk_genotype) for i in 1:n_tests-ifelse(gxe_omit_main_eff, 0, 1)]
                    @time for i in 1:n_tests-ifelse(gxe_omit_main_eff, 0, 1)
                        chunk_SIQ_list[i] .= glo_EA.vectors' .* _struct_ITERM.ITERM[:, i]' * chunk_genotype
                    end
                    GC.gc()
                    if _args_preadj_covar
                        chunk_yQ = glo_EA.vectors' * get_lm_residuals(X_MME, phenotype[:, _gene_annot.index], Xt_X)
                    else
                        chunk_yQ = glo_EA.vectors' * phenotype[:, _gene_annot.index] 
                    end
                end
            end
        end
        if with_group
            @runif calcu_variant_threshold _df_tops_group[:, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]] .= NaN
            _n_groups = size(_df_tops_group, 1)
            @runif !nominal_only _df_perm_group = hcat(DataFrame(group_id=unique(_gene_annot.group_id), chrom=chrom), DataFrame(Dict(Symbol(lpad(i, 5, '0')) => repeat([NAN], _n_groups) for i in 1:((n_tests+(1*(n_tests-1)))*(n_perms+1)))))
            @runif !nominal_only rename!(_df_perm_group, ["group_id", "chrom", string.("X", 1:((n_tests+(1*(n_tests-1)))*(n_perms+1)))...])
            _df_tops_group[:, ["pval_g1_group", "pval_g2_group"][n_tests]] .= NaN
            _df_tops_group_assigned_colnames = intersect(names(_df_tops_group), names(_df_tops))[3:end]
            _df_tops_group_assigned_indices = findall(names(_df_tops_group) .∈ (_df_tops_group_assigned_colnames,))
            _df_tops_assigned_indices = findall(names(_df_tops) .∈ (_df_tops_group_assigned_colnames,))
        end
        @runif multiple_testing_method == "beta_approx" begin
            _df_tops.pval_beta .= NaN
            _df_tops.beta_shape1 .= NAN
            _df_tops.beta_shape2 .= NAN
            _df_tops.true_dof .= NAN
            _df_tops.pval_true_dof .= NaN
        end
        @runif !nominal_only begin
            if (n_tests == 2) & (_args_run_mode == "cis_interaction") & (permutation_method == "standard") & (qtl_map_algo == "idul")
                chunk_xtx_SQt_xQ = zeros(FloatT, _c2, size(chunk_SQ, 2))
                chunk_xtx_SQ2t_xQ = zeros(FloatT, _c2, size(chunk_SQ, 2))
                chunk_xtx_SQt_xQ[1, :] .= vec(sum(abs2, chunk_SQ, dims=1))
                chunk_xtx_SQ2t_xQ[2, :] = vec(sum(abs2, chunk_SIQ, dims=1))
                chunk_xtx_SQt_xQ[2, :] .= chunk_xtx_SQ2t_xQ[1, :] = vec(sum(chunk_SQ .* chunk_SIQ, dims=1))
                chunk_xtx_SQt_xQ[3:end, :] = mat_mul(xQ', chunk_SQ)
                chunk_xtx_SQ2t_xQ[3:end, :] = mat_mul(xQ', chunk_SIQ)
                stat_perm = fill(NAN, n_perms * (n_tests + 1))
            end
            if (n_tests == 2) & (perm_geno_type == "residualized") & (_args_run_mode == "cis_interaction")
                if qtl_map_algo == "standard_fast"
                    chunk_SQ1_pcc = similar(chunk_genotype)
                    @runif _args_run_mode == "cis_interaction" get_lm_residuals!(chunk_SQ1_pcc, X_MME, chunk_genotype, chunk_genotype .* _struct_ITERM.ITERM)
                    @runif _args_run_mode == "cis" get_lm_residuals!(chunk_SQ1_pcc, X_MME, chunk_genotype, chunk_genotype_2nd)
                    chunk_SQ1_pcc ./= sqrt.(sum(abs2, chunk_SQ1_pcc, dims=1))
                    chunk_SQ1_pcc .-= Statistics.mean(chunk_SQ1_pcc, dims=1)
                end
                chunk_SQ2_pcc = similar(chunk_genotype)
                @runif _args_run_mode == "cis_interaction" get_lm_residuals!(chunk_SQ2_pcc, X_MME, chunk_genotype .* _struct_ITERM.ITERM, chunk_genotype)
                @runif _args_run_mode == "cis" get_lm_residuals!(chunk_SQ2_pcc, X_MME, chunk_genotype_2nd, chunk_genotype)
                chunk_SQ2_pcc ./= sqrt.(sum(abs2, chunk_SQ2_pcc, dims=1))
                chunk_SQ2_pcc .-= Statistics.mean(chunk_SQ2_pcc, dims=1)
                stat_perm = fill(NAN, n_perms * (n_tests + 1))
            end
        end
        GC.gc()
        for i in 1:_n_phenos
            gene = _gene_annot.pheno_id[i]
            println_to_file(string("    PHENO: ", i, "/", _n_phenos, " <", gene,">"), log_file)
            if test_phen_type == "residualized"
                if !USE_GPU
                    testpheno = chunk_phenotype_res[:, i]
                else
                    testpheno = chunk_phenotype_res[:, i] |> CuArray
                end
            end
            if test_phen_type == "rotate"
                testpheno = chunk_yQ[:, i]
            end
            pheno_index = findfirst(pheno_annotation.pheno_id .== gene)
            if test_phen_type == "raw"
                if !USE_GPU
                    testpheno = vec(phenotype[:, pheno_index]) 
                else
                    testpheno = vec(phenotype[:, pheno_index]) |> CuArray
                end
            end
            @timeit to "Keep SNPs within cis-region" cissnps_annot, n_cis_snps = get_cis_snp_info(_gene_annot, _snp_annot, gene, _args_cis_window; window_type=_args_window_type)
            if with_group 
                is_pre_data = _gene_annot.group_id_first[i]
            else
                is_pre_data = true
            end
            if !with_group | is_pre_data
                @timeit to "Pull cis-SNP genotypes" begin
                    if qtl_map_algo == "standard"
                        if !USE_GPU
                            cis_genotype = chunk_genotype[:, cissnps_annot.reindex]
                        else
                            cis_genotype = chunk_genotype[:, cissnps_annot.reindex] |> CuArray
                        end
                    elseif qtl_map_algo == "standard_fast" 
                        cis_g1_std = chunk_g1_std[cissnps_annot.reindex]
                        if _args_run_mode == "cis"
                            if !USE_GPU
                                cis_genotype = chunk_genotype_res[:, cissnps_annot.reindex]
                            else
                                cis_genotype = chunk_genotype_res[:, cissnps_annot.reindex] |> CuArray
                            end
                        elseif _args_run_mode == "cis_interaction"
                            if !USE_GPU
                                cis_genotype = chunk_SQ1_pcc[:, cissnps_annot.reindex]
                            else
                                cis_genotype = chunk_SQ1_pcc[:, cissnps_annot.reindex] |> CuArray
                            end
                        end
                        @runif _args_run_mode == "cis_interaction" cis_genotype_iterm = chunk_SQ2_pcc[:, cissnps_annot.reindex]
                        @runif _args_run_mode == "cis_interaction" cis_g2_std = chunk_g2_std[cissnps_annot.reindex]
                    elseif qtl_map_algo == "idul"
                        if n_grms < 2
                            @runif !gxe_omit_main_eff if !USE_GPU
                                cis_genotype = chunk_SQ[:, cissnps_annot.reindex] 
                            else
                                cis_genotype = chunk_SQ[:, cissnps_annot.reindex] |> Array 
                            end
                            if _args_run_mode == "cis_interaction" 
                                if n_tests > 2
                                    cis_genotype_iterm = [chunk_SIQ_list[ii][:, cissnps_annot.reindex] for ii in 1:n_tests-ifelse(gxe_omit_main_eff, 0, 1)] 
                                else
                                    cis_genotype_iterm = chunk_SIQ[:, cissnps_annot.reindex] 
                                end
                            end
                        elseif n_grms == 2 
                            cis_genotype = chunk_genotype[:, cissnps_annot.reindex] 
                        end
                    end
                    if (qtl_map_model in ["a+d+A+D", "a+d"]) & (qtl_map_algo == "idul")  
                        cis_genotype_2nd = chunk_genotype_2nd[:, cissnps_annot.reindex]
                    end
                end
            end
            if !_args_linear_model
                if (qtl_map_algo == "standard") | (n_grms == 2)
                    @timeit to "Global heritability" if !isnothing(prior_h2)
                        index_prior_h2 = findfirst(prior_h2.pheno_id .== gene)
                        Σ_i = [prior_h2.vg1[index_prior_h2], prior_h2.vg2[index_prior_h2], prior_h2.ve[index_prior_h2]]
                    else
                        if _args_h2_algo == "minmax"
                            glo_vc = getMRVCModel(glo_KS, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                        elseif _args_h2_algo == "idul"
                            glo_vc = get_IDUL_VarianceComponent(glo_EA, testpheno, X_MME, xQ=xQ)
                        end
                        Σ_i = [glo_vc[:ΣG][1], NAN, glo_vc[:Σe]]
                    end
                end
                if qtl_map_algo == "standard"
                    @timeit to "Vi" if qtl_map_model in ["a+A", "d+D", "d+A"]
                        getVinv!(Vi, Di, glo_EA, Σ_i[[1, 3]], EAvec_Di)
                    elseif qtl_map_model in ["a+A+D", "d+A+D", "a+d+A+D"]
                        getVinv!(Vi, Σ_i, _struct_PartialV, Di, EAvec_Di, EAvec)
                    end
                end
                @timeit to "Prepare for IDUL" if (qtl_map_algo == "idul") & (n_grms == 2)
                    ratio = maximum([Σ_i[1] / Σ_i[2], Σ_i[2] / Σ_i[1]])
                    @runif _args_debug println("Ratio between vg1 and vg2: ", ratio)
                    if ratio > maximum(_struct_PartialV.vec_ratio)
                        if ratio > 1000
                            if Σ_i[2] > Σ_i[1] 
                                glo_EA.vectors .= _struct_KIN.domEA.vectors
                                glo_EA.values .= _struct_KIN.domEA.values
                            else 
                                glo_EA.vectors .= _struct_KIN.EA.vectors
                                glo_EA.values .= _struct_KIN.EA.values
                            end
                        else
                            cGRM .= Σ_i[1] / sum(Σ_i[1:2]) * _struct_KIN.GRM + Σ_i[2] / sum(Σ_i[1:2]) * _struct_KIN.domGRM
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
                        @timeit to "p2" glo_EA.values .= _struct_PartialV.eigvals[vals_inds]
                    end
                    if !is_pre_data
                        cis_genotype .= chunk_genotype[:, cissnps_annot.reindex]
                        cis_genotype_2nd .= chunk_genotype_2nd[:, cissnps_annot.reindex]
                    end
                    @timeit to "p3" cis_genotype .= glo_EA.vectors' * cis_genotype
                    @runif qtl_map_model == "a+d+A+D" cis_genotype_2nd .= glo_EA.vectors' * cis_genotype_2nd 
                    @timeit to "p4" mul!(yQ, glo_EA.vectors', testpheno)
                    @timeit to "p5" mul!(xQ, glo_EA.vectors', X_MME)
                end
            end
            @timeit to "Nominal association" if _args_linear_model
                @timeit to "Get associations" if _args_run_mode == "cis"
                    if !is_pre_data 
                        cis_genotype .= chunk_genotype_res[:, cissnps_annot.reindex]
                    end
                    if !USE_GPU
                        df_test = fill(NAN, n_cis_snps, 3)
                        assoc_r = fill(NAN, n_cis_snps)
                    else
                        df_test = CUDA.fill(NAN, n_cis_snps, 3)
                        assoc_r = CUDA.fill(NAN, n_cis_snps)
                    end
                    mul!(assoc_r, cis_genotype', testpheno)
                    get_approx_slope_se_from_r!(df_test, assoc_r, DOF, cis_g1_std, chunk_phenotype_std[i]; is_calcu_pv=false)
                elseif _args_run_mode == "cis_interaction"
                    if qtl_map_algo == "idul"
                        if !is_pre_data 
                            cis_genotype .= chunk_SQ[:, cissnps_annot.reindex]
                            if n_tests > 2
                                [cis_genotype_iterm[ii] .= chunk_SIQ_list[ii][:, cissnps_annot.reindex] for ii in 1:n_tests-1] 
                            else
                                cis_genotype_iterm .= chunk_SIQ[:, cissnps_annot.reindex] 
                            end
                        end
                        df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                        if n_tests > 2
                            idul_multi_assoc_test_0eta!(df_test, cis_genotype, testpheno, xQ, cis_genotype_iterm; is_calcu_pv=false)
                        else
                            idul_two_assoc_test_0eta!(df_test, cis_genotype, testpheno, xQ, cis_genotype_iterm; is_calcu_pv=false)
                        end
                    elseif qtl_map_algo == "standard_fast"
                        if !is_pre_data 
                            cis_genotype .= chunk_SQ1_pcc[:, cissnps_annot.reindex]
                            cis_genotype_iterm .= chunk_SQ2_pcc[:, cissnps_annot.reindex]
                        end
                        df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                        assoc_r = fill(NAN, n_cis_snps)
                        yadj = get_phenotype_resid(testpheno, I_minus_Q_Qt) 
                        mul!(assoc_r, cis_genotype', yadj)
                        get_approx_slope_se_from_r!(view(df_test, :, 1:3), assoc_r, DOF, cis_g1_std, std(testpheno); is_calcu_pv=false)
                        mul!(assoc_r, cis_genotype_iterm', yadj)
                        get_approx_slope_se_from_r!(view(df_test, :, 4:6), assoc_r, DOF, cis_g2_std, std(testpheno); is_calcu_pv=false)
                    end
                end
            else
                if qtl_map_algo == "idul"
                    if _args_run_mode == "cis"
                        if n_tests == 1
                            if !is_pre_data 
                                @runif test_geno_type == "rotate" @runif !gxe_omit_main_eff cis_genotype .= chunk_SQ[:, cissnps_annot.reindex]
                            end
                            @runif qtl_map_model == "d+A+D" testpheno .= yQ
                            df_test = fill(NAN, n_cis_snps, 3)
                            Random.seed!(_args_seed + n_cis_snps)
                            idul_prior_index = sample(1:n_cis_snps, min(5, n_cis_snps); replace=false)
                            fill!(df_subs, NAN)
                            @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, cis_genotype[:, idul_prior_index], testpheno, xQ, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=FloatT(1.5), return_detail=true, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision))
                            init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                            @runif _args_verbose println(string("    η: ",init_eta))
                            @runif isnan(init_eta) begin
                                init_eta = FloatT(1e-5)
                                println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                                @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                            end
                            if !_args_exact_map
                                @timeit to "Idul test" idul_assoc_test_approx!(df_test, cis_genotype, testpheno, xQ, glo_EA.values; init_eta=init_eta, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision)) 
                            else
                                @timeit to "Idul test" idul_assoc_test_plus!(df_test, cis_genotype, testpheno, xQ, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                            end
                        elseif n_tests == 2
                            df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                            Random.seed!(_args_seed + n_cis_snps)
                            idul_prior_index = sample(1:n_cis_snps, min(5, n_cis_snps); replace=false)
                            fill!(df_subs, NAN)
                            @timeit to "Prior idul" idul_two_assoc_test!(df_subs, cis_genotype[:, idul_prior_index], testpheno, xQ, glo_EA.values, cis_genotype_2nd[:, idul_prior_index]; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=FloatT(1.5), return_detail=true, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision))
                            init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                            @runif _args_verbose println(string("    η: ",init_eta))
                            @runif isnan(init_eta) begin
                                init_eta = FloatT(1e-5)
                                println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                                @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                            end
                            if !_args_exact_map
                                @timeit to "Idul test" idul_two_assoc_test_approx!(df_test, cis_genotype, testpheno, xQ, glo_EA.values, cis_genotype_2nd; init_eta=init_eta, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision)) 
                            else
                                idul_two_assoc_test!(df_test, cis_genotype, testpheno, xQ, glo_EA.values, cis_genotype_2nd; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                            end
                        end
                    elseif _args_run_mode == "cis_interaction"
                        if !is_pre_data 
                            @runif !gxe_omit_main_eff cis_genotype .= chunk_SQ[:, cissnps_annot.reindex]
                            if n_tests > 2
                                [cis_genotype_iterm[ii] .= chunk_SIQ_list[ii][:, cissnps_annot.reindex] for ii in 1:n_tests-ifelse(gxe_omit_main_eff, 0, 1)]
                            else
                                cis_genotype_iterm .= chunk_SIQ[:, cissnps_annot.reindex]
                            end
                        end
                        df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                        Random.seed!(_args_seed + n_cis_snps)
                        idul_prior_index = sample(1:n_cis_snps, min(5, n_cis_snps); replace=false)
                        fill!(df_subs, NAN)
                        if n_tests > 2
                            if gxe_omit_main_eff
                                @timeit to "Prior idul" idul_multi_assoc_test!(df_subs, testpheno, xQ, glo_EA.values, [cis_genotype_iterm[ii][:, idul_prior_index] for ii in 1:n_tests]; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=FloatT(1.5), return_detail=true, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision))
                            else
                                @timeit to "Prior idul" idul_multi_assoc_test!(df_subs, cis_genotype[:, idul_prior_index], testpheno, xQ, glo_EA.values, [cis_genotype_iterm[ii][:, idul_prior_index] for ii in 1:n_tests-1]; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=FloatT(1.5), return_detail=true, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision))
                            end
                        else
                            @timeit to "Prior idul" idul_two_assoc_test!(df_subs, cis_genotype[:, idul_prior_index], testpheno, xQ, glo_EA.values, cis_genotype_iterm[:, idul_prior_index]; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=FloatT(1.5), return_detail=true, is_calcu_pv=false, add_to_diag=FloatT(_args_inv_precision))    
                        end
                        init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                        @runif _args_verbose println(string("    η: ",init_eta))
                        @runif isnan(init_eta) begin
                            init_eta = FloatT(1e-5)
                            println_to_file(string(" * [INFO] The XtX matrix is non-positive definite for ", gene, ". Linear model will be used instead."), log_file)
                            @runif _args_inv_precision == 1e-10 println_to_file(string(" * [WARN] A small value (e.g., 1e-5) can be added to the diagonal by specify the option '--inv-precision'."), log_file)
                        end
                        if !_args_exact_map
                            if n_tests > 2
                                if gxe_omit_main_eff
                                    @timeit to "Idul test" idul_multi_assoc_test_approx!(df_test, testpheno, xQ, glo_EA.values, cis_genotype_iterm, thread_local_storage_for_idul_multi; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                                else
                                    @timeit to "Idul test" idul_multi_assoc_test_approx!(df_test, cis_genotype, testpheno, xQ, glo_EA.values, cis_genotype_iterm; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                                end
                            else
                                @timeit to "Idul test" idul_two_assoc_test_approx!(df_test, cis_genotype, testpheno, xQ, glo_EA.values, cis_genotype_iterm; init_eta=init_eta, is_calcu_pv=false,add_to_diag=FloatT(_args_inv_precision)) 
                            end
                        else
                            idul_two_assoc_test!(df_test, cis_genotype, testpheno, xQ, glo_EA.values, cis_genotype_iterm; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), init_eta=init_eta, is_calcu_pv=false)
                        end
                    end
                elseif qtl_map_algo == "standard"
                    @timeit to "Adjusted phenotype" yadj = testpheno 
                    @timeit to "Get associations" if _args_run_mode == "cis"
                        if qtl_map_model in ["a+d+A+D"]
                            df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                            calcu_two_assoc_test!(df_test, Vi, cis_genotype, cis_genotype_2nd, yadj)
                        else
                            df_test = fill(NAN, n_cis_snps, 3)
                            calcu_one_assoc_test!(df_test, Vi, cis_genotype, yadj)
                        end
                    elseif _args_run_mode == "cis_interaction"
                        df_test = fill(NAN, n_cis_snps, 3 * n_tests + 1)
                        calcu_two_assoc_test!(df_test, Vi, cis_genotype, _struct_ITERM.ITERM, yadj)
                    end
                end
            end
            @timeit to "Temporary DF to store associations" begin
                _df_full_index = _gene_start_index[i]:_gene_end_index[i]
                df_full.pheno_id[_df_full_index] .= gene
                df_full.variant_id[_df_full_index] .= cissnps_annot.variant
                df_full.start_distance[_df_full_index] .= cissnps_annot.start_distance
                @runif with_strand df_full.end_distance[_df_full_index] .= cissnps_annot.end_distance
                if USE_GPU
                    df_full.af[_df_full_index] .= cissnps_annot.af |> Array
                else
                    df_full.af[_df_full_index] .= cissnps_annot.af
                end
                @runif USE_GPU df_test = df_test |> Array
                if n_tests == 1
                    df_full.beta_g1[_df_full_index] .= df_test[:, 1]
                    df_full.beta_se_g1[_df_full_index] .= df_test[:, 2]
                    if _args_linear_model & (qtl_map_algo != "idul") 
                        df_full.pval_g1[_df_full_index] .= cdf(TDist(DOF), -df_test[:, 3]) * 2
                    else
                        df_full.pval_g1[_df_full_index] .= ccdf(WaldTest(1, DOF), df_test[:, 3])
                    end
                end
                if n_tests == 2
                    df_full.beta_g1[_df_full_index] .= df_test[:, 1]
                    df_full.beta_se_g1[_df_full_index] .= df_test[:, 2]
                    if _args_linear_model & (qtl_map_algo != "idul") 
                        df_full.pval_g1[_df_full_index] .= ccdf(TDist(DOF - 1), df_test[:, 3]) * 2
                    else
                        df_full.pval_g1[_df_full_index] .= ccdf(WaldTest(1, DOF - 1), df_test[:, 3])
                    end
                    df_full.beta_g2[_df_full_index] .= df_test[:, 4]
                    df_full.beta_se_g2[_df_full_index] .= df_test[:, 5]
                    if _args_linear_model & (qtl_map_algo != "idul") 
                        df_full.pval_g2[_df_full_index] .= ccdf(TDist(DOF - 1), df_test[:, 6]) * 2
                    else
                        df_full.pval_g2[_df_full_index] .= ccdf(WaldTest(1, DOF - 1), df_test[:, 6])
                    end
                end
                if n_tests > 2
                    pval_cols_indices = findall(names(df_full) .∈ (string.("pval_g", 1:n_tests),))
                    beta_cols_indices = findall(names(df_full) .∈ (string.("beta_g", 1:n_tests),))
                    beta_se_cols_indices = findall(names(df_full) .∈ (string.("beta_se_g", 1:n_tests),))
                    @timeit to "1" df_full[_df_full_index, beta_cols_indices] .= df_test[:, 1:3:n_tests*3]
                    @timeit to "2" df_full[_df_full_index, beta_se_cols_indices] .= df_test[:, 2:3:n_tests*3]
                    @timeit to "3" if _args_linear_model & (qtl_map_algo != "idul") 
                        df_full[_df_full_index, pval_cols_indices] .= ccdf(TDist(DOF - 1), df_test[:, 3:3:n_tests*3]) * 2
                    else
                        df_full[_df_full_index, pval_cols_indices] .= ccdf(WaldTest(1, DOF - 1), df_test[:, 3:3:n_tests*3])
                    end
                    @timeit to "4" @inbounds @threads for bb in eachindex(_df_full_index)
                        df_full.pval_joint[_df_full_index[bb]] = ACATest(df_full[_df_full_index[bb], pval_cols_indices] |> Vector; is_check=false)
                    end
                end
            end
            top_absr_exper = append_top_summary!(df_full[_df_full_index, :], _df_tops, i, n_tests, DOF) 
            if !with_group
                @runif !nominal_only @timeit to "Permutation" if n_tests == 1 
                    @timeit to "Perm 1.0" begin
                        Ainds = CartesianIndices(cis_genotype)
                        Binds = CartesianIndices((1:size(chunk_genotype_res, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                        copyto!(cis_genotype, Ainds, chunk_genotype_res, Binds)
                    end
                    @timeit to "Perm clipper" if permutation_method == "clipper"
                        permpheno = chunk_phenotype_res[:, i]
                        pval_nominal_threshold = cis_permutation_one_assoc!("clipper", permpheno, nothing, cis_genotype, absr_perm, DOF, fdr=FloatT(_args_fdr))
                        _df_perm[i, 3:end] .= [top_absr_exper; absr_perm]
                    end
                    @timeit to "Perm standard_fast" if permutation_method == "standard_fast"
                        permpheno = phenotype[:, pheno_index]
                        pval_nominal_threshold = cis_permutation_one_assoc!("standard_fast", permpheno, Q, cis_genotype, absr_perm, DOF; y_perms=y_perms, fdr=FloatT(_args_fdr))
                        _df_perm[i, 3:end] .= [top_absr_exper; absr_perm]
                    end
                    _df_tops.pval_g1_threshold[i] = pval_nominal_threshold
                elseif (n_tests == 2) & (_args_run_mode == "cis")
                    if permutation_method == "standard"
                        permpheno = phenotype[:, pheno_index]
                        cis_permutation_two_assoc!("standard", i, permpheno, cis_genotype, cis_genotype_2nd, _df_perm, _df_tops, df_full[_df_full_index, :], y_perms, cissnps_annot; X_MME=nothing, Xt_X=nothing, glo_EA=glo_EA, chunk_xtx_SQt_xQ=chunk_xtx_SQt_xQ, chunk_xtx_SQ2t_xQ=chunk_xtx_SQ2t_xQ, xW0=xW0, xty_init=xty_init, xtx_init=xtx_init, fdr=FloatT(_args_fdr))
                    elseif permutation_method == "standard_fast"
                        Ainds = CartesianIndices(cis_genotype_2nd)
                        Binds = CartesianIndices((1:size(chunk_SQ2_pcc, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                        copyto!(cis_genotype_2nd, Ainds, chunk_SQ2_pcc, Binds)
                        permpheno = phenotype[:, pheno_index] 
                        cis_permutation_two_assoc!("standard_fast", i, permpheno, cis_genotype, cis_genotype_2nd, _df_perm, _df_tops, df_full[_df_full_index, :], y_perms, cissnps_annot, stat_perm; X_MME=X_MME, Xt_X=Xt_X, fdr=FloatT(_args_fdr))
                    end
                elseif (n_tests == 2) & (_args_run_mode == "cis_interaction")
                    if permutation_method == "standard"
                        permpheno = phenotype[:, pheno_index]
                        pval_nominal_threshold = cis_permutation_two_assoc!("standard", permpheno, cis_genotype, cis_genotype_iterm, y_perms, cissnps_annot, stat_perm; X_MME=nothing, Xt_X=nothing, glo_EA=glo_EA, chunk_xtx_SQt_xQ=chunk_xtx_SQt_xQ, chunk_xtx_SQ2t_xQ=chunk_xtx_SQ2t_xQ, xW0=xW0, xty_init=xty_init, xtx_init=xtx_init, fdr=FloatT(_args_fdr))
                        _df_perm[i, 3:end] .= [NAN; top_absr_exper; NAN; stat_perm]
                    elseif permutation_method == "standard_fast"
                        Ainds = CartesianIndices(cis_genotype_iterm)
                        Binds = CartesianIndices((1:size(chunk_SQ2_pcc, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                        copyto!(cis_genotype_iterm, Ainds, chunk_SQ2_pcc, Binds)
                        permpheno = phenotype[:, pheno_index] 
                        pval_nominal_threshold = cis_permutation_two_assoc!("standard_fast", permpheno, cis_genotype, cis_genotype_iterm, y_perms, cissnps_annot, stat_perm; X_MME=X_MME, Xt_X=Xt_X, fdr=FloatT(_args_fdr))
                        _df_perm[i, 3:end] .= [NAN; top_absr_exper; NAN; stat_perm]
                    end
                    _df_tops.pval_g2_threshold[i] = pval_nominal_threshold
                else
                    error("!")
                end
            end
            @timeit to "ACATest" if multiple_testing_method == "acat"
                if n_tests > 2
                    pvals_for_acat = df_full[_df_full_index, "pval_joint"]
                else
                    pvals_for_acat = df_full[_df_full_index, ["pval_g1", "pval_g2"][n_tests]]
                end
                if _args_dtss_weight
                    exp_weight_val = 5e-6
                    _df_tops[i, ["pval_g1_pheno", "pval_g2_pheno"][n_tests]] = ACATest(pvals_for_acat; Weights=exp.(.-abs.(exp_weight_val * df_full.start_distance[_df_full_index])), is_check=false)
                else
                    if n_tests > 2
                        _df_tops.pval_joint_pheno[i] = ACATest(pvals_for_acat; is_check=false)
                    else
                        _df_tops[i, ["pval_g1_pheno", "pval_g2_pheno"][n_tests]] = ACATest(pvals_for_acat; is_check=false)
                    end
                end
            end
            if multiple_testing_method == "beta_approx"
                append_top_summary!(df_full[_df_full_index, :], _df_tops, i, n_tests)
                if n_tests == 1
                    r2_nominal = abs2(_df_perm.X1[i])
                elseif n_tests == 2
                    r2_nominal = abs2(_df_perm.X2[i])
                    absr_perm = _df_perm[i, ((2-1)+6):3:size(_df_perm, 2)] |> Vector
                end
                _df_tops[i, ["pval_beta", "beta_shape1", "beta_shape2", "true_dof", "pval_true_dof"]] .= calculate_beta_approx_pval(abs2.(absr_perm), r2_nominal, DOF, 1e-4)
            end
            if with_group
                if _gene_annot.group_id_first[i]
                    top_absr_exper_group = top_absr_exper
                    _group_start_index = _gene_start_index[i]
                else
                    if top_absr_exper > top_absr_exper_group
                        top_absr_exper_group = top_absr_exper
                    end
                end
                if _gene_annot.group_id_last[i]
                    gid = _gene_annot.group_id[i]
                    _group_i = findfirst(_df_tops_group.group_id .== gid)
                    _group_end_index = _gene_end_index[i]
                    group_size = _group_annot.group_size[findfirst(_group_annot.group_id .== _gene_annot.group_id[i])]
                    group_pheno_indices = _gene_annot.index[_gene_annot.group_id.==gid]
                    _pval_group = reshape(df_full[:, ["pval_g1", "pval_g2"][n_tests]][_group_start_index:_group_end_index], :, group_size)
                    best_grpi = argmin(_df_tops[_df_tops.group_id.==gid, ["pval_g1", "pval_g2"][n_tests]])
                    _df_tops_group[_group_i, _df_tops_group_assigned_indices] = _df_tops[_df_tops.group_id.==gid, _df_tops_assigned_indices][best_grpi, :]
                    @runif !nominal_only begin
                        _pval_group_min_indices = [argmin(x) for x in eachrow(_pval_group)]
                        perm_r_mat = fill(NAN, n_perms, length(_pval_group_min_indices))
                        for _grpi in 1:group_size
                            pheno_index = group_pheno_indices[_grpi]
                            variant_subset_indices = findall(_pval_group_min_indices .== _grpi)
                            @runif !nominal_only @timeit to "Permutation" if n_tests == 1
                                @timeit to "Perm 1.0" begin 
                                    Ainds = CartesianIndices(cis_genotype)
                                    Binds = CartesianIndices((1:size(chunk_genotype_res, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                                    copyto!(cis_genotype, Ainds, chunk_genotype_res, Binds)
                                end
                                @timeit to "Perm clipper" if permutation_method == "clipper"
                                    permpheno = chunk_phenotype_res[:, i]
                                    perm_r_mat[:, variant_subset_indices] = cis_permutation_one_assoc!("clipper", permpheno, cis_genotype, absr_perm, variant_subset_indices; Q=nothing, y_perms=nothing)
                                end
                                @timeit to "Perm standard_fast" if permutation_method == "standard_fast"
                                    permpheno = phenotype[:, pheno_index]
                                    perm_r_mat[:, variant_subset_indices] = cis_permutation_one_assoc!("standard_fast", permpheno, cis_genotype, absr_perm, variant_subset_indices; Q=Q, y_perms=y_perms)
                                end
                            elseif (n_tests == 2) & (_args_run_mode == "cis")
                                if permutation_method == "standard"
                                    permpheno = phenotype[:, pheno_index]
                                    cis_permutation_two_assoc!("standard", i, permpheno, cis_genotype, cis_genotype_2nd, _df_perm, _df_tops, df_full[_df_full_index, :], y_perms, cissnps_annot;
                                        X_MME=nothing, Xt_X=nothing, glo_EA=glo_EA, chunk_xtx_SQt_xQ=chunk_xtx_SQt_xQ, chunk_xtx_SQ2t_xQ=chunk_xtx_SQ2t_xQ, xW0=xW0, xty_init=xty_init,
                                        xtx_init=xtx_init, fdr=FloatT(_args_fdr))
                                elseif permutation_method == "standard_fast"
                                    Ainds = CartesianIndices(cis_genotype_2nd)
                                    Binds = CartesianIndices((1:size(chunk_SQ2_pcc, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                                    copyto!(cis_genotype_2nd, Ainds, chunk_SQ2_pcc, Binds)
                                    permpheno = phenotype[:, pheno_index] 
                                    cis_permutation_two_assoc!("standard_fast", i, permpheno, cis_genotype, cis_genotype_2nd, _df_perm, _df_tops, df_full[_df_full_index, :], y_perms, cissnps_annot, stat_perm; X_MME=X_MME, Xt_X=Xt_X, fdr=FloatT(_args_fdr))
                                end
                            elseif (n_tests == 2) & (_args_run_mode == "cis_interaction")
                                if permutation_method == "standard"
                                    permpheno = phenotype[:, pheno_index]
                                    cis_permutation_two_assoc!("standard", i, permpheno, cis_genotype, cis_genotype_iterm, _df_perm, _df_tops, df_full[_df_full_index, :], y_perms, cissnps_annot;
                                        X_MME=nothing, Xt_X=nothing, glo_EA=glo_EA, chunk_xtx_SQt_xQ=chunk_xtx_SQt_xQ, chunk_xtx_SQ2t_xQ=chunk_xtx_SQ2t_xQ, xW0=xW0, xty_init=xty_init,
                                        xtx_init=xtx_init, fdr=FloatT(_args_fdr))
                                elseif permutation_method == "standard_fast"
                                    Ainds = CartesianIndices(cis_genotype_iterm)
                                    Binds = CartesianIndices((1:size(chunk_SQ2_pcc, 1), minimum(cissnps_annot.reindex):maximum(cissnps_annot.reindex)))
                                    copyto!(cis_genotype_iterm, Ainds, chunk_SQ2_pcc, Binds)
                                    permpheno = phenotype[:, pheno_index] 
                                    perm_r_mat[:, variant_subset_indices] = cis_permutation_two_assoc!("standard_fast", permpheno, cis_genotype_iterm, y_perms, X_MME, Xt_X, variant_subset_indices)
                                end
                            else
                                error("!")
                            end
                        end
                        if n_tests == 1
                            absr_perm .= vec(maximum(abs, perm_r_mat, dims=2))
                            _df_perm_group[_group_i, 3:end] .= [top_absr_exper_group; absr_perm] 
                            pval_nominal_threshold = get_approx_p_from_r(sort!(absr_perm; rev=true)[ceil(Int, n_perms * _args_fdr)], DOF)
                        elseif n_tests == 2
                            stat_perm[2:3:length(stat_perm)] .= vec(maximum(abs, perm_r_mat, dims=2))
                            _df_perm_group[_group_i, 3:end] .= [NAN; top_absr_exper_group; NAN; stat_perm]
                            pval_nominal_threshold = get_approx_p_from_r(sort(stat_perm[2:3:length(stat_perm)]; rev=true)[ceil(Int, n_perms * _args_fdr)], DOF - 1)
                        end
                        top_absr_exper_group = nothing
                        _df_tops_group[_group_i, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]] = pval_nominal_threshold
                    end
                    if multiple_testing_method == "acat"
                        if !use_acat2
                            best_p_group_for_acat = vec(minimum(_pval_group, dims=2)) 
                            _df_tops_group[_group_i, ["pval_g1_group", "pval_g2_group"][n_tests]] = ACATest(best_p_group_for_acat; is_check=false)
                        elseif use_acat2
                            if group_size == 1
                                _df_tops_group[_group_i, ["pval_g1_group", "pval_g2_group"][n_tests]] = ACATest(_pval_group[:,1]; is_check=false)
                            else
                                _df_tops_group[_group_i, ["pval_g1_group", "pval_g2_group"][n_tests]] = ACATest(_pval_group; dims=1, is_check=false)
                            end
                        end
                    end
                end
            end
            @runif _args_debug if i % 10 == 0
                println_to_file(string("  *** Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
                println_to_file(string(to), log_file)
                println()
            end
            @runif USE_GPU & _args_linear_model CUDA.unsafe_free!(cis_genotype)
            @runif USE_GPU & _args_linear_model CUDA.unsafe_free!(cis_g1_std)
        end 
        @runif USE_GPU & _args_linear_model CUDA.unsafe_free!(chunk_genotype_res)
        @runif USE_GPU CUDA.unsafe_free!(_2p)
        @runif USE_GPU CUDA.unsafe_free!(_2pq)
        @runif USE_GPU & _args_linear_model CUDA.unsafe_free!(chunk_g1_std)
        @runif !nominal_only | (multiple_testing_method == "acat") begin
            if with_group
                append!(df_tops, _df_tops_group)
                append!(df_tops_detail, _df_tops)
                @runif !nominal_only append!(df_perm, _df_perm_group)
            else
                append!(df_tops, _df_tops)
                @runif !nominal_only append!(df_perm, _df_perm)
            end
        end
        @timeit to "Write the full summary to disk" if !_args_write_top 
            if chid != chunk_map.chunkindex[1]
                wait(_task_write_full_pairs) 
            end
            if nrow(df_full) > 0
                task_df_full = copy(df_full)
                task_chunk = chrom_mode ? chrom : string("chunk", chid, ".", chrom)
                _task_write_full_pairs = @spawn begin
                    task_cols_float = [task_df_full[1, x] isa AbstractFloat for x in range(1, ncol(task_df_full))]
                    task_df_full[:, task_cols_float] .= round.(task_df_full[:, task_cols_float], sigdigits=6)
                    task_df_full.af .= round.(task_df_full.af, digits=3)
                    if isnothing(_args_output_format)
                        chunk_full_out_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl_pairs.", task_chunk, ".txt.gz"))
                        push!(list_full_summary_files, chunk_full_out_file)
                        @runif _args_use_gzip chunk_full_out_file = replace(chunk_full_out_file, r".gz$" => "")
                        CSV.write(chunk_full_out_file, task_df_full, delim="\t", compress=!_args_use_gzip)
                        @runif _args_use_gzip run(`gzip -f $chunk_full_out_file`; wait=chrom == chroms[end])
                    else
                        pairs_file_suffix = _args_output_format in ["jld", "jld_compress"] ? ".jld" : _args_output_format in ["jld2", "jld2_compress"] ? ".jld2" : _args_output_format in ["arrow", "arrow_zstd" , "arrow_lz4"] ? ".arrow" : ".txt.gz"
                        chunk_full_out_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl_pairs.", task_chunk, pairs_file_suffix))
                        push!(list_full_summary_files, chunk_full_out_file)
                        write_qtl_pairs_file(chunk_full_out_file, task_df_full; format=_args_output_format)
                    end
                end
            end
            if chid == chunk_map.chunkindex[end]
                wait(_task_write_full_pairs)
            end
        end
        @runif !nominal_only | (multiple_testing_method == "acat") begin
            cols_float = [_df_tops[1, x] isa AbstractFloat for x in range(1, ncol(_df_tops))]
            _df_tops[:, cols_float] .= round.(_df_tops[:, cols_float], sigdigits=6)
            if with_group
                cols_float = [_df_tops_group[1, x] isa AbstractFloat for x in range(1, ncol(_df_tops_group))]
                _df_tops_group[:, cols_float] .= round.(_df_tops_group[:, cols_float], sigdigits=6)
            end
            if length(_args_chrom) > 0
                task_chunk = chrom_mode ? chrom : string("chunk", chid, ".", chrom)
                is_append = false
                cis_qtl_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.", task_chunk, ".txt.gz"))
                @runif with_group cis_qtl_detail_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.detail.", task_chunk, ".txt.gz"))
                @runif !nominal_only cis_perm_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.perm.", task_chunk, ".arrow")) 
            else
                is_append = chid != chunk_map.chunkindex[1]
                cis_qtl_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.txt.gz"))
                @runif with_group cis_qtl_detail_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.detail.txt.gz"))
                @runif !nominal_only cis_perm_file = joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.perm.arrow"))
            end
            if with_group
                if multiple_testing_method != "acat"
                    _df_tops_detail = select(_df_tops, Not(["pval_g1_threshold", "pval_g2_threshold"][n_tests]))
                else
                    _df_tops_detail = _df_tops
                end
                CSV.write(cis_qtl_file, _df_tops_group, delim="\t", compress=endswith(cis_qtl_file, ".gz"), append=is_append)
                CSV.write(cis_qtl_detail_file, _df_tops_detail, delim="\t", compress=endswith(cis_qtl_detail_file, ".gz"), append=is_append)
                @runif !nominal_only if !is_append
                    Arrow.write(cis_perm_file, _df_perm_group, compress=:zstd, file=false)
                else
                    Arrow.append(cis_perm_file, _df_perm_group)
                end
            else
                CSV.write(cis_qtl_file, _df_tops, delim="\t", compress=endswith(cis_qtl_file, ".gz"), append=is_append)
                @runif !nominal_only if !is_append
                    Arrow.write(cis_perm_file, _df_perm, compress=:zstd, file=false)
                else
                    Arrow.append(cis_perm_file, _df_perm)
                end
            end
        end
        println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
        if _args_debug
            println_to_file(string(to), log_file)
            println()
        end
    end
    @runif multiple_testing_method == "beta_approx" if length(_args_chrom) == 0
        mt_keep_index = .!isnan.(df_tops.pval_g1)
        df_tops_mt = df_tops[mt_keep_index, :]
        pval_mt = df_tops_mt.pval_beta
        qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
        if sum(qval_mt .<= _args_fdr) > 0
            set0_indices = findall(qval_mt .<= _args_fdr)
            set1_indices = findall(qval_mt .> _args_fdr)
            pthreshold = (sort(df_tops_mt.pval_beta[set1_indices])[1] - sort(-1.0 * df_tops_mt.pval_beta[set0_indices])[1]) / 2
            @runif n_tests == 1 df_tops_mt.pval_g1_threshold .= quantile.(Beta.(df_tops_mt.beta_shape1, df_tops_mt.beta_shape2), pthreshold)
            @runif n_tests == 2 df_tops_mt.pval_g2_threshold .= quantile.(Beta.(df_tops_mt.beta_shape1, df_tops_mt.beta_shape2), pthreshold)
        end
        df_tops_mt[:, ["qval_g1", "qval_g2"]] .= qval_mt
        cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
        df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.txt.gz")), df_tops_mt, delim="\t", compress=true)
        n_significant_phenotypes = sum(qval_mt .< FloatT(_args_fdr))
        logtxt = string("  *** ", n_significant_phenotypes, " out of ", length(qval_mt), " tested phenotypes with qval_g", n_tests, "<", _args_fdr, ".")
        println_to_file(logtxt, log_file)
    end
    @runif multiple_testing_method == "acat" if length(_args_chrom) == 0
        mt_col = n_tests > 2 ? "joint" : string("g", n_tests)
        pval_colname = string("pval_", mt_col)
        mt_keep_index = .!isnan.(df_tops[:, with_group ? string("pval_", mt_col, "_group") : string("pval_", mt_col, "_pheno")])
        df_tops_keep = df_tops[mt_keep_index, :]
        df_tops_mt = df_tops_keep
        pval_mt = df_tops_mt[:, with_group ? string("pval_", mt_col, "_group") : string("pval_", mt_col, "_pheno")]
        qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochberg())
        df_tops_mt[:, string("qval_", mt_col)] .= qval_mt
        @timeit to "Get final acat threshold" begin
            println_to_file(string("    Calculating the ", string("pval_", mt_col, "_threshold"), " ..."), log_file)
            fdr_closest_i = find_closest_below_fast(qval_mt, _args_fdr)
            @timeit to "Get 1" if fdr_closest_i != -1
                fdr_closest_i_up = find_closest_up_fast(qval_mt, _args_fdr)
                if with_group
                    critical_acat_p = sqrt(df_tops_mt[fdr_closest_i, string("pval_", mt_col, "_group")] * df_tops_mt[fdr_closest_i_up, string("pval_", mt_col, "_group")])
                else
                    critical_acat_p = sqrt(df_tops_mt[fdr_closest_i, string("pval_", mt_col, "_pheno")] * df_tops_mt[fdr_closest_i_up, string("pval_", mt_col, "_pheno")])
                end
                full_summary_files = list_full_summary_files
                for fn in full_summary_files
                    println_to_file(string("    *** Processing ", fn), log_file)
                    _cis_qtl_pairs = read_qtl_pairs_file(fn)
                    pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                    group_ids = nothing
                    println_to_file(string("        Number of phenotypes: ", length(pheno_ids)), log_file)
                    if with_group 
                        group_ids = unique(df_tops_detail.group_id[df_tops_detail.pheno_id.∈(pheno_ids,)])
                        println_to_file(string("        Number of groups: ", length(group_ids)), log_file)
                    end
                    @threads for _pid in pheno_ids
                        _pi = findfirst(df_test_gene_annot.pheno_id .== _pid)
                        _df_full_index = df_test_gene_annot.summary_start[_pi]:df_test_gene_annot.summary_end[_pi]
                        _pval_mt = _cis_qtl_pairs[_df_full_index, pval_colname]
                        nonan_indices = .!isnan.(_pval_mt)
                        if any(.!nonan_indices)
                            _pval_mt = _pval_mt[nonan_indices]
                        end
                        if length(_pval_mt) > 0
                            _, pval_nominal_threshold, _ = calcu_target_pval_threshold_acat!(_pval_mt, critical_acat_p)
                            if !with_group 
                                df_tops_mt[findfirst(df_tops_mt.pheno_id.==_pid), string("pval_", mt_col, "_threshold")] = pval_nominal_threshold
                            else
                                df_tops_detail[findfirst(df_tops_detail.pheno_id.==_pid), string("pval_", mt_col, "_threshold")] = pval_nominal_threshold
                            end
                        end
                    end
                    if with_group
                        matchidx = vmatch(df_tops_detail.pheno_id, df_tops_mt.pheno_id)
                        df_tops_mt[.!isnothing.(matchidx), string("pval_", mt_col, "_threshold")] .= df_tops_detail[matchidx[.!isnothing.(matchidx)], string("pval_", mt_col, "_threshold")]
                        cols_float = [df_tops_detail[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_detail))]
                        df_tops_detail[:, cols_float] = round.(df_tops_detail[:, cols_float], sigdigits=6)
                        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.detail.txt.gz")), df_tops_detail, delim="\t", compress=true)
                    end
                end
            end
            println_to_file(string("    Done."), log_file)
        end
        cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
        df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.txt.gz")), df_tops_mt, delim="\t", compress=true)
        n_significant_phenotypes = sum(qval_mt .< FloatT(_args_fdr))
        logtxt = string("  *** ", n_significant_phenotypes, " out of ", length(qval_mt), " tested ", with_group ? "phenotype groups" : "phenotypes", " with qval_", mt_col, "<", _args_fdr, ".")
        println_to_file(logtxt, log_file)
    end
    @runif multiple_testing_method == "clipper" if length(_args_chrom) == 0
        mt_keep_index = .!isnan.(df_tops.pval_g1)
        df_tops_keep = df_tops[mt_keep_index, :]
        df_perm_keep = df_perm[mt_keep_index, :]
        df_perm_mt = df_perm_keep
        df_tops_mt = df_tops_keep
        if n_tests == 1
            df_tops_mt.qval_g1 .= NaN
            if n_perms < 1000
                res_clipper = Clipper(Matrix(df_perm_mt[:, 3:3]), Matrix(df_perm_mt[:, 4:end]), analysis="enrichment", procedure="GZ", contrast_score="max", FDR=[FloatT(_args_fdr)])
                qval_mt = res_clipper["q"]
                if !with_group
                    select!(df_tops_mt, Not("pval_g1_pheno"))
                else
                    select!(df_tops_mt, Not("pval_g1_group"))
                end
            else
                numsOfPermsMoreSig = sum(Matrix(df_perm_mt[:, 3:3]) .- Matrix(df_perm_mt[:, 4:end]) .<= 0, dims=2)
                pval_mt = vec((numsOfPermsMoreSig .+ 1) ./ (n_perms + 1))
                if !with_group 
                    df_tops_mt.pval_g1_pheno .= pval_mt
                else
                    df_tops_mt.pval_g1_group .= pval_mt
                end
                qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
            end
            df_tops_mt.qval_g1 .= qval_mt
            n_significant_phenotypes = sum(df_tops_mt.qval_g1 .< _args_fdr)
            logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(mt_keep_index), " tested ", with_group ? "phenotype groups" : "phenotypes", " with qval_g", n_tests, "<", _args_fdr, ".")
            println_to_file(logtxt, log_file)
        elseif n_tests == 2
            for i in 2:2
                idx_exp = ((i-1)+3):((i-1)+3)
                idx_back = ((i-1)+6):3:size(df_perm_mt, 2) |> Vector
                if n_perms < 1000
                    res_clipper = Clipper(Matrix(df_perm_mt[:, idx_exp]), Matrix(df_perm_mt[:, idx_back]), analysis="enrichment", procedure="GZ", contrast_score="max", FDR=[FloatT(_args_fdr)])
                    qval_mt = res_clipper["q"]
                    if !with_group
                        select!(df_tops_mt, Not("pval_g2_pheno"))
                    else
                        select!(df_tops_mt, Not("pval_g2_group"))
                    end
                else
                    numsOfPermsMoreSig = sum(Matrix(df_perm_mt[:, idx_exp]) .- Matrix(df_perm_mt[:, idx_back]) .<= 0, dims=2)
                    pval_mt = vec((numsOfPermsMoreSig .+ 1) ./ (n_perms + 1))
                    if !with_group 
                        df_tops_mt.pval_g2_pheno .= pval_mt
                    else
                        df_tops_mt.pval_g2_group .= pval_mt
                    end
                    qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
                end
                mt_col_name = ["qval_g1", "qval_g2", "qval_joint"][i]
                mt_threshold_col_name = ["pval_g1_threshold", "pval_g2_threshold", "pval_joint_threshold"][i]
                df_tops_mt[:, mt_col_name] .= NaN
                df_tops_mt[:, mt_col_name] .= qval_mt
                n_significant_phenotypes = sum(df_tops_mt[:, mt_col_name] .< FloatT(_args_fdr))
                logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(mt_keep_index), " tested ", with_group ? "phenotype groups" : "phenotypes", " with qval_", ["g1", "g2", "g_joint"][i], "<", _args_fdr, ".")
                println_to_file(logtxt, log_file)
            end
        end
        cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
        df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_qtl.txt.gz")), df_tops_mt, delim="\t", compress=true)
    end
end
