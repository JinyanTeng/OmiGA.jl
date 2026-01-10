function runOmiGA_gwas(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM::Union{Nothing,Dominance}=nothing, is_exclude_cis_chrom::Bool=false, exclude_window::Union{Nothing,Int}=nothing)
    is_write_text = _args_run_mode == "gwas"
    task_id = 0
    task_num = 0
    if length(_args_multi_task) == 2
        task_id = minimum(_args_multi_task) 
        task_num = maximum(_args_multi_task) 
        println_to_file(string("*** The analysis was split into [", task_num, "] tasks."), log_file)
        println_to_file(string("*** The current task ID is [", task_id, "]."), log_file)
        task_collects = collect(Iterators.partition(1:_struct_PHENO.n_phenotypes, ceil(Int, _struct_PHENO.n_phenotypes / task_num)))
        pheno_annotation = _struct_PHENO.annotation[task_collects[task_id], :]
        phenotype = _struct_PHENO.phenotype[:, task_collects[task_id]]
    elseif !isnothing(_args_lambda_phenos) & _args_lambda_only
        keep_pheno_index = floor.(Int, Vector(range(1, _struct_PHENO.n_phenotypes, length=_args_lambda_phenos)))
        pheno_annotation = _struct_PHENO.annotation[keep_pheno_index,:]
        phenotype = _struct_PHENO.phenotype[:,keep_pheno_index]
    else
        pheno_annotation = _struct_PHENO.annotation
        phenotype = _struct_PHENO.phenotype
    end
    n_phenotypes = size(phenotype, 2)
    if !isnothing(_struct_PRIOR_HerB)
        prior_h2 = _struct_PRIOR_HerB.h2
        prior_B = _struct_PRIOR_HerB.B
    else
        prior_h2 = prior_B = nothing
    end
    snp_annotation = _struct_GENO.annotation
    X_MME = _struct_COVAR.X_MME
    if isnothing(_args_h2_algo)
        h2_algo = "idul"
    else
        h2_algo = _args_h2_algo
    end
    _n, _X_c = size(X_MME)
    is_calcu_lambda = _args_run_mode == "gwas"
    if isnothing(_args_pheno_annot_file) & (_args_run_mode == "gwas")
        reserved_pthreshold = FloatT(1.0)
    else
        reserved_pthreshold = FloatT(_args_pval_threshold)
        if reserved_pthreshold < 1
            println_to_file(string(" * Only store results with P-value < ", reserved_pthreshold, " to disk."), log_file)
        end
    end
    tstat_threshold = quantile(WaldTest(1, _n - _X_c - 1), 1 - reserved_pthreshold)
    if n_phenotypes > 1000
        @info string("A total of ", n_phenotypes, " phenotypes used for analysis, this may causa memory overflow!\nRecommend to split the analysis to a sets of small tasks using the '--multi-task' option.")
    end
    qtl_map_model = _args_qtl_map_model
    qtl_map_algo = _args_qtl_map_algo
    n_tests = 1
    n_grms = 1
    if qtl_map_model == "a"
        glo_EA = eigen(diagm(ones(FloatT, _n)))
        genotype = _struct_GENO.genotype
    end
    if qtl_map_model == "a+A"
        glo_EA = _struct_KIN.EA
        genotype = _struct_GENO.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
            glo_KS = _struct_KIN.GRM
        end
    end
    if qtl_map_model == "d+D"
        glo_EA = _struct_KIN.domEA
        genotype = _struct_DOM.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
            glo_KS = _struct_KIN.domGRM
        end
    end
    if qtl_map_model == "d+A"
        glo_EA = _struct_KIN.EA
        genotype = _struct_DOM.genotype
        if isnothing(prior_h2) & (h2_algo == "minmax")
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
    @runif n_grms > 1 error("[Exit] Current version supported only one GRM used in GWAS mode.")
    n_snps_est_lambda = _args_lambda_snps
    if n_snps_est_lambda > size(snp_annotation, 1)
        n_snps_est_lambda = size(snp_annotation, 1)
    end
    lambda_snp_index = floor.(Int, Vector(range(1, _struct_GENO.n_snps, length=n_snps_est_lambda)))
    if _args_lambda_only
        genotype = copy(genotype[:, lambda_snp_index])
        snp_annotation = copy(snp_annotation[lambda_snp_index, :])
        snp_annotation.index .= 1:size(snp_annotation, 1)
        lambda_snp_index .= 1:length(lambda_snp_index)
    end
    if is_calcu_lambda
        println_to_file(string(" * GC lambda will be estimated using P-values from ", n_snps_est_lambda, " variants uniformly extracted from whole chromosomes being tested."), log_file)
        lambda_p_mat = fill(NAN, n_snps_est_lambda, n_phenotypes)
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
    if (isnothing(prior_h2) & (h2_algo == "idul")) | (qtl_map_algo == "idul")
        if _args_preadj_covar
            xQ = glo_EA.vectors' * X_MME[:, 1:1]
        else
            xQ = glo_EA.vectors' * X_MME 
        end
    end
    Vi = zeros(FloatT, _n, _n)
    Di = Diagonal(zeros(FloatT, _n, _n))
    EAvec_Di = zeros(FloatT, _n, _n)
    if qtl_map_model in ["a+A+D", "d+A+D", "a+d+A+D"]
        EAvec = zeros(FloatT, _n, _n)
    end
    n_eff_varaint = 1000000
    df_tops = DataFrame()
    pheno_annotation.index = 1:nrow(pheno_annotation)
    chroms = string.(unique(snp_annotation.chromosome))
    if length(chroms) == 0
        error("Non-overlapping chromosomes between genotype and phenotype data!")
    end
    _task_write_full_pairs = @task @info "Using Tasks"
    for chrom in chroms 
        snp_index_chrom = snp_annotation.chromosome .== chrom
        _snp_annot = snp_annotation[snp_index_chrom, :]
        _snp_annot.reindex = 1:size(_snp_annot, 1) 
        if is_calcu_lambda
            chrom_lambda_snp_reindex = _snp_annot.reindex[_snp_annot.index.∈(lambda_snp_index,)]
            chrom_lambda_snp_matindex = findall(lambda_snp_index .∈ (_snp_annot.index,))
        end
        n_chrom_snps = sum(snp_index_chrom)
        _gene_annot = pheno_annotation 
        _n_phenos = length(_gene_annot.pheno_id) 
        @timeit to "Pull gwas-SNP genotypes" begin
            _2p = 2 .* _snp_annot.af
            _2pq = _2p .* (1 .- _snp_annot.af)
            if is_sample_byrow
                _2p = _2p' |> Matrix
                _2pq = _2pq' |> Matrix
            end
        end
        if qtl_map_algo == "idul"
            if _args_preadj_covar
                rotate_yQ = glo_EA.vectors' * get_lm_residuals(X_MME, phenotype)
            else
                rotate_yQ = glo_EA.vectors' * phenotype 
            end
            _2p = 2 .* _snp_annot.af
            _2pq = _2p .* (1 .- _snp_annot.af)
            if is_sample_byrow
                _2p = _2p' |> Matrix
                _2pq = _2pq' |> Matrix
                chrom_genotype = similar(genotype[:, snp_index_chrom], FloatT)
                chrom_genotype .= genotype[:, snp_index_chrom]
                chrom_genotype .-= _2p
                chrom_SQ = glo_EA.vectors' * chrom_genotype
                chrom_genotype .= chrom_SQ
                GC.gc()
            end
        end
        @timeit to "temporary DataFrame to store associations" begin
            _chrom_df_full = spzeros(FloatT, n_chrom_snps, 3 * _n_phenos)
            _n_genes_batch = ceil(Int, _n_phenos / 10)
            iter_collects = collect(Iterators.partition(1:_n_phenos, _n_genes_batch))
            _df_test = zeros(FloatT, n_chrom_snps, length(iter_collects[1]) * 3) 
        end
        begin
            SQ_h1inv_init_t_y1_init = zeros(FloatT, n_chrom_snps)
            X_init_SQ_h1inv_init = zeros(FloatT, size(xQ, 2) + 1, n_chrom_snps)
        end
        GC.gc()
        time_start = now()
        println_to_file(string("Chromosome: ", chrom, ", start at: ", time_start), log_file)
        for batchi in eachindex(iter_collects)
            fill!(_df_test, 0)
            for i in iter_collects[batchi]
                gene = _gene_annot.pheno_id[i]
                println_to_file(string("    PHENO: ", i, "/", _n_phenos, " <", gene,">"), log_file)
                @runif is_exclude_cis_chrom (_gene_annot.chrom[i] == chrom) && continue
                exppheno = phenotype[:, findfirst(pheno_annotation.pheno_id .== gene)] 
                if qtl_map_algo == "standard"
                    @timeit to "Global heritability" if !isnothing(prior_h2)
                        index_prior_h2 = findfirst(prior_h2.pheno_id .== gene)
                        Σ_i = [prior_h2.glo_vg1[index_prior_h2], prior_h2.glo_vg2[index_prior_h2], prior_h2.glo_ve[index_prior_h2]]
                        glo_vc_B = prior_B[index_prior_h2, 3:end] |> Vector
                    else
                        if h2_algo == "minmax"
                            glo_vc = getMRVCModel(glo_KS, exppheno; X=X_MME)
                        elseif h2_algo == "idul"
                            glo_vc = get_IDUL_VarianceComponent(glo_EA, exppheno, X_MME, xQ=xQ)
                        end
                        Σ_i = [glo_vc[:ΣG][1], glo_vc[:Σe]]
                        glo_vc_B = glo_vc[:Fix_eff]
                    end
                    @timeit to "Vi" if qtl_map_model in ["a+A", "d+D", "d+A"]
                        getVinv!(Vi, Di, glo_EA, Σ_i[[1, 3]], EAvec_Di)
                    elseif qtl_map_model in ["a+A+D", "d+A+D", "a+d+A+D"]
                        getVinv!(Vi, Σ_i, _struct_PartialV, Di, EAvec_Di, EAvec)
                    end
                end
                current_df_cols = (((i-1)*3+1):((i-1)*3+3)) .- ((batchi - 1) * size(_df_test, 2))
                if qtl_map_algo == "idul"
                    yQ = rotate_yQ[:, i]
                    if _args_linear_model
                        init_eta = zero(FloatT)
                    else
                        idul_prior_index = sample(1:n_chrom_snps, 10) 
                        df_subs = zeros(FloatT, 10, 3 * n_tests + 2)
                        @timeit to "Prior idul" idul_assoc_test_plus!(df_subs, chrom_SQ[:, idul_prior_index], yQ, xQ, glo_EA.values; max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), return_detail=true, is_calcu_pv=false)
                        init_eta = Statistics.mean(df_subs[.!isnan.(df_subs[:, 2]), 2])
                        @runif _args_verbose println(string("    η: ",init_eta))
                    end
                    if !_args_exact_map
                        @timeit to "Idul test" idul_global_test_approx!(view(_df_test, :, current_df_cols), chrom_SQ, yQ, xQ, glo_EA.values, SQ_h1inv_init_t_y1_init, X_init_SQ_h1inv_init; init_eta=init_eta, is_calcu_pv=false) 
                        chrom_SQ .= chrom_genotype
                    else
                        @timeit to "Idul test" idul_assoc_test_plus!(view(_df_test, :, current_df_cols), chrom_SQ, yQ, xQ, glo_EA.values; init_eta=init_eta, max_iter=_args_idul_max_iter, thre=FloatT.(_args_idul_converge), is_calcu_pv=false)
                        chrom_SQ .= chrom_genotype
                    end
                    if is_calcu_lambda
                        lambda_p_mat[chrom_lambda_snp_matindex, i] .= view(_df_test, :, current_df_cols)[chrom_lambda_snp_reindex, 3]
                    end
                    if !_args_lambda_only 
                        if tstat_threshold > 0
                            df_test_omit_rows = view(_df_test, :, current_df_cols)[:, 3] .< tstat_threshold
                            fill!(view(_df_test, df_test_omit_rows, current_df_cols), 0)
                        elseif tstat_threshold == 0 && !is_write_text 
                            df_test_omit_rows = iszero.(view(_df_test, :, current_df_cols)[:, 3])
                            fill!(view(_df_test, df_test_omit_rows, current_df_cols), 0)
                        end
                        if !isnothing(exclude_window) && chrom == _gene_annot.chrom[i]
                            exclude_snp_indices = get_cis_snp_info(_gene_annot, _snp_annot, gene, exclude_window; window_type=_args_window_type)[1].reindex
                            fill!(view(_df_test, exclude_snp_indices, current_df_cols), 0)
                        end
                    end
                elseif qtl_map_algo == "standard"
                    @timeit to "Adjusted phenotype" yadj = exppheno 
                    @timeit to "Get associations" if qtl_map_model in ["a+d+A+D"]
                        calcu_two_assoc_test!(df_test, Vi, genotype, domGenotype, yadj)
                    else
                        calcu_one_assoc_test!(df_test, Vi, genotype, _2pq, yadj, Vi_yc, batch_collects, S0, Vi_S0, σ0, S0_Vi_yc, S1, Vi_S1, σ1, S1_Vi_yc, is_sample_byrow)
                    end
                end
                @runif _args_debug if i % 10 == 0
                    println_to_file(string("*** Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
                    println_to_file(string(to), log_file)
                    println()
                end
            end
            if !_args_lambda_only
                start_coli = (batchi - 1) * size(_df_test, 2) + 1
                end_coli = start_coli + length(iter_collects[batchi]) * 3 - 1
                @timeit to "Append results to SparseArrays" _chrom_df_full[:, start_coli:end_coli] .= _df_test[:, 1:(length(iter_collects[batchi])*3)]
            end
        end
        @runif !_args_lambda_only @timeit to "Write the full summary to disk" begin
            if chrom != chroms[1]
                wait(_task_write_full_pairs) 
            end
            task_df_full = copy(_chrom_df_full)
            task_snp_annot = copy(_snp_annot)
            task_chrom = chrom
            _task_write_full_pairs = @spawn begin
                write_summary_file(task_df_full, task_snp_annot, pheno_annotation, task_chrom, is_write_text, _n, _X_c)
            end
            if chrom == chroms[end]
                wait(_task_write_full_pairs)
            end
        end
        println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
        if _args_debug
            println_to_file(string(to), log_file)
            println()
        end
    end
    if is_calcu_lambda
        lambda_est = vec(get_GC_lambda.(ccdf(WaldTest(1, _n - size(X_MME,2) - 1), median(lambda_p_mat, dims=1))))
        lambda_DF = DataFrame(:pheno_id => pheno_annotation.pheno_id, :lambda => lambda_est)
        if length(_args_multi_task) == 2
            out_text = joinpath(_args_output_dir, string(_args_out_prefix, ".GC_lambda.task_", task_id, ".txt"))
        else
            out_text = joinpath(_args_output_dir, string(_args_out_prefix, ".GC_lambda.txt"))
        end
        CSV.write(out_text, lambda_DF, delim="\t", compress=false)
    end
end
function write_summary_file(task_df_full, _snp_annot, pheno_annotation, chrom, is_write_text, _n, _X_c)
    if !is_write_text
        if length(_args_multi_task) == 2
            out_jld = joinpath(_args_output_dir, string(_args_out_prefix, ".assoc_pairs.", chrom, ".task_", minimum(_args_multi_task), ".jld2"))
        else
            out_jld = joinpath(_args_output_dir, string(_args_out_prefix, ".assoc_pairs.", chrom, ".jld2"))
        end
        save(out_jld, "data", task_df_full)
    else
        task_df_full_mat = task_df_full |> Matrix
        zeros_index = findall(iszero, task_df_full)
        task_df_full_mat[zeros_index] .= NAN
        colnames = ["chrom", "position", "variant_id", "alt", "ref", "af", "het_rate", string.(repeat(pheno_annotation.pheno_id, inner=3), repeat(["(beta)", "(se)", "(pval)"], Int(size(task_df_full_mat, 2) / 3)))...]
        task_df_full_df = hcat(_snp_annot[:, [1, 4, 2, 5, 6, 7, 8]], DataFrame(task_df_full_mat, :auto))
        rename!(task_df_full_df, colnames)
        pval_cols = findall(.!isnothing.(match.(r"pval", colnames)))
        for pcol in pval_cols
            task_df_full_df[!,pcol] = ccdf(WaldTest(1, _n - _X_c - 1), task_df_full_df[:,pcol])
        end
        cols_float = [task_df_full_df[1, x] isa AbstractFloat for x in range(1, ncol(task_df_full_df))]
        task_df_full_df[:, cols_float] .= round.(task_df_full_df[:, cols_float], sigdigits=6)
        task_df_full_df.af .= round.(task_df_full_df.af, digits=3)
        task_df_full_df.het_rate .= round.(task_df_full_df.het_rate, digits=3)
        if length(_args_multi_task) == 2
            out_text = joinpath(_args_output_dir, string(_args_out_prefix, ".assoc_pairs.", chrom, ".task_", task_id, ".txt.gz"))
        else
            out_text = joinpath(_args_output_dir, string(_args_out_prefix, ".assoc_pairs.", chrom, ".txt.gz"))
        end
        @runif _args_use_gzip out_text = replace(out_text, r".gz$"=>"")
        CSV.write(out_text, task_df_full_df, delim="\t", compress=!_args_use_gzip)
        @runif _args_use_gzip run(`gzip -f $out_text`; wait = chrom == chroms[end])
        task_df_full_df = nothing
        GC.gc()
    end
end
