function runOmiGA_her(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR; _struct_DOM::Union{Nothing,Dominance}=nothing)
    phenotype = _struct_PHENO.phenotype
    pheno_annotation = _struct_PHENO.annotation
    X_MME = copy(_struct_COVAR.X_MME)
    if !isnothing(_args_geno_file_prefix)
        genotype = _struct_GENO.genotype
        snp_annotation = _struct_GENO.annotation
        if !isnothing(_struct_DOM)
            domGenotype = _struct_DOM.genotype
        end
        @runif USE_Float32 snp_annotation.af .= Float64.(snp_annotation.af)
    else
        genotype = nothing
        domGenotype = nothing
        snp_annotation = nothing
    end
    convergence_parameters = FloatT(_args_dpars)
    convergence_dlogL = FloatT(NAN)
    if isnothing(genotype)
        n_glo_snps = 0
    else
        n_glo_snps = _struct_GENO.n_snps
    end
    _n, _X_c = size(X_MME)
    _args_glo_h2 = false
    _args_trans_h2 = false
    _args_cis_h2 = false
    h2_model = isnothing(_args_h2_model) ? "undef" : _args_h2_model
    glo_h2_model = h2_model
    trans_h2_model = h2_model
    cis_h2_model = h2_model
    if h2_model in ["Ag", "Dg", "Ag+Dg"]
        _args_glo_h2 = true
    elseif h2_model in ["At", "Dt", "Ag+Dt"]
        _args_trans_h2 = true
    elseif h2_model in ["Ac", "Dc", "Ag+Ac", "Ag+Dc", "At+Ac"]
        _args_cis_h2 = true
    else
        @info "Use user-defined heritability model."
    end
    if isnothing(_args_h2_algo)
        h2_algo = "em_aireml"
    else
        h2_algo = _args_h2_algo
    end
    if h2_model == "undef"
        h2_algo = "minmax" 
    end
    glo_h2_algo = h2_algo
    @runif _args_glo_h2 if glo_h2_model in ["Ag+Dg"]
        if !(h2_algo in ["em_aireml", "minmax"])
            h2_algo = glo_h2_algo = "em_aireml"
            @info string("--h2-algo was set to 'em_aireml' for ", glo_h2_model)
        end
    end
    trans_h2_algo = h2_algo
    @runif _args_trans_h2 if trans_h2_model in ["Ag+Dt"]
        if !(h2_algo in ["em_aireml", "minmax"])
            h2_algo = trans_h2_algo = "em_aireml"
            @info string("--h2-algo was set to 'em_aireml' for ", trans_h2_model)
        end
    end
    cis_h2_algo = h2_algo
    @runif _args_cis_h2 if cis_h2_model in ["Ag+Ac", "Ag+Dc", "At+Ac"]
        if !(h2_algo in ["em_aireml", "minmax"])
            h2_algo = cis_h2_algo = "em_aireml"
            @info string("--h2-algo was set to 'em_aireml' for ", cis_h2_model)
        end
    end
    In = diagm(ones(FloatT, _n)) 
    pre_calcu_xQ_glo = glo_h2_algo in ["idul", "idul_aireml"] 
    pre_calcu_xQ_trans = trans_h2_algo in ["idul", "idul_aireml"] 
    pre_calcu_eigen = true
    if pre_calcu_xQ_glo
        if glo_h2_model == "Ag"
            glo_xQ = _struct_KIN.EA.vectors' * X_MME 
        elseif glo_h2_model == "Dg"
            glo_xQ = _struct_KIN.domEA.vectors' * X_MME 
        end
    end
    if _args_glo_h2
        if glo_h2_model == "Ag"
            glo_KSs = [_struct_KIN.GRM, In]
            glo_EA = ifelse(pre_calcu_eigen, _struct_KIN.EA, nothing)
        elseif glo_h2_model == "Dg"
            glo_KSs = [_struct_KIN.domGRM, In]
            glo_EA = ifelse(pre_calcu_eigen, _struct_KIN.domEA, nothing)
        elseif glo_h2_model == "Ag+Dg"
            glo_KSs = [_struct_KIN.GRM, _struct_KIN.domGRM, In]
            glo_EA = nothing
        end
        glo_n_cmp = length(glo_KSs)
    end
    _args_PreV = _args_glo_h2 & (glo_h2_model == "Ag+Dg") & (glo_h2_algo == "em_aireml") | _args_export_paed_file
    @timeit to "Pre-calcu V" if _args_PreV
        path_partialv = _args_path_partialv
        if !isnothing(path_partialv)
            _struct_PartialV = load(path_partialv)["data"]
        else
            println_to_file("Performing PAED ...", log_file)
            @runif !_args_export_paed_file println_to_file("Note: Please specify '--export-paed-file' if you want to map QTL using model with two GRMs (i.e., d+A+D) later.", log_file)
            G1 = _struct_KIN.GRM |> Matrix{Float32}
            G2 = _struct_KIN.domGRM |> Matrix{Float32}
            ratio_step = Float32(_args_paed_step)
            ratio_max = Float32(_args_paed_max)
            vec_ratio = range(1, ratio_max, step=ratio_step) |> Vector{Float32}
            eig_rows = length(vec_ratio) * 2 * _n
            cGRM_eigvals = zeros(Float32, eig_rows)
            cGRM_eigvecs = zeros(Float32, eig_rows, _n)
            cGRM = similar(G1, Float32)
            for i in eachindex(vec_ratio)
                v1_v2_ratio = vec_ratio[i]
                v1 = v1_v2_ratio / (1 + v1_v2_ratio)
                v2 = 1 - v1
                cGRM .= v1 * G1 + v2 * G2
                cval, cvec = eigen!(cGRM)
                j = i * 2 - 1
                range_inds = _n*(j-1)+1:_n*j
                vals_inds = CartesianIndices((range_inds,))
                vecs_inds = CartesianIndices((range_inds, 1:_n))
                copyto!(cGRM_eigvals, vals_inds, cval, CartesianIndices(cval))
                copyto!(cGRM_eigvecs, vecs_inds, cvec, CartesianIndices(cvec))
                cGRM .= v2 * G1 + v1 * G2
                cval, cvec = eigen!(cGRM)
                j = i * 2
                range_inds = _n*(j-1)+1:_n*j
                vals_inds = CartesianIndices((range_inds,))
                vecs_inds = CartesianIndices((range_inds, 1:_n))
                copyto!(cGRM_eigvals, vals_inds, cval, CartesianIndices(cval))
                copyto!(cGRM_eigvecs, vecs_inds, cvec, CartesianIndices(cvec))
            end
            _struct_PartialV = PartialV{Float32}(vec_ratio, cGRM_eigvals, cGRM_eigvecs)
            if _args_export_paed_file
                save(joinpath(_args_output_dir, string(_args_out_prefix, ".paed.jld2")), "data", _struct_PartialV)
            end
            println_to_file("PAED done.", log_file)
            cGRM_eigvals = nothing
            cGRM_eigvecs = nothing
            G1 = nothing
            G2 = nothing
            cGRM = nothing
            GC.gc()
        end
    end
    if _args_trans_h2
        if trans_h2_model == "At"
            trans_KSs = [zeros(FloatT, _n, _n), In]
        elseif trans_h2_model == "Dt"
            trans_KSs = [zeros(FloatT, _n, _n), In]
        elseif trans_h2_model == "Ag+Dt"
            trans_KSs = [_struct_KIN.GRM, zeros(FloatT, _n, _n), In]
        end
        trans_n_cmp = length(trans_KSs)
    end
    if _args_cis_h2
        if cis_h2_model == "Ac"
            cis_KSs = [zeros(FloatT, _n, _n), In]
        elseif cis_h2_model == "Dc"
            cis_KSs = [zeros(FloatT, _n, _n), In]
        elseif cis_h2_model == "Ag+Ac"
            cis_KSs = [_struct_KIN.GRM, zeros(FloatT, _n, _n), In]
        elseif cis_h2_model == "Ag+Dc"
            cis_KSs = [_struct_KIN.GRM, zeros(FloatT, _n, _n), In]
        elseif cis_h2_model == "At+Ac"
            cis_KSs = [zeros(FloatT, _n, _n), zeros(FloatT, _n, _n), In]
        end
        cis_n_cmp = length(cis_KSs)
        if _args_use_low_rank_approx
            if _args_approx_rank isa AbstractFloat
                approx_rank = Int(round(_n * _args_approx_rank))
            elseif _args_approx_rank isa Integer
                approx_rank = _args_approx_rank
            end
            if isnothing(_args_max_approx_rank)
                max_approx_rank = 1.0
            else
                max_approx_rank = _args_max_approx_rank
            end
        end
    end
    if h2_model == "undef"
        multi_KSs = [_struct_KIN.GRM..., In]
        multi_n_cmp = length(multi_KSs)
    end
    if _args_glo_h2
        glo_Σ_i = fill(NAN, glo_n_cmp)
        glo_h2 = fill(NAN, glo_n_cmp - 1)
        glo_h2_se = fill(NAN, glo_n_cmp - 1)
        glo_fix_eff = zeros(FloatT, _X_c)
    end
    if _args_cis_h2
        cis_Σ_i = fill(NAN, cis_n_cmp)
        cis_h2 = fill(NAN, cis_n_cmp - 1)
        cis_h2_se = fill(NAN, cis_n_cmp - 1)
    end
    if _args_trans_h2
        trans_Σ_i = fill(NAN, trans_n_cmp)
        trans_h2 = fill(NAN, trans_n_cmp - 1)
        trans_h2_se = fill(NAN, trans_n_cmp - 1)
    end
    if h2_model == "undef"
        multi_Σ_i = fill(NAN, multi_n_cmp)
        multi_h2 = fill(NAN, multi_n_cmp - 1)
        multi_h2_se = fill(NAN, multi_n_cmp - 1)
    end
    n_trans_snps = 0
    n_cis_snps = 0
    _args_save_memory = h2_algo in ["em_aireml", "idul_aireml"]
    if _args_save_memory
        _Vi_X = zeros(FloatT, _n, _X_c)
        _Xt_Vi_X_i = zeros(FloatT, _X_c, _X_c)
        _Py = zeros(FloatT, _n)
        _P = zeros(FloatT, _n, _n)
        _PA = zeros(FloatT, _n, _n)
        _Vi = zeros(FloatT, _n, _n) 
        _cVi = zeros(FloatT, _n, _n) 
        _Di = Diagonal(zeros(FloatT, _n, _n))
        _EAvec_Di = zeros(FloatT, _n, _n)
        _EAvec = zeros(FloatT, _n, _n)
        if _args_glo_h2
            glo_Hi = zeros(FloatT, glo_n_cmp, glo_n_cmp)
        end
        if _args_trans_h2
            trans_Hi = zeros(FloatT, trans_n_cmp, trans_n_cmp)
        end
        if _args_cis_h2
            cis_Hi = zeros(FloatT, cis_n_cmp, cis_n_cmp)
        end
    end
    if h2_model == "undef"
        n_grm = _struct_KIN.n_grm
    else
        n_grm = length(split(h2_model, "+"))
    end
    chroms = string.(unique(pheno_annotation.chrom))
    for chrom in chroms
        @runif _args_cis_h2 | _args_trans_h2 snp_index_chrom = snp_annotation.chromosome .== chrom
        _gene_annot = pheno_annotation[pheno_annotation.chrom.==chrom, :]
        _n_genes = length(_gene_annot.pheno_id) 
        _df_hers = DataFrame(pheno_id=_gene_annot.pheno_id,
                chrom=chrom,
            )
        @runif _args_glo_h2 _df_hers.glo_snps .= n_glo_snps
        @runif _args_trans_h2 | _args_cis_h2 _df_hers.trans_snps .= 0
        @runif _args_cis_h2 _df_hers.cis_snps .= 0
        _df_hers.algo .= h2_algo
        _df_hers.model .= h2_model
        [_df_hers[:,string("vg",i)] .= NAN for i in 1:n_grm]
        _df_hers.ve .= NAN
        [
            begin
                _df_hers[:, string("h2_g", i)] .= NAN
                _df_hers[:, string("h2se_g", i)] .= NAN
                _df_hers[:, string("h2p_g", i)] .= NaN 
            end for i in 1:n_grm
        ]
        if _args_use_low_rank_approx & (cis_h2_algo == "em_aireml") & (cis_h2_model in ["Ag+Ac", "Ag+Dc", "At+Ac"])
            _df_hers.cis_rank .= 0
        end
        if (h2_algo in ["em_aireml", "idul_aireml", "minmax"]) | (h2_model == "undef")
            converged_flag = false
            _df_hers.converged .= converged_flag
        end
        _df_covar_B = hcat(DataFrame(pheno_id=_gene_annot.pheno_id, chrom=chrom), DataFrame(Dict(Symbol(lpad(i, 5, '0')) => repeat([NAN], _n_genes) for i in 0:(_X_c-1))))
        rename!(_df_covar_B, ["pheno_id", "chrom", string.("B", 0:(_X_c-1))...])
        if _args_export_rand_eff
            _df_u_hat = hcat(DataFrame(pheno_id=repeat(_gene_annot.pheno_id, inner=n_grm + 1), chrom=repeat([chrom], inner=_n_genes * (n_grm + 1)), eff=repeat([string.("u", 1:n_grm)..., "e"], outer=_n_genes)), DataFrame(fill(NAN, _n_genes * (n_grm + 1), _n), _struct_PHENO.iid))
        end
        @timeit to "Pull chrom-genotypes" if _args_cis_h2
            if !is_pre_adj_genotype
                _2p = 2 * snp_annotation.af[snp_index_chrom] |> Vector{FloatT}
                if cis_h2_model in ["Ac", "Ag+Ac", "At+Ac"]
                    chrom_genotype = similar(genotype[:, snp_index_chrom], FloatT)
                    chrom_genotype .= genotype[:, snp_index_chrom] .- _2p' 
                end
                if cis_h2_model in ["Dc", "Ag+Dc"]
                    _2pq = _2p .* (1 .- snp_annotation.af[snp_index_chrom]) |> Vector{FloatT}
                    chrom_domGenotype = similar(domGenotype[:, snp_index_chrom], FloatT)
                    chrom_domGenotype .= domGenotype[:, snp_index_chrom] .- _2pq' 
                end
            else
                error("No support pre-adjusted genotype data!")
            end
            _chrom_snp_annot = snp_annotation[snp_index_chrom, :]
            _chrom_snp_annot.reindex = 1:size(_chrom_snp_annot, 1)
        end
        if _args_trans_h2 | (cis_h2_model == "At+Ac")
            @timeit to "(trans) Pull tested Chromosome" begin 
                _trans_annot = snp_annotation[.!snp_index_chrom, :]
                n_trans_snps = nrow(_trans_annot) 
                trans_afs = snp_annotation.af[_trans_annot.index]
            end
            if cis_h2_model in ["At+Ac"] 
                @timeit to "Build At" cis_KSs[1] .= getG_fast(genotype[:, _trans_annot.index]; alpha=_args_grm_alpha, afs=trans_afs, return_eltype=FloatT)
            end
            if trans_h2_model in ["At"] 
                @timeit to "Build At" trans_KSs[1] .= getG_fast(genotype[:, _trans_annot.index]; alpha=_args_grm_alpha, afs=trans_afs, return_eltype=FloatT)
                trans_EA = ifelse(pre_calcu_eigen, eigen(trans_KSs[1]), nothing)
                if pre_calcu_xQ_trans
                    trans_xQ = mat_mul(trans_EA.vectors', X_MME) 
                end
            end
            if trans_h2_model in ["Dt", "Ag+Dt"] 
                @timeit to "Build Dt" trans_KSs[trans_n_cmp-1] .= getG_dominance(domGenotype[:, _trans_annot.index], code_type="dominance", afs=trans_afs, byrow=true)
                if trans_h2_model in ["Dt"]
                    trans_EA = ifelse(pre_calcu_eigen, eigen(trans_KSs[trans_n_cmp-1]), nothing)
                    if pre_calcu_xQ_trans
                        trans_xQ = mat_mul(trans_EA.vectors', X_MME) 
                    end
                end
            end
        end
        GC.gc()
        time_start = now()
        println_to_file(string("Chromosome: ", chrom, ", start at: ", time_start), log_file)
        for i in 1:_n_genes
            gene = _gene_annot.pheno_id[i]
            println_to_file(string("    PHENO: ", i, "/", _n_genes, " <", gene,">"), log_file)
            exppheno = phenotype[:, findfirst(pheno_annotation.pheno_id .== gene)] 
            if _args_glo_h2
                @timeit to "Global heritability" if glo_h2_model in ["Ag", "Dg"]
                    if glo_h2_algo == "minmax"
                        glo_vc = getMRVCModel(glo_KSs, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                    elseif glo_h2_algo == "em_aireml"
                        @notimeit to glo_vc = reml(glo_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, glo_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=glo_EA, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    elseif glo_h2_algo == "idul_aireml"
                        @timeit to "idul prior" prior_glo_vc = get_IDUL_VarianceComponent(glo_EA, exppheno, X_MME, xQ=glo_xQ, max_iter=5, thre=FloatT(1e-4))
                        @notimeit to glo_vc = reml(glo_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, glo_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=glo_EA, varcmp=[prior_glo_vc[:ΣG][1], prior_glo_vc[:Σe]], em_update=true, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    elseif glo_h2_algo == "idul"
                        glo_vc = get_IDUL_VarianceComponent(glo_EA, exppheno, X_MME, xQ=glo_xQ)
                    end
                    glo_Σ_i .= [glo_vc[:ΣG][1], glo_vc[:Σe]]
                elseif glo_h2_model in ["Ag+Dg"]
                    if glo_h2_algo == "em_aireml"
                        @notimeit to glo_vc = reml(glo_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, glo_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, partialV=_struct_PartialV, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    end
                    glo_Σ_i .= [glo_vc[:ΣG][1], glo_vc[:ΣG][2], glo_vc[:Σe]]
                end
                @runif glo_h2_algo == "minmax" println(glo_Σ_i)
                glo_h2 .= glo_vc[:h2]
                glo_h2_se .= glo_vc[:h2_se]
                glo_fix_eff .= glo_vc[:Fix_eff]
                @timeit to "Append df_covar_B" _df_covar_B[i, 3:end] = glo_fix_eff'
                _df_hers.glo_snps[i] = n_glo_snps
                _df_hers.vg1[i] = glo_Σ_i[1]
                _df_hers.ve[i] = glo_Σ_i[glo_n_cmp]
                _df_hers.h2_g1[i] = glo_h2[1]
                _df_hers.h2se_g1[i] = glo_h2_se[1]
                _df_hers.h2p_g1[i] = 0.5 * ccdf(Chisq(1), (glo_h2[1] / glo_h2_se[1])^2)
                if glo_n_cmp == 3
                    _df_hers.vg2[i] = glo_Σ_i[2]
                    _df_hers.h2_g2[i] = glo_h2[2]
                    _df_hers.h2se_g2[i] = glo_h2_se[2]
                    _df_hers.h2p_g2[i] = 0.5 * ccdf(Chisq(1), (glo_h2[2] / glo_h2_se[2])^2)
                end
                if h2_algo in ["em_aireml", "idul_aireml", "minmax"]
                    _df_hers.converged[i] = glo_vc[:converged]
                end
                if _args_export_rand_eff
                    _df_u_hat[findall(gene .== _df_u_hat.pheno_id), 4:end] .= glo_vc[:u]'
                end
            end
            if _args_cis_h2 
                @timeit to "Keep SNPs within cis-region" cissnps_annot, n_cis_snps = get_cis_snp_info(_gene_annot, _chrom_snp_annot, gene, _args_cis_window)
                fill!(cis_h2, NAN)
                fill!(cis_h2_se, NAN)
                fill!(cis_Σ_i, NAN)
                approx_ranki = 0
                if n_cis_snps > 0
                    if cis_h2_model in ["Ac", "Ag+Ac", "At+Ac"]
                        @timeit to "Build Ac" cis_KSs[cis_n_cmp-1] .= getG_fast(chrom_genotype[:, cissnps_annot.reindex]; alpha=_args_grm_alpha, afs=cissnps_annot.af, return_eltype=FloatT)
                    end
                    if cis_h2_model in ["Dc", "Ag+Dc"]
                        @timeit to "Build Dc" cis_KSs[cis_n_cmp-1] .= getG_dominance(chrom_domGenotype[:, cissnps_annot.reindex], code_type="dominance", adjusted=true, afs=cissnps_annot.af, byrow=true)
                    end
                    @timeit to "Cis heritability" if cis_h2_model in ["Ac", "Dc"]
                        cis_EA = eigen(cis_KSs[1])
                        if cis_h2_algo == "minmax"
                            cis_vc = getMRVCModel(cis_KSs, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                        elseif cis_h2_algo == "em_aireml"
                            @notimeit to cis_vc = reml(cis_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, cis_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=cis_EA, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                        elseif cis_h2_algo == "idul_aireml"
                            prior_cis_vc = get_IDUL_VarianceComponent(cis_EA, exppheno, X_MME, max_iter=5, thre=FloatT(1e-4))
                            @notimeit to cis_vc = reml(cis_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, cis_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=cis_EA, varcmp=[prior_cis_vc[:ΣG][1], prior_cis_vc[:Σe]], em_update=true, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                        elseif cis_h2_algo == "idul"
                            cis_vc = get_IDUL_VarianceComponent(cis_EA, exppheno, X_MME)
                        end
                        cis_Σ_i .= [cis_vc[:ΣG][1], cis_vc[:Σe]]
                    elseif cis_h2_model in ["Ag+Ac", "Ag+Dc", "At+Ac"]
                        if cis_h2_algo == "em_aireml"
                            if _args_use_low_rank_approx
                                println("** LowRankApprox **")
                                @timeit to "LowRankApproxUDV" LRA_factors = LowRankApproxUDV(cis_KSs[2], _args_approx_method, approx_rank, _args_avg_2norm_threshold, max_approx_rank)
                                @notimeit to cis_vc = reml([cis_KSs[1], LRA_factors.M, cis_KSs[3]], exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, cis_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=_struct_KIN.EA, LRA_factors=LRA_factors, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                                approx_ranki = LRA_factors.rank
                            else
                                @notimeit to cis_vc = reml(cis_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, cis_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; 
                                pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                            end
                        elseif cis_h2_algo == "minmax"
                            cis_vc = getMRVCModel(cis_KSs, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                        end
                        cis_Σ_i .= [cis_vc[:ΣG][1], cis_vc[:ΣG][2], cis_vc[:Σe]]
                        @runif cis_h2_algo == "minmax" println(cis_Σ_i)
                    end
                    cis_h2 .= cis_vc[:h2]
                    cis_h2_se .= cis_vc[:h2_se]
                    if h2_algo in ["em_aireml", "idul_aireml", "minmax"]
                        _df_hers.converged[i] = cis_vc[:converged]
                    end
                    _df_hers.vg1[i] = cis_Σ_i[1]
                    _df_hers.ve[i] = cis_Σ_i[cis_n_cmp]
                    _df_hers.h2_g1[i] = cis_h2[1]
                    _df_hers.h2se_g1[i] = cis_h2_se[1]
                    _df_hers.h2p_g1[i] = ccdf(Chisq(1), (cis_h2[1] / cis_h2_se[1])^2) / 2
                    if _args_use_low_rank_approx & (cis_h2_algo == "em_aireml") & (cis_h2_model in ["Ag+Ac", "Ag+Dc", "At+Ac"])
                        _df_hers.cis_rank[i] = approx_ranki
                    end
                    if cis_n_cmp == 3
                        _df_hers.vg2[i] = cis_Σ_i[2]
                        _df_hers.h2_g2[i] = cis_h2[2]
                        _df_hers.h2se_g2[i] = cis_h2_se[2]
                        _df_hers.h2p_g2[i] = ccdf(Chisq(1), (cis_h2[2] / cis_h2_se[2])^2) / 2
                    end
                    _df_hers.cis_snps[i] = n_cis_snps
                    if cis_h2_model in ["At+Ac"]
                        _df_hers.trans_snps[i] = n_trans_snps
                    end
                    if _args_export_rand_eff
                        _df_u_hat[findall(gene .== _df_u_hat.pheno_id), 4:end] .= cis_vc[:u]'
                    end
                end
            end
            if _args_trans_h2
                @timeit to "Trans heritability" if trans_h2_model in ["At", "Dt"]
                    if trans_h2_algo == "minmax"
                        trans_vc = getMRVCModel(trans_KSs, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                    elseif trans_h2_algo == "em_aireml"
                        @notimeit to trans_vc = reml(trans_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, trans_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=trans_EA, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    elseif trans_h2_algo == "idul_aireml"
                        prior_trans_vc = get_IDUL_VarianceComponent(trans_EA, exppheno, X_MME, xQ=trans_xQ, max_iter=5, thre=FloatT(1e-3))
                        @notimeit to trans_vc = reml(trans_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, trans_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, EA=trans_EA, varcmp=[prior_trans_vc[:ΣG][1], prior_trans_vc[:Σe]], em_update=true, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    elseif trans_h2_algo == "idul"
                        trans_vc = get_IDUL_VarianceComponent(trans_EA, exppheno, X_MME, xQ=trans_xQ)
                    end
                    trans_Σ_i .= [trans_vc[:ΣG][1], trans_vc[:Σe]]
                elseif trans_h2_model in ["Ag+Dt"]
                    if trans_h2_algo == "em_aireml"
                        @notimeit to trans_vc = reml(trans_KSs, exppheno, X_MME, _Vi_X, _Xt_Vi_X_i, trans_Hi, _Py, _P, _PA, _Vi, _cVi, _Di, _EAvec_Di, _EAvec; pred_rand_eff=_args_export_rand_eff, reml_max_iter=_args_reml_max_iter, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
                    end
                    trans_Σ_i .= [trans_vc[:ΣG][1], trans_vc[:ΣG][2], trans_vc[:Σe]]
                end
                @runif trans_h2_algo == "minmax" println(trans_Σ_i)
                trans_h2 .= trans_vc[:h2]
                trans_h2_se .= trans_vc[:h2_se]
                _df_hers.vg1[i] = trans_Σ_i[1]
                _df_hers.ve[i] = trans_Σ_i[trans_n_cmp]
                _df_hers.h2_g1[i] = trans_h2[1]
                _df_hers.h2se_g1[i] = trans_h2_se[1]
                _df_hers.h2p_g1[i] = ccdf(Chisq(1), (trans_h2[1] / trans_h2_se[1])^2) / 2
                if trans_n_cmp == 3
                    _df_hers.vg2[i] = trans_Σ_i[2]
                    _df_hers.h2_g2[i] = trans_h2[2]
                    _df_hers.h2se_g2[i] = trans_h2_se[2]
                    _df_hers.h2p_g2[i] = ccdf(Chisq(1), (trans_h2[2] / trans_h2_se[2])^2) / 2
                end
                _df_hers.trans_snps[i] = n_trans_snps
                if h2_algo in ["em_aireml", "idul_aireml", "minmax"]
                    _df_hers.converged[i] = trans_vc[:converged]
                end
                if _args_export_rand_eff
                    _df_u_hat[findall(gene .== _df_u_hat.pheno_id), 4:end] .= trans_vc[:u]'
                end
            end
            if h2_model == "undef"
                multi_vc = getMRVCModel(multi_KSs, exppheno; X=X_MME, maxiter=_args_mm_max_iter, verbose=_args_verbose)
                multi_Σ_i .= [multi_vc[:ΣG]; multi_vc[:Σe]]
                _df_hers[i, [string.("vg", 1:n_grm); "ve"]] .= multi_Σ_i
                _df_hers[i, string.("h2_g", 1:n_grm)] .= multi_vc[:h2]
                _df_hers[i, string.("h2se_g", 1:n_grm)] .= multi_vc[:h2_se]
                _df_hers[i, string.("h2p_g", 1:n_grm)] .= ccdf(Chisq(1), (multi_vc[:h2] ./ multi_vc[:h2_se]) .^ 2) / 2
                _df_hers.converged[i] = multi_vc[:converged]
                println(multi_Σ_i)
            end
            @runif _args_debug if i % 10 == 0
                println_to_file(string("*** Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
                show(to)
            end
        end
        @timeit to "Output results" begin
            cols_float = [_df_hers[1, x] isa AbstractFloat for x in range(1, ncol(_df_hers))]
            _df_hers[:, cols_float] = round.(_df_hers[:, cols_float], sigdigits=6)
            if length(_args_chrom) > 0
                CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".her_est.txt")), _df_hers, delim="\t")
                if _args_glo_h2
                    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".her_est.", chrom, ".B")), _df_covar_B, delim="\t")
                end
                if _args_export_rand_eff
                    cols_float = [_df_u_hat[1, x] isa AbstractFloat for x in range(1, ncol(_df_u_hat))]
                    _df_u_hat[:, cols_float] = round.(_df_u_hat[:, cols_float], sigdigits=6)
                    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".rand_eff.", chrom, ".gz")), _df_u_hat, delim="\t", compress=true, append=is_append)
                end
            else
                is_append = chrom != chroms[1]
                CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".her_est.txt.gz")), _df_hers, delim="\t", compress=true, append=is_append)
                if _args_glo_h2
                    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".her_est.B.gz")), _df_covar_B, delim="\t", compress=true, append=is_append)
                end
                if _args_export_rand_eff
                    cols_float = [_df_u_hat[1, x] isa AbstractFloat for x in range(1, ncol(_df_u_hat))]
                    _df_u_hat[:, cols_float] = round.(_df_u_hat[:, cols_float], sigdigits=6)
                    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".rand_eff.gz")), _df_u_hat, delim="\t", compress=true, append=is_append)
                end
            end
        end
        println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
        if _args_debug
            println_to_file(string(to), log_file)
        end
    end
end
