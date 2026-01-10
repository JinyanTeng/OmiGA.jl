module OmiGA
using DataFrames
using Statistics
using StatsBase
using CSV
using LinearAlgebra
using Random
using Dates
using SparseArrays
using Glob
using Combinatorics
using Base.Threads
using StatsModels
using GLM 
using IterativeSolvers 
using InvertedIndices, SweepOperator 
using CodecZlib, CodecXz, CodecBzip2, CodecZstd, TranscodingStreams 
using Adapt, Missings, Printf 
using Requires, VectorizationBase, VariantCallFormat 
using ArgParse
using Distributions 
using Optim 
using HypothesisTests 
using SpecialFunctions 
using MultipleTesting 
using TimerOutputs 
using ACAFact 
using LowRankApprox 
using StaticArrays 
using JLD2
using MultivariateStats 
using Mmap 
using OffsetArrays 
using Arrow
using PooledArrays 
using Polyester
using LoopVectorization
using Tables 
using GenomicFeatures 
using IntervalTrees 
using MKL 
using CUDA
using Suppressor
if false
    mat_mul = matmul
    mat_mul! = matmul! 
else
    mat_mul = *
    mat_mul! = mul!
end
global relased_omiga_version = "1.8.8"
const to = TimerOutput()
dir_src = @__DIR__
include(joinpath(dir_src, "Modified_MultiResponseVarianceComponentModels/MultiResponseVarianceComponentModels.jl"))
using .MultiResponseVarianceComponentModels
include(joinpath(dir_src, "Modified_SnpArrays/SnpArrays.jl"))
using .SnpArrays
include(joinpath(dir_src, "SNP_matrix.jl"))
include(joinpath(dir_src, "Utils.jl"))
include(joinpath(dir_src, "ARGS.jl"))
include(joinpath(dir_src, "Freq.jl"))
include(joinpath(dir_src, "Init.jl"))
include(joinpath(dir_src, "QTL_cis.jl"))
include(joinpath(dir_src, "Permutation.jl"))
include(joinpath(dir_src, "QTL_independent.jl"))
include(joinpath(dir_src, "QTL_trans.jl"))
include(joinpath(dir_src, "QTL_global.jl"))
include(joinpath(dir_src, "QTL_multiple_test.jl"))
include(joinpath(dir_src, "Heritability.jl"))
include(joinpath(dir_src, "GRM.jl"))
include(joinpath(dir_src, "Clipper.jl"))
include(joinpath(dir_src, "IDUL.jl"))
include(joinpath(dir_src, "REML.jl"))
include(joinpath(dir_src, "MOLQTL.jl"))
include(joinpath(dir_src, "Plot.jl"))
include(joinpath(dir_src, "Enrichment.jl"))
function (@main)(ARGS)
    try
        real_main(ARGS)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
function real_main(ARGS)::Cint
    MAIN_ARGS = ARGS
    reset_timer!(to)
    if length(MAIN_ARGS) == 0
        println(" * No arguments are provided!")
        return 0
    end
    global parsed_args = get_parsed_args(MAIN_ARGS)
    str_cmd = join(MAIN_ARGS, " ")
    if parsed_args["debug"]
        println(" * Runing OmiGA with debug mode!")
        println(parsed_args)
    end
    begin
        global _args_geno_file_prefix = parsed_args["genotype"]
        global _args_phenotype_files = parsed_args["phenotype"]
        global _args_pheno_annot_file = parsed_args["pheno-annot"]
        global _args_pheno_group_file = parsed_args["pheno-group"]
        global _args_id_map_file = parsed_args["id-map"]
        global _args_out_prefix = parsed_args["prefix"]
        global _args_run_mode = parsed_args["mode"]
        global _args_covar_files = parsed_args["covariates"]
        global _args_covar_transpose = parsed_args["covar-transpose"]
        global _args_extract_covar_name = parsed_args["extract-covar-name"]
        global _args_exclude_covar_name = parsed_args["exclude-covar-name"]
        global _args_dcovar_name = parsed_args["dcovar-name"]
        global _args_rm_collinear_covar = parsed_args["rm-collinear-covar"]
        global _args_vif_threshold = parsed_args["vif-threshold"]
        global _args_normalized_covar = parsed_args["normalized-covar"]
        global _args_normalized_interaction = parsed_args["normalized-interaction"]
        global _args_normalized_phenotype = parsed_args["normalized-phenotype"]
        global _args_phen_pc_covar = parsed_args["phen-pc-covar"]
        global _args_geno_pc_covar = parsed_args["geno-pc-covar"]
        global _args_dprop_pc_covar = parsed_args["dprop-pc-covar"]
        global _args_PCAForQTL = parsed_args["PCAForQTL"]
        global _args_no_intercept = parsed_args["no-intercept"]
        global _args_interaction_file = parsed_args["interaction"]
        global _args_interaction_name = parsed_args["interaction-name"]
        global _args_no_add_interaction_to_covar = parsed_args["no-add-interaction-to-covar"]
        global _args_cis_output = parsed_args["cis-file"] 
        global _args_cis_pairs_output = parsed_args["cis-qtl-pairs"] 
        global _args_her_output = parsed_args["her-file"]
        global _args_multiple_testing = parsed_args["multiple-testing"]
        global _args_n_perms = parsed_args["permutations"]
        global _args_storey_lambda = parsed_args["storey-lambda"]
        global _args_dtss_weight = parsed_args["dtss-weight"]
        global _args_calcu_variant_threshold = parsed_args["calcu-variant-threshold"]
        global _args_output_significant_qtls = parsed_args["extract-significant-qtls"]
        global _args_nominal_only = parsed_args["nominal-only"]
        global _args_base_index = parsed_args["base-index"]
        global _args_cis_window = parsed_args["cis-window"]
        global _args_window_type = parsed_args["window-type"]
        global _args_trans_window = parsed_args["trans-window"]
        global _args_mac_threshold = parsed_args["mac-threshold"]
        global _args_maf_threshold = parsed_args["maf-threshold"]
        global _args_maf_threshold_interaction = parsed_args["maf-threshold-interaction"]
        global _args_het_threshold_interaction = parsed_args["het-threshold-interaction"]
        global _args_nofilter = parsed_args["nofilter"]
        global _args_write_top = parsed_args["write-top"]
        global _args_output_format = parsed_args["output-format"]
        global _args_fdr = parsed_args["fdr"]
        global _args_seed = parsed_args["seed"] 
        global _args_multi_task = parsed_args["multi-task"]
        global _args_chrom = parsed_args["chrom"]
        global _args_output_dir = parsed_args["output-dir"]
        global _args_sparse_grm = parsed_args["sparse-grm"] 
        global _args_regional_grm = parsed_args["regional-grm"] 
        global _args_save_grm = parsed_args["save-grm"] 
        global _args_grm_prefix = parsed_args["grm"]
        global _args_grm_snps = parsed_args["grm-snps"]
        global _args_grm_alpha = parsed_args["grm-alpha"]
        global _args_extract_pheno_files = parsed_args["extract-pheno"]
        global _args_extract_pheno_name = parsed_args["extract-pheno-name"]
        global _args_exclude_pheno = parsed_args["exclude-pheno"]
        global _args_exclude_pheno_name = parsed_args["exclude-pheno-name"]
        global _args_dpars = parsed_args["dpars"] 
        global _args_h2_model = parsed_args["h2-model"] 
        global _args_h2_algo = parsed_args["h2-algo"] 
        global _args_export_rand_eff = parsed_args["export-random-eff"]
        global _args_est_ind_h2 = parsed_args["est-ind-h2"] 
        global _args_ind_r2 = parsed_args["ind-r2"] 
        global _args_rm_pheno_threshold = parsed_args["rm-pheno-threshold"]
        global _args_het_rate_threshold = parsed_args["het-rate-threshold"]
        global _args_qtl_map_algo = parsed_args["qtl-map-algo"] 
        global _args_qtl_map_model = parsed_args["qtl-map-model"] 
        global _args_exact_map = parsed_args["exact-map"]
        global _args_idul_max_iter = parsed_args["idul-max-iter"]
        global _args_reml_max_iter = parsed_args["reml-max-iter"]
        global _args_mm_max_iter = parsed_args["mm-max-iter"]
        global _args_idul_converge = parsed_args["idul-converge"]
        global _args_inv_precision = parsed_args["inv-precision"]
        global _args_wald_test = parsed_args["wald-test"]
        global _args_eigen_algo = parsed_args["eigen-algo"]
        global _args_paed_step = parsed_args["paed-step"]
        global _args_paed_max = parsed_args["paed-max"]
        global _args_path_partialv = parsed_args["paed-file"]
        global _args_export_paed_file = parsed_args["export-paed-file"]
        global _args_preadj_covar = parsed_args["preadj-covar"]
        global _args_pval_threshold = parsed_args["pval-threshold"]
        global _args_lambda_only = parsed_args["lambda-only"]
        global _args_lambda_snps = parsed_args["lambda-snps"]
        global _args_lambda_phenos = parsed_args["lambda-phenos"]
        global _args_use_low_rank_approx = parsed_args["low-rank-approx"]
        global _args_approx_method = parsed_args["approx-method"]
        global _args_approx_rank = parsed_args["approx-rank"]
        global _args_max_approx_rank = parsed_args["max-approx-rank"]
        global _args_avg_2norm_threshold = parsed_args["avg-2norm-threshold"]
        global _args_pheno_pairs = parsed_args["pheno-pairs"]
        global _args_maxiter = parsed_args["maxiter"]
        global _args_cc_par = parsed_args["cc-par"]
        global _args_cc_gra = parsed_args["cc-gra"]
        global _args_emwa = parsed_args["emwa"]
        global _args_opt_stri = parsed_args["opt-stri"]
        global _args_opt_mode = parsed_args["opt-mode"]
        global _args_pvpair_file = parsed_args["pvpair-file"]
        global _args_maf_match = parsed_args["maf-match"]
        global _args_ld_match = parsed_args["ld-match"]
        global _args_ldscore_file = parsed_args["ldscore-file"]
        global _args_enrich_target = parsed_args["enrich-target"]
        global _args_enrich_annot = parsed_args["enrich-annot"]
        global _args_gpu = parsed_args["gpu"]
        global _args_chunk_size = parsed_args["chunk-size"]
        global _args_use_gzip = parsed_args["use-gzip"]
        global _args_low_mem = parsed_args["low-mem"]
        global _args_force_double_precision = parsed_args["force-double-precision"]
        global _args_threads = isnothing(parsed_args["threads"]) ? nthreads() : parsed_args["threads"]
        global _args_mkl_threads = parsed_args["mkl-threads"]
        global _args_verbose = parsed_args["verbose"]
        global _args_debug = parsed_args["debug"]
        global _args_experimental = parsed_args["experimental"]
        @runif _args_debug ENV["JULIA_DEBUG"] = Main
        global _args_low_mem = (isnothing(_args_chunk_size) ? false : true) || _args_low_mem
        global getG = _args_low_mem ? getG_fast : getG_fast
        if isnothing(_args_mkl_threads)
            threads_mkl = _args_threads
        else
            threads_mkl = _args_mkl_threads
        end
        BLAS.set_num_threads(threads_mkl) 
        @runif _args_debug println(string("Assigned MKL threads: ", threads_mkl))
        @runif _args_debug println(string("Used MKL threads: ", BLAS.get_num_threads()))
        Random.seed!(_args_seed)
        global _args_linear_model = _args_qtl_map_model in ["a", "d", "a+d"]
        global gxe_omit_main_eff = false
        global is_pre_adj_genotype = false 
        global is_sample_byrow = true 
        if (_args_run_mode in ["cis", "cis_interaction", "cis_independent"]) & (_args_qtl_map_algo != "idul")
            global is_sample_byrow = false 
        end
        if _args_run_mode in ["cis", "cis_interaction", "cis_independent", "trans", "gwas", "her_est", "birg", "xbirg"]
            global FloatT = Float32
        else
            global FloatT = Float64
        end
        if _args_run_mode in ["cis_independent"]
            global FloatT2 = Float64
        else
            global FloatT2 = Float32
        end
        if (_args_run_mode == "her_est") & !isnothing(_args_h2_algo)
            @runif _args_h2_algo == "minmax" global FloatT = Float64
        end
        @runif _args_force_double_precision global FloatT = Float64
        @runif _args_force_double_precision global FloatT2 = Float64
        global USE_Float32 = FloatT == Float32
        global NAN = USE_Float32 ? NaN32 : NaN
        global USE_GPU = _args_gpu
        @runif USE_GPU begin
            global USE_GPU = CUDA.functional()
        end
        if isnothing(_args_out_prefix)
            error(string("The '--prefix' option need to be specified."))
        end
    end
    global log_file = create_OUTDIR_and_LOG(_args_output_dir, _args_run_mode, _args_chrom, _args_out_prefix)
    println_to_file(string("****************************************************************"), log_file)
    println_to_file(string("* OmiGA: a toolkit for Omics Genetic Analysis"), log_file)
    println_to_file(string("* Version: ", relased_omiga_version), log_file)
    println_to_file(string("* Help Pages: https://omiga.bio/"), log_file)
    println_to_file(string("* Please report bugs to Jinyan Teng <jinyan.teng@scau.edu.cn>"), log_file)
    println_to_file(string("* or post a disscusion in OmiGA BBS (https://bbs.omiga.bio/)"), log_file)
    println_to_file(string("****************************************************************"), log_file)
    println_to_file("", log_file)
    println_to_file("Command Line Input:", log_file)
    println_to_file(replace(replace(str_cmd, " --" => " \\\n--"), "--" => "  --"), log_file) 
    println_to_file("", log_file)
    @runif _args_debug println_to_file(string(BLAS.get_config()), log_file)
    println_to_file(string("Hostname: ", gethostname(), ", CPU.threads: ", nthreads(), ", ", match(r"mkl", string(BLAS.get_config())) != nothing ? "MKL" : "BLAS", ".threads: ", BLAS.get_num_threads()), log_file)
    println_to_file("", log_file)
    time_start_omiga = now()
    println_to_file(string("Program start at: ", time_start_omiga), log_file)
    if !(_args_run_mode in ["cis_mt", "enrich"])
        if _args_run_mode != "xbirg"
            @time _struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_ITERM, _struct_PRIOR_HerB, _struct_DOM = prepare_Data()
        else
            _struct_PHENO, _struct_PHENO_2, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_COVAR_2 = prepare_Data()
        end
    end
    if _args_run_mode == "her_est"
        runOmiGA_her(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "cis"
        runOmiGA_cis(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "cis_mt"
        runOmiGA_cis_mt()
    elseif _args_run_mode == "cis_independent"
        runOmiGA_independent(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "cis_interaction"
        runOmiGA_cis(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM, _struct_ITERM=_struct_ITERM)
    elseif _args_run_mode == "trans"
        runOmiGA_trans(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "gwas"
        runOmiGA_gwas(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "plot"
        runOmiGA_plot(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM)
    elseif _args_run_mode == "birg"
        runOmiGA_birg(_struct_PHENO, _struct_KIN, _struct_COVAR)
    elseif _args_run_mode == "xbirg"
        runOmiGA_xbirg(_struct_PHENO, _struct_PHENO_2, _struct_KIN, _struct_COVAR, _struct_COVAR_2)
    elseif _args_run_mode == "enrich"
        runOmiGA_enrich(_args_enrich_target, _args_enrich_annot; bkg_plink_prefix=_args_geno_file_prefix, ldscore_file=_args_ldscore_file, maf_match=_args_maf_match, ld_match=_args_ld_match, threshold=nothing, n_perms=_args_n_perms)
    end
    time_end_omiga = now()
    println_to_file(string("Program end at: ", time_end_omiga), log_file)
    println_to_file(string("Elapsed time (h:m:s:ms): ", format_milliseconds(time_end_omiga - time_start_omiga)), log_file)
    return 0
end
if abspath(PROGRAM_FILE) == @__FILE__
    real_main(ARGS)
end
end 
