function get_parsed_args(args)
    settings = ArgParseSettings()
    settings.version = string("OmiGA v", relased_omiga_version)
    settings.add_version = true
    settings.prog = "OmiGA: a toolkit for Omics Genetic Analysis"
    settings.usage = "OmiGA: a toolkit for Omics Genetic Analysis\nVersion: " * relased_omiga_version
    settings.description = "****************************************************************\n* Help pages at https://omiga.bio/\n* Please report bugs to Jinyan Teng <jinyan.teng@scau.edu.cn>\n****************************************************************\n"
    settings.preformatted_description = true
    settings.help_width = 120
    add_arg_group!(settings, "required arguments")
    @add_arg_table! settings begin
        "--prefix"
            help = "prefix for output file names"
            arg_type = String
        "--mode"
            help = "mode selection"
            arg_type = String
    end;
    add_arg_group!(settings, "input arguments")
    @add_arg_table! settings begin
        "--genotype", "--geno"
            help = "genotypes in PLINK format"
            arg_type = String
        "--phenotype", "--pheno"
            help = "phenotypes file"
            arg_type = String
            nargs = '*'
        "--pheno-annot"
            help = "annotation file for phenotypes"
            arg_type = String
        "--pheno-group"
            help = "group file for phenotypes"
            arg_type = String
        "--pheno-pairs"
            help = "non-definition"
            arg_type = String
        "--id-map"
            help = "map file for sample ID"
            arg_type = String
        "--covariates", "--covariate"
            help = "covariates file"
            arg_type = String
            nargs = '*'
        "--covar-transpose"
            help = "transposed covariates file"
            action = :store_true
        "--extract-covar-name"
            help = "covariates name list extracted"
            arg_type = String
            nargs = '*'
        "--exclude-covar-name"
            help = "covariates name list excluded"
            arg_type = String
            nargs = '*'
        "--dcovar-name"
            help = "name of discrete covariates"
            arg_type = String
            nargs = '*'
    end;
    add_arg_group!(settings, "optional arguments")
    @add_arg_table settings begin
        "--phen-pc-covar"
            help = "number of phenotype PCs"
            arg_type = Int
        "--geno-pc-covar"
            help = "number of genotype PCs"
            arg_type = Int
        "--dprop-pc-covar"
            help = "parameter for automatic elbow detection"
            arg_type = Float64
            nargs = '*'
        "--PCAForQTL", "--pca-for-qtl"
            help = "method of PCAForQTL"
            arg_type = String
        "--rm-collinear-covar"
            help = "collinear threshold (absolute correlation) of covariates"
            arg_type = Float64
            default = 0.95
        "--vif-threshold"
            help = "collinear threshold (VIF) of covariates"
            arg_type = Float64
            default = 10.0
        "--no-intercept"
            help = "do not use the intercept term"
            action = :store_true
        "--normalized-covar"
            help = "normalized the covariates"
            action = :store_true
        "--interaction"
            help = "interaction file"
            arg_type = String
        "--interaction-name"
            help = "name of interaction term used"
            arg_type = String
        "--no-add-interaction-to-covar"
            help = "do not use interaction term as covariate"
            action = :store_true
        "--cis-file"
            help = "output file(s) from '--mode cis'"
            nargs = '*'
            arg_type = String
        "--cis-qtl-pairs"
            help = "output file(s) from '--mode cis'"
            nargs = '*'
            arg_type = String
        "--her-file"
            help = "output file(s) from '--mode her_est'"
            nargs = '*'
            arg_type = String
        "--multiple-testing", "--multiple-testing-first"
            help = "method for multiple testing"
            arg_type = String
        "--multiple-testing-second"
            help = "method for multiple testing"
            arg_type = String
        "--storey-lambda"
            help = "lambda for Storey's qvalues"
            arg_type = Float64
            default = 0.1 
        "--permutations"
            help = "number of permutations"
            arg_type = Int
        "--calcu-variant-threshold"
            help = "compute variant-level nominal P-value threshold"
            action = :store_true
        "--extract-significant-qtls"
            help = "extract vairant-phenotype pairs bellow the P-value threshold"
            action = :store_true
        "--maf-match"
            help = "Approximate range of MAF matched for target variants"
            arg_type = Float64
        "--ld-match"
            help = "Approximate s.d. of LD matched for target variants"
            arg_type = Float64
        "--ldscore-file"
            help = "LD score file"
            arg_type = String
        "--enrich-target"
            help = "bed file(s) of target set used for enrichment analysis"
            arg_type = String
        "--enrich-annot"
            help = "bed file(s) of annotation used for enrichment analysis"
            nargs = '*'
            arg_type = String
        "--dtss-weight"
            help = "dTSS weight for calculate ACAT P-values"
            action = :store_true
        "--nominal-only"
            help = "permutation-free"
            action = :store_true
        "--base-index"
            help = "0-based or 1-based genomic coordinate"
            arg_type = Int
        "--cis-window"
            help = "cis-window size in bp"
            arg_type = Int
            default = 1000000
        "--trans-window"
            help = "trans-window size in bp"
            arg_type = Int
        "--window-type"
            help = "type of cis-region definition"
            arg_type = String
            default = "tss"
        "--mac-threshold"
            help = "minor allele count (MAC) threshold for genotypes"
            arg_type = Int
            default = 6
        "--maf-threshold"
            help = "minor allele frequency (MAF) threshold for genotypes"
            arg_type = Float64
            default = 0.01
        "--maf-threshold-interaction"
            help = "MAF threshold for '--mode cis_interaction'"
            arg_type = Float64
            default = 0.05
        "--het-threshold-interaction"
            help = "Her rate threshold for '--mode cis_interaction'"
            arg_type = Float64
            default = 0.99
        "--nofilter"
            help = "mute this default filtering"
            action = :store_true
        "--normalized-phenotype"
            help = "inverse normal transformed the phenotypic values"
            action = :store_true
        "--normalized-interaction"
            help = "inverse normal transformed the interaction term"
            action = :store_true
        "--write-top"
            help = "write only summary statistics of top variants"
            action = :store_true
        "--output-format"
            help = "the output format of the full summary statistics files"
            arg_type = String
        "--fdr"
            help = "false discovery rate (FDR) for QTL discovery"
            arg_type = Float64
            default = 0.05
        "--seed"
            help = "random seed"
            arg_type = Int
            default = 666
        "--chrom"
            help = "chromosome(s) used for analysis"
            arg_type = String
            nargs = '*'
        "--output-dir", "-O"
            help = "output directory"
            arg_type = String
            default = "."
        "--sparse-grm"
            help = "cut-off threshold for sparse GRM"
            arg_type = Float64
        "--grm"
            help = "prefix of GRM file(s)"
            arg_type = String
            nargs = '*'
        "--regional-grm"
            help = "regional file to build multiple GRMs"
            arg_type = String
        "--grm-alpha"
            help = "alpha for compute GRM"
            arg_type = Int
            default = 0
        "--grm-snps"
            help = "number of variants used to compute GRM"
            arg_type = Int
        "--save-grm"
            help = "save GRM to disk"
            action = :store_true
        "--extract-pheno"
            help = "file including phenotypes used for analysis"
            arg_type = String
            nargs = '*'
        "--extract-pheno-name"
            help = "name of phenotypes used for analysis"
            arg_type = String
            nargs = '*'
        "--exclude-pheno"
            help = "file including phenotypes excluded for analysis"
            arg_type = String
        "--exclude-pheno-name"
            help = "name of phenotypes excluded for analysis"
            arg_type = String
            nargs = '*'
        "--dpars"
            help = "converge parameter of AIREML"
            arg_type = Float64
            default = 1e-6
        "--enable-lrt"
            help = "enable LRT for heritability"
            action = :store_true
        "--h2-model"
            help = "statistical model for heritability estimation"
            arg_type = String
        "--h2-algo"
            help = "algorithm for heritability estimation"
            arg_type = String
        "--export-random-eff"
            help = "export the esimates of random effects"
            action = :store_true
        "--qtl-map-algo"
            help = "algorithm for QTL mapping"
            arg_type = String
            default = "idul"
        "--qtl-map-model"
            help = "statistical model for QTL mapping"
            arg_type = String
            default = "a+A" 
        "--exact-map"
            help = "exact QTL mapping"
            action = :store_true
        "--idul-max-iter"
            help = "maximum iterations for IDUL"
            arg_type = Int64
            default = 20
        "--reml-max-iter"
            help = "maximum iterations for REML"
            arg_type = Int64
            default = 100
        "--mm-max-iter"
            help = "maximum iterations for MM algorithm"
            arg_type = Int64
            default = 1000
        "--idul-converge"
            help = "converge threshold for IDUL"
            arg_type = Float64
            default = 1e-6
        "--inv-precision"
            help = "add a small value to the diagonal of XtX"
            arg_type = Float64
            default = 1e-10
        "--eigen-algo"
            help = "algorithm used for eigen decomposition"
            arg_type = String
            default = "gesdd"
        "--wald-test"
            help = "distribution (FDist/Chisq) used for Wald test"
            arg_type = String
            default = "FDist"
        "--est-ind-h2"
            help = "perform SNP heritability estimation for '--mode cis_independent'"
            action = :store_true
        "--ind-r2"
            help = "r2 threshold for independent cis-QTL mapping"
            arg_type = Float64
            default = 0.5
        "--paed-step"
            help = "step length for Pre-Allocated Eigen Decomposition (PAED)"
            arg_type = Float64
            default = 0.01
        "--paed-max"
            help = "upper boundary for Pre-Allocated Eigen Decomposition (PAED)"
            arg_type = Float64
            default = 25
        "--paed-file"
            help = "file path (*.paed) of PAED"
            arg_type = String
        "--export-paed-file"
            help = "save PAED file to disk"
            action = :store_true
        "--low-rank-approx"
            help = "use cisLRA for cis-heritability estimation."
            action = :store_true
        "--approx-method"
            help = "method for low-rank matrix approximation"
            arg_type = String
            default = "psvd"
        "--approx-rank"
            help = "approximate rank for cisLRA"
            arg_type = Number
            default = 0.05
        "--max-approx-rank"
            help = "maximum approximate rank for cisLRA"
            arg_type = Number
        "--avg-2norm-threshold"
            help = "threshold of average 2-norm for cisLRA"
            arg_type = Float64
            default = 2e-4
        "--rm-pheno-threshold"
            help = "threshold of median proportion for excluding phenotypes"
            arg_type = Float64
            default = 0.01
        "--het-rate-threshold", "--het-threshold"
            help = "threshold of heterozygote rate for genotype QC"
            arg_type = Float64
            default = 0.99
        "--pval-threshold"
            help = "output only significant phenotype-variant pairs with a P-value below threshold"
            arg_type = Float64
            default = 1e-5
        "--preadj-covar"
            help = "pre-adjusting the phenotypes"
            action = :store_true
        "--lambda-only"
            help = "only compute GC lambda"
            action = :store_true
        "--lambda-snps"
            help = "number of variants used for computing GC lambda"
            arg_type = Int
            default = 100000
        "--lambda-phenos"
            help = "number of phenotypes used for computing GC lambda"
            arg_type = Int
        "--pvpair-file"
            help = "file including phenotype-variant pairs for extracting plot data"
            arg_type = String
        "--maxiter"
            help = "non-definition"
            arg_type = Int
            default = 20
        "--cc-par"
            help = "non-definition"
            arg_type = Float64
            default = 1e-6
        "--cc-gra"
            help = "non-definition"
            arg_type = Float64
            default = 1e-4
        "--emwa"
            help = "non-definition"
            arg_type = Float64
            default = 0.05
        "--opt-stri"
            help = "non-definition"
            arg_type = String
        "--opt-mode"
            help = "non-definition"
            arg_type = String
            nargs = '*'
        "--multi-task"
            help = "split the analysis into multiple independent tasks and run the one of them"
            nargs = 2
            arg_type = Int
        "--chunk-size"
            help = "load genotypes into memory in chunks of approximate 'chunk-size' variants in QTL mapping"
            arg_type = Int
        "--batch-size"
            help = "load phenotypes into memory by batch size"
            arg_type = Int
        "--low-mem"
            help = "memory-saving mode"
            action = :store_true
        "--gpu"
            help = "using GPU"
            action = :store_true
        "--force-double-precision"
            help = "use double precision for computation"
            action = :store_true
        "--use-gzip"
            help = "use gzip to compress output files"
            action = :store_true
        "--threads"
            help = "number of threads"
            arg_type = Int
        "--mkl-threads"
            help = "number of MKL threads"
            arg_type = Int
        "--verbose"
            help = "verbose logs output"
            action = :store_true
        "--debug"
            help = "debug"
            action = :store_true
        "--experimental"
            help = "enable experimental functions"
            action = :store_true
    end
    parsed_args = parse_args(args, settings)
    return parsed_args
end
