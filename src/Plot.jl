function runOmiGA_plot(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM::Union{Nothing,Dominance}=nothing, _struct_ITERM::Union{Nothing,InteractionTerm}=nothing)
    time_start = now()
    phenotype = _struct_PHENO.phenotype
    pheno_annotation = _struct_PHENO.annotation
    X_MME = _struct_COVAR.X_MME
    snp_annotation = _struct_GENO.annotation
    genotype = _struct_GENO.genotype
    df_pv_pairs = CSV.File(_args_pvpair_file) |> DataFrame
    if issubset(["pval_g1"], names(df_pv_pairs))
        df_pv_pairs = df_pv_pairs[.!isnan.(df_pv_pairs.pval_g1),:]
    end
    req_colnames = ["pheno_id","variant_id"]
    if !issubset(req_colnames, names(df_pv_pairs))
        error(req_colnames, " are required for `--mode plot`")
    end
    variant_subannot = subset(snp_annotation, :variant => x -> x .∈ Ref(Set(df_pv_pairs.variant_id))) 
    pheno_subannot = subset(pheno_annotation, :pheno_id => x -> x .∈ Ref(Set(df_pv_pairs.pheno_id))) 
    df_pv_pairs.pheno_reindex = vmatch(pheno_subannot.pheno_id, df_pv_pairs.pheno_id)
    df_pv_pairs.variant_reindex = vmatch(variant_subannot.variant, df_pv_pairs.variant_id)
    println_to_file(string("\t* Extracting phenotype vs. genotype data ..."), log_file)
    subGENO = genotype[:,variant_subannot.index]
    transform!(variant_subannot, [:a1, :a1] => ByRow((x, y) -> join([x, y], "/")) => :a1a1) 
    transform!(variant_subannot, [:a1, :a2] => ByRow((x, y) -> join([x, y], "/")) => :a1a2) 
    transform!(variant_subannot, [:a2, :a2] => ByRow((x, y) -> join([x, y], "/")) => :a2a2) 
    variant_gt = variant_subannot[:,[:a2a2,:a1a2,:a1a1]] |> Matrix
    subGENO_str = similar(subGENO,String)
    for x in axes(subGENO, 2)
        subGENO_str[:,x] = variant_gt[x,subGENO[:,x] .+ 1]
    end
    Xt_X = X_MME' * X_MME
    subPHENO = phenotype[:,pheno_subannot.index]
    subPHENO_residuals = get_lm_residuals(X_MME, subPHENO, Xt_X)
    output_type = "combine"
    if output_type == "combine"
        df_plot = repeat(df_pv_pairs[:, ["pheno_id", "variant_id"]], inner=_struct_PHENO.n_samples)
        df_plot.pheno_val = vec(subPHENO[:, df_pv_pairs.pheno_reindex])
        df_plot.pheno_adj = vec(subPHENO_residuals[:, df_pv_pairs.pheno_reindex])
        df_plot.gt_base = vec(subGENO_str[:, df_pv_pairs.variant_reindex])
        df_plot.gt_dosage = vec(subGENO[:, df_pv_pairs.variant_reindex])
        cols_float = [df_plot[1, x] isa AbstractFloat for x in range(1, ncol(df_plot))]
        df_plot[:, cols_float] = round.(df_plot[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".pheno_geno.plot.txt.gz")), df_plot, delim="\t", compress=true)
    else
        subGENO_str_df = DataFrame(subGENO_str, variant_subannot.variant)
        subPHENO_df = DataFrame(round.(subPHENO, sigdigits=6), pheno_subannot.pheno_id)
        subPHENO_residuals_df = DataFrame(round.(subPHENO_residuals, sigdigits=6), pheno_subannot.pheno_id)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".gt.plot.txt.gz")), subGENO_str_df, delim="\t", compress=true)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".pheno_val.plot.txt.gz")), subPHENO_df, delim="\t", compress=true)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".pheno_adj.plot.txt.gz")), subPHENO_residuals_df, delim="\t", compress=true)
    end
    println_to_file(string("\t* DONE."), log_file)
    println_to_file(string("\t* Calculating Pearson's r of variants for local Manhattan plot ..."), log_file)
    chroms = intersect(unique(pheno_subannot.chrom), unique(snp_annotation.chromosome))
    for ch in chroms
        is_append = ch != chroms[1]
        ch_snp_annotation = snp_annotation[snp_annotation.chromosome .== ch,:]
        ch_pheno_annotation = pheno_subannot[pheno_subannot.chrom .== ch,:]
        ch_df_pv_pairs = subset(df_pv_pairs, :pheno_id => x -> x .∈ Ref(Set(ch_pheno_annotation.pheno_id))) 
        _gene_snps = [get_cis_snp_info(ch_pheno_annotation, ch_snp_annotation, gene, _args_cis_window; nsnp_only=true) for gene in ch_df_pv_pairs.pheno_id]
        _gene_end_index = accumulate(+, _gene_snps)
        _gene_start_index = _gene_end_index[1:end-1] .+ 1
        prepend!(_gene_start_index, 1)
        _df_snpr = DataFrame(:pheno_id => repeat([""],sum(_gene_snps)),:variant_id => repeat([""],sum(_gene_snps)),:chrom => repeat([""],sum(_gene_snps)),:pos => repeat([0],sum(_gene_snps)),:corr => repeat([NAN],sum(_gene_snps)))
        _n_pairs = size(ch_df_pv_pairs,1) 
        for i in 1:_n_pairs
            gene = ch_df_pv_pairs.pheno_id[i]
            variant = ch_df_pv_pairs.variant_id[i]
            cissnps_annot, n_cis_snps = get_cis_snp_info(ch_pheno_annotation, ch_snp_annotation, gene, _args_cis_window)
            variant_index = cissnps_annot.index[findfirst(cissnps_annot.variant .== variant)]
            snp_rs = vec(cor(genotype[:,cissnps_annot.index], genotype[:,variant_index]))
            snp_rs .= round.(snp_rs, sigdigits=3)
            _df_snpr[_gene_start_index[i]:_gene_end_index[i],:pheno_id] .= gene
            _df_snpr[_gene_start_index[i]:_gene_end_index[i],:variant_id] .= variant
            _df_snpr[_gene_start_index[i]:_gene_end_index[i],:chrom] .= cissnps_annot.chromosome
            _df_snpr[_gene_start_index[i]:_gene_end_index[i],:pos] .= cissnps_annot.position
            _df_snpr[_gene_start_index[i]:_gene_end_index[i],:corr] .= snp_rs
        end
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".snpr.plot.txt.gz")), _df_snpr, delim="\t", compress=true, append=is_append)
    end
    println_to_file(string("\t* DONE."), log_file)
    println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
end
