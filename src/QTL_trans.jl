function runOmiGA_trans(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM::Union{Nothing,Dominance}=nothing)
    _n, _X_c = size(_struct_COVAR.X_MME)
    time_start = now()
    runOmiGA_gwas(_struct_PHENO, _struct_GENO, _struct_KIN, _struct_COVAR, _struct_PRIOR_HerB; _struct_DOM=_struct_DOM, is_exclude_cis_chrom=true)
    if length(_args_multi_task) == 2
        task_id = _args_multi_task[2]
        task_num = _args_multi_task[1]
        task_collects = collect(Iterators.partition(1:_struct_PHENO.n_phenotypes,ceil(Int,_struct_PHENO.n_phenotypes / task_num)))
        gene_annotation = _struct_PHENO.annotation[task_collects[task_id],:]
    else
        gene_annotation = _struct_PHENO.annotation
    end
    snp_annotation = _struct_GENO.annotation
    df_full = DataFrame()
    gene_annotation.index = 1:nrow(gene_annotation)
    _n_genes = length(gene_annotation.pheno_id) 
    chroms = string.(unique(snp_annotation.chromosome))
    if length(_args_multi_task) == 2
        out_text = joinpath(_args_output_dir,string(_args_out_prefix,".trans_qtl_pairs.task_",_args_multi_task[2],".txt.gz"))
    else
        out_text = joinpath(_args_output_dir,string(_args_out_prefix,".trans_qtl_pairs.txt.gz"))
    end
    _task_write_full_pairs = @task @info "Using Tasks"
    write_header = true
    is_append = false
    for chrom in chroms
        snp_index_chrom = snp_annotation.chromosome .== chrom
        _snp_annot = snp_annotation[snp_index_chrom,:]
        if length(_args_multi_task) == 2
            out_jld = joinpath(_args_output_dir,string(_args_out_prefix,".assoc_pairs.",chrom,".task_",_args_multi_task[2],".jld2"))
        else
            out_jld = joinpath(_args_output_dir,string(_args_out_prefix,".assoc_pairs.",chrom,".jld2"))
        end
        chrom_jld = load(out_jld)["data"]
        nonzero_rows, nonzero_cols, nonzero_vals = findnz(chrom_jld)
        if length(nonzero_vals) == 0
            continue
        end
        snp_index = nonzero_rows[nonzero_cols .∈ (1:3:size(chrom_jld,2),)]
        gene_index = Int.((nonzero_cols[nonzero_cols .∈ (1:3:size(chrom_jld,2),)] .- 1) ./ 3 .+ 1)
        _df_full = DataFrame(
            pheno_id = gene_annotation.pheno_id[gene_index],
            variant_id = _snp_annot.variant[snp_index], 
            af=_snp_annot.af[snp_index],
            beta_g1 = nonzero_vals[nonzero_cols .∈ (1:3:size(chrom_jld,2),)],
            beta_se_g1 = nonzero_vals[nonzero_cols .∈ (2:3:size(chrom_jld,2),)], 
            pval_g1 = ccdf(WaldTest(1, _n - _X_c - 1), nonzero_vals[nonzero_cols .∈ (3:3:size(chrom_jld,2),)])
        )
        append!(df_full,_df_full)
        cols_float = [_df_full[1,x] isa AbstractFloat for x in range(1,ncol(_df_full))]
        _df_full[:,cols_float] = round.(_df_full[:,cols_float],sigdigits=6)
        _df_full.af .= round.(_df_full.af, sigdigits=3)
        CSV.write(out_text,_df_full,delim="\t",compress=true,append=is_append,writeheader=write_header)
        write_header = false
        is_append = true
    end
    if isfile(out_text)
        if length(_args_multi_task) == 2
            rm.(glob(string(_args_out_prefix,".assoc_pairs.*.task_",_args_multi_task[2],".jld2"),_args_output_dir))
        else
            rm.(glob(string(_args_out_prefix,".assoc_pairs.*.jld2"),_args_output_dir))
        end
    end
    println_to_file(string("+++ Total elapsed time (h:m:s:ms): ",format_milliseconds(now()-time_start)), log_file)
    if _args_debug
        println_to_file(string(to), log_file)
        println()
    end
end
