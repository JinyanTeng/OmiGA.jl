function runOmiGA_cis_mt()
    time_start = now()
    println_to_file(string("  cis_mt mode start at: ", time_start), log_file)
    df_tops = DataFrame()
    df_tops_detail = nothing
    if isfile(_args_cis_output[1])
        file_cis_qtl_txt = _args_cis_output
        df_tops = CSV.File(file_cis_qtl_txt, header=true, buffer_in_memory=true) |> DataFrame
        if issubset(["group_id"], names(df_tops))
            df_tops_detail = CSV.File(replace.(file_cis_qtl_txt,".txt.gz"=>".detail.txt.gz"), header=true, buffer_in_memory=true) |> DataFrame
        end
    else
        error("Please provide correct file for `--cis-file`!")
    end
    with_group = !isnothing(df_tops_detail)
    n_tests = any(.!isnothing.(match.(r"g2", names(df_tops)))) ? 2 : 1
    cis_summary_prefix = replace(_args_cis_output[1], ".cis_qtl.txt.gz" => "")
    summary_dir = dirname(cis_summary_prefix)
    file_prefix = basename(cis_summary_prefix)
    if _args_output_significant_qtls
        significant_phenotypes = df_tops.pheno_id[df_tops[:, ["qval_g1", "qval_g2"][n_tests]].<_args_fdr]
        all_files = readdir(summary_dir)
        full_summary_files = all_files[.!isnothing.(match.(Regex(string(file_prefix, ".cis_qtl_pairs")), all_files))]
        println_to_file(string("Found ", length(full_summary_files), " *.cis_qtl_pairs.txt.gz files."), log_file)
        if length(full_summary_files) == 0
            error("*.cis_qtl_pairs.txt.gz files must be provided the same dir as ", basename(_args_cis_output[1]))
        else
            sorted_indices = sortperm(replace.(full_summary_files, ".txt.gz" => "", r"(.)*.cis_qtl_pairs." => ""), lt=compare_chromosomes)
            full_summary_files = full_summary_files[sorted_indices]
            n_significant_qtls = 0
            out_file = joinpath(_args_output_dir, string(_args_out_prefix, ".significant_qtl_pairs.txt.gz")) 
            for fn in full_summary_files
                _df_sig_qtls_df = DataFrame()
                println_to_file(string("    * Processing ", fn), log_file)
                _cis_qtl_pairs = CSV.File(joinpath(summary_dir, fn)) |> DataFrame
                _cis_qtl_pairs = _cis_qtl_pairs[_cis_qtl_pairs[:, ["pval_g1", "pval_g2"][n_tests]].<=maximum(df_tops[:, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]), :]
                pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                pheno_ids = pheno_ids[findall(pheno_ids .∈ (significant_phenotypes,))]
                for _pid in pheno_ids
                    _pval_mt = _cis_qtl_pairs[_cis_qtl_pairs.pheno_id.==_pid, :]
                    append!(_df_sig_qtls_df, _pval_mt[_pval_mt[:, ["pval_g1", "pval_g2"][n_tests]].<=df_tops[findfirst(df_tops.pheno_id .== _pid), ["pval_g1_threshold", "pval_g2_threshold"][n_tests]], :])
                end
                n_significant_qtls += size(_df_sig_qtls_df, 1)
                CSV.write(out_file, _df_sig_qtls_df, delim="\t", compress=true, append=fn != full_summary_files[1])
            end
            println_to_file(string(" *** Found ", n_significant_qtls, " significant variant-phenotype pairs."), log_file)
        end
        return 0
    end
    if isnothing(_args_multiple_testing)
        error("--multiple-testing must be specified!")
    else
        multiple_testing_method = _args_multiple_testing
    end
    storey_lambda = _args_storey_lambda
    out_file_middle_part = multiple_testing_method
    df_tops_mt = copy(df_tops)
    df_tops_detail_mt = !with_group ? nothing : df_tops_detail
    chroms = nothing
    if issubset(["chrom"], names(df_tops_mt))
        chroms = unique(df_tops_mt.chrom)
    end
    mt_method_list = Dict("acat" => "ACAT", "storey" => "Storey", "bh" => "BenjaminiHochberg", "bonf" => "Bonferroni", "by" => "BenjaminiYekutieli", "clipper" => "ClipperQTL", "beta_approx" => "Beta approximation")
    if multiple_testing_method in ["acat", "storey", "bh", "bonf", "by"]
        println_to_file(string("    * Re-perform multiple-testing correction: first-level correction -> ", mt_method_list[multiple_testing_method], ", second-level correction -> ", mt_method_list["bh"]), log_file)
        if multiple_testing_method == "storey"
            println_to_file(string("    * lambda=",storey_lambda," was used for Storey."), log_file)
            out_file_middle_part = string(multiple_testing_method, "_", storey_lambda)
        end
        pval_colname = ["pval_g1", "pval_g2"][n_tests]
        rm_colnames_list = pval_colname * "_" .* ["acat", "storey", "bh", "bonf", "by"]
        pval_adj_colname = ["pval_g1_" * multiple_testing_method, "pval_g2_" * multiple_testing_method][n_tests]
        qval_colname = ["qval_g1", "qval_g2"][n_tests]
        rm_colnames = rm_colnames_list[rm_colnames_list .∈ (names(df_tops_mt),)]
        @runif length(rm_colnames) > 0 select!(df_tops_mt, Not(rm_colnames))
        @runif issubset([qval_colname], names(df_tops_mt)) select!(df_tops_mt, Not(qval_colname))
        df_tops_mt[:, pval_adj_colname] .= NaN
        df_tops_mt[:, qval_colname] .= NaN
        if with_group
            rm_colnames = rm_colnames_list[rm_colnames_list .∈ (names(df_tops_detail_mt),)]
            @runif length(rm_colnames) > 0 select!(df_tops_detail_mt, Not(rm_colnames))
            df_tops_detail_mt[:, pval_adj_colname] .= NaN
        end
        all_files = readdir(summary_dir)
        full_summary_files = all_files[.!isnothing.(match.(Regex(string(file_prefix, ".cis_qtl_pairs")), all_files))]
        println_to_file(string("    ** Found ", length(full_summary_files), " *.cis_qtl_pairs.txt.gz files"), log_file)
        if length(full_summary_files) == 0
            error("*.cis_qtl_pairs.txt.gz files must be provided the same dir as ", basename(_args_cis_output[1]))
        else
            sorted_indices = sortperm(replace.(full_summary_files, ".txt.gz" => "", r"(.)*.cis_qtl_pairs." => ""), lt=compare_chromosomes) 
            full_summary_files = full_summary_files[sorted_indices]
            if isnothing(chroms)
                chroms = replace.(full_summary_files, ".txt.gz" => "", r"(.)*.cis_qtl_pairs." => "")
            end
            for fn in full_summary_files
                println_to_file(string("    *** Processing ", fn), log_file)
                _cis_qtl_pairs = CSV.File(joinpath(summary_dir, fn)) |> DataFrame
                pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                group_ids = nothing
                @runif with_group group_ids = unique(df_tops_detail_mt.group_id[df_tops_detail_mt.pheno_id .∈ (pheno_ids,)])
                println_to_file(string("        Number of phenotypes: ", length(pheno_ids)), log_file)
                if !with_group
                    @threads for _pid in pheno_ids
                        _pval_mt = _cis_qtl_pairs[_cis_qtl_pairs.pheno_id.==_pid, pval_colname]
                        pval_adj = NaN
                        _pval_mt = _pval_mt[.!isnan.(_pval_mt)]
                        @runif length(_pval_mt)>0 if multiple_testing_method == "acat"
                            pval_adj = ACATest(_pval_mt; is_check=false)
                        elseif multiple_testing_method == "storey"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochbergAdaptive(Storey(storey_lambda)))) 
                        elseif multiple_testing_method == "bh"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochberg()))
                        elseif multiple_testing_method == "bonf"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, Bonferroni()))
                        elseif multiple_testing_method == "by"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiYekutieli()))
                        end
                        df_tops_mt[df_tops_mt.pheno_id.==_pid, pval_adj_colname] .= pval_adj
                    end
                else
                    println_to_file(string("        Number of groups: ", length(group_ids)), log_file)
                    @threads for _gid in group_ids
                        _pids = df_tops_detail_mt.pheno_id[df_tops_detail_mt.group_id .== _gid]
                        group_size = df_tops_mt.group_size[findfirst(df_tops_mt.group_id .== _gid)]
                        _pval_each_mt = reshape(_cis_qtl_pairs[_cis_qtl_pairs.pheno_id .∈ (_pids,), pval_colname],:,length(_pids))
                        _pval_each_mt = _pval_each_mt[sum(isnan.(_pval_each_mt),dims=2) .== 0,:]
                        _pval_mt = vec(minimum(_pval_each_mt,dims=2))
                        pval_adj = NaN
                        pval_adj_each = repeat([NaN],group_size)
                        _pval_mt = _pval_mt[.!isnan.(_pval_mt)]
                        @runif length(_pval_mt)>0 if multiple_testing_method == "acat"
                            pval_adj = ACATest(_pval_mt; is_check=false)
                            pval_adj_each .= [ACATest(x; is_check=false) for x in eachcol(_pval_each_mt)]
                        elseif multiple_testing_method == "storey"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochbergAdaptive(Storey(storey_lambda)))) 
                            pval_adj_each .= [minimum(MultipleTesting.adjust(x, BenjaminiHochbergAdaptive(Storey(storey_lambda)))) for x in eachcol(_pval_each_mt)]
                        elseif multiple_testing_method == "bh"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochberg()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochberg())) for x in eachcol(_pval_each_mt)]
                        elseif multiple_testing_method == "bonf"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, Bonferroni()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, Bonferroni())) for x in eachcol(_pval_each_mt)]
                        elseif multiple_testing_method == "by"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiYekutieli()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, BenjaminiYekutieli())) for x in eachcol(_pval_each_mt)]
                        end
                        df_tops_mt[df_tops_mt.group_id.==_gid, pval_adj_colname] .= pval_adj
                        df_tops_detail_mt[df_tops_detail_mt.group_id.==_gid, pval_adj_colname] .= pval_adj_each
                    end
                end
            end
            pval_mt = df_tops_mt[:, pval_adj_colname]
            if any(isnan.(pval_mt))
                error(string("Cannot found ", df_tops_mt.pheno_id[isnan.(pval_mt)], " in *.cis_qtl_pairs.txt.gz files."))
            end
            println_to_file(string("    *** Note: '", qval_colname, "' will be re-calculated based on '", pval_adj_colname, "'."), log_file)
            qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochberg())
            df_tops_mt[:, qval_colname] .= qval_mt
        end
    end
    if multiple_testing_method in ["beta_approx", "clipper"]
        DOF = nothing
        if !isfile(replace(_args_cis_output[1], r"txt.gz$" => "info"))
            req_colnames = ["dof"]
            if !issubset(req_colnames, names(df_tops))
                error(req_colnames, " are required for `--mode cis_mt`")
            end
            DOF = df_tops.dof[1]
        else
            df_info = CSV.File(replace(_args_cis_output[1], r"txt.gz$" => "info")) |> DataFrame
            DOF = df_info.dof[1]
        end
        perm_file = joinpath(summary_dir, file_prefix) * ".perm.txt.gz"
        if !isfile(perm_file)
            error("Cannot found " * perm_file)
        end
        df_perm = CSV.File(perm_file) |> DataFrame
    end
    if multiple_testing_method == "beta_approx"
        df_tops_mt.pval_beta .= NaN
        df_tops_mt.beta_shape1 .= NAN
        df_tops_mt.beta_shape2 .= NAN
        df_tops_mt.true_dof .= NAN
        df_tops_mt.pval_true_dof .= NaN
        for i in axes(df_tops_mt, 1)
            r2_nominal = abs2(df_perm[i, ["X1", "X2"][n_tests]])
            if n_tests == 1
                absr_perm = df_perm[i, 4:end] |> Vector
            elseif n_tests == 2
                absr_perm = df_perm[i, ((2-1)+6):3:size(df_perm, 2)] |> Vector
            end
            df_tops_mt[i, ["pval_beta", "beta_shape1", "beta_shape2", "true_dof", "pval_true_dof"]] .= calculate_beta_approx_pval(abs2.(absr_perm), r2_nominal, DOF, 1e-4)
        end
        pval_mt = df_tops_mt.pval_beta
        qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
        if sum(qval_mt .<= _args_fdr) > 0
            set0_indices = findall(qval_mt .<= _args_fdr)
            set1_indices = findall(qval_mt .> _args_fdr)
            pthreshold = (sort(df_tops_mt.pval_beta[set1_indices])[1] - sort(-1.0 * df_tops_mt.pval_beta[set0_indices])[1]) / 2
            funbeta = Beta.(df_tops_mt.beta_shape1, df_tops_mt.beta_shape2)
            df_tops_mt[:, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]] .= quantile.(funbeta, pthreshold)
        end
        df_tops_mt[:, ["qval_g1", "qval_g2"][n_tests]] .= qval_mt
    end
    if multiple_testing_method == "clipper"
        if !isnothing(_args_pheno_group_file)
            reduced_df_tops = DataFrames.combine(groupby(df_tops, :group_id)) do grouped
                min_row_index = first_nonnan_argmin(grouped.pval_g1)
                return grouped[min_row_index, :]
            end
            reduced_df_tops.group_size .= DataFrames.combine(groupby(df_tops, :group_id), nrow).nrow
            reduced_df_perm = copy(df_perm)
            reduced_df_perm.pheno_id .= df_tops.group_id
            df_perm_gdf = groupby(reduced_df_perm, :pheno_id)
            reduced_df_perm = DataFrames.combine(df_perm_gdf, :pheno_id, valuecols(df_perm_gdf) .=> maximum)
            reduced_df_perm = unique(reduced_df_perm, [:pheno_id])
            df_perm_clipper = reduced_df_perm
            df_tops_mt = reduced_df_tops
        else
            df_perm_clipper = df_perm
            df_tops_mt = df_tops
        end
        clipper_keep_index = .!isnan.(df_perm_clipper[:, 3])
        if n_tests == 1
            df_tops_mt.qval_g1 .= NaN
            if n_perms < 1000
                res_clipper = Clipper(Matrix(df_perm_clipper[clipper_keep_index, 3:3]), Matrix(df_perm_clipper[clipper_keep_index, 4:end]), analysis="enrichment", procedure="GZ", contrast_score="max", FDR=[FloatT(_args_fdr)])
                qval_mt = res_clipper["q"]
            else
                numsOfPermsMoreSig = sum(Matrix(df_perm_clipper[clipper_keep_index, 3:3]) .- Matrix(df_perm_clipper[clipper_keep_index, 4:end]) .<= 0, dims=2)
                pval_mt = vec((numsOfPermsMoreSig .+ 1) ./ (n_perms + 1))
                qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
            end
            df_tops_mt.qval_g1[clipper_keep_index] .= qval_mt
            n_significant_phenotypes = sum((df_tops_mt.qval_g1 .< FloatT(_args_fdr)) .& (df_tops_mt.pval_g1 .< df_tops_mt.pval_g1_threshold))
            logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(clipper_keep_index), " tested phenotypes with significant cis-QTL.")
            println_to_file(logtxt, log_file)
        elseif n_tests == 2
            for i in 2:2
                idx_exp = ((i-1)+3):((i-1)+3)
                idx_back = ((i-1)+6):3:size(df_perm_clipper, 2) |> Vector
                n_perms = length(idx_back)
                if n_perms < 1000
                    res_clipper = Clipper(Matrix(df_perm_clipper[clipper_keep_index, idx_exp]), Matrix(df_perm_clipper[clipper_keep_index, idx_back]), analysis="enrichment", procedure="GZ", contrast_score="max", FDR=[FloatT(_args_fdr)])
                    qval_mt = res_clipper["q"]
                else
                    numsOfPermsMoreSig = sum(Matrix(df_perm_clipper[clipper_keep_index, idx_exp]) .- Matrix(df_perm_clipper[clipper_keep_index, idx_back]) .<= 0, dims=2)
                    pval_mt = vec((numsOfPermsMoreSig .+ 1) ./ (n_perms + 1))
                    qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochbergAdaptive(Storey(_args_storey_lambda)))
                end
                if i == 1
                    df_tops_mt.qval_g1 .= NaN
                    df_tops_mt.qval_g1[clipper_keep_index] .= qval_mt
                    logtxt = string("  *** ", sum(qval_mt .< FloatT(_args_fdr)), " out of ", sum(clipper_keep_index), " tested phenotypes with significant cis-QTL(g1).")
                    println_to_file(logtxt, log_file)
                elseif i == 2
                    pval_g2_threshold = [sort(get_approx_p_from_r(Vector(df_perm_clipper[pix, idx_back]), DOF - 1))[ceil(Int, n_perms * _args_fdr)] for pix in axes(df_perm_clipper, 1)] 
                    df_tops_mt.pval_g2_threshold .= pval_g2_threshold
                    df_tops_mt.qval_g2 .= NaN
                    df_tops_mt.qval_g2[clipper_keep_index] .= qval_mt
                    n_significant_phenotypes = sum((df_tops_mt.qval_g2 .< FloatT(_args_fdr)) .& (df_tops_mt.pval_g2 .< df_tops_mt.pval_g2_threshold))
                    logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(clipper_keep_index), " tested phenotypes with significant cis-QTL(g2).")
                    println_to_file(logtxt, log_file)
                elseif i == 3
                    df_tops_mt.qval_joint .= NaN
                    df_tops_mt.qval_joint[clipper_keep_index] .= qval_mt
                end
            end
        end
    end
    cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
    df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt.txt.gz")), df_tops_mt, delim="\t", compress=true)
    if with_group
        cols_float = [df_tops_detail_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_detail_mt))]
        df_tops_detail_mt[:, cols_float] = round.(df_tops_detail_mt[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt_detail.txt.gz")), df_tops_detail_mt, delim="\t", compress=true)
    end
    println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
    if _args_debug
        println_to_file(string(to), log_file)
        println()
    end
    return 0
end
