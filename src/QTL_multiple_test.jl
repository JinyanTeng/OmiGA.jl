function runOmiGA_cis_mt()
    time_start = now()
    df_tops = DataFrame()
    df_tops_detail = nothing
    if length(_args_cis_output) > 0
        file_cis_qtl_txt = _args_cis_output
        df_tops = CSV.File(file_cis_qtl_txt, header=true, buffer_in_memory=true) |> DataFrame
        if issubset(["group_id"], names(df_tops))
            if length(file_cis_qtl_txt) == 1
                file_cis_qtl_detail_txt = replace(file_cis_qtl_txt[1], ".txt.gz" => ".detail.txt.gz")
            else
                match_suffixes = [match(r".[0-9]+.txt|.[0-9]+.txt.gz", fn).match for fn in file_cis_qtl_txt]
                match_prefixes = replace.(file_cis_qtl_txt, r"[0-9]+.txt|[0-9]+.txt.gz" => "")
                file_cis_qtl_detail_txt = string.(match_prefixes, "detail", match_suffixes)
            end
            df_tops_detail = CSV.File(file_cis_qtl_detail_txt, header=true, buffer_in_memory=true) |> DataFrame
        end
    else
        error("Please provide correct file for `--cis-file`!")
    end
    with_group = !isnothing(df_tops_detail)
    n_tests = any(.!isnothing.(match.(r"g2", names(df_tops)))) ? 2 : 1
    if endswith(_args_cis_output[1], ".cis_mt.txt.gz")
        cis_summary_prefix = replace(_args_cis_output[1], r".[a-z,A-Z,0-9]+.cis_mt.txt.gz" => "")
    elseif endswith(_args_cis_output[1], r".cis_qtl.txt.gz|.cis_qtl.[0-9]+.txt|.cis_qtl.[0-9]+.txt.gz")
        cis_summary_prefix = replace(_args_cis_output[1], r".cis_qtl.txt.gz|.cis_qtl.[0-9]+.txt|.cis_qtl.[0-9]+.txt.gz" => "")
    else
        error("Incorrect summary file was provided.")
    end
    summary_dir = dirname(cis_summary_prefix)
    file_prefix = basename(cis_summary_prefix)
    if _args_output_significant_qtls
        use_detail_thresholds = false
        if !with_group
            significant_phenotypes = df_tops.pheno_id[df_tops[:, ["qval_g1", "qval_g2"][n_tests]].<_args_fdr]
        else
            significant_groups = df_tops.group_id[df_tops[:, ["qval_g1", "qval_g2"][n_tests]].<_args_fdr]
            significant_phenotypes = df_tops_detail.pheno_id[df_tops_detail.group_id.∈(significant_groups,)]
            use_detail_thresholds = any(.!isnothing.(match.(r"threshold", names(df_tops_detail))))
        end
        full_summary_files = _args_cis_pairs_output 
        if length(full_summary_files) == 0
            all_files = readdir(summary_dir)
            full_summary_files = all_files[.!isnothing.(findfirst.(string(file_prefix, ".cis_qtl_pairs"), all_files))]
        else
            summary_dir = dirname(cis_summary_prefix)
        end
        println_to_file(string("    ** Found ", length(full_summary_files), " *.cis_qtl_pairs.* files"), log_file)
        if length(full_summary_files) == 0
            error("All required *.cis_qtl_pairs.* files must be located in the same directory as '", basename(_args_cis_output[1]), "'. Alternatively, you can provide their paths using the '--cis-qtl-pairs' option.")
        else
            if length(unique(replace.(full_summary_files, r".txt.gz|.arrow|.jld|.jld2" => ""))) != length(full_summary_files)
                error("Found duplicated *.cis_qtl_pairs.* files with different suffixes, please specify the '--cis-qtl-pairs' option to prevent error.")
            end
            sorted_indices = sortperm(replace.(full_summary_files, r".txt.gz|.arrow|.jld|.jld2" => "", r"(.)*.cis_qtl_pairs." => ""), lt=compare_chromosomes) 
            full_summary_files = full_summary_files[sorted_indices]
            n_significant_qtls = 0
            out_file = joinpath(_args_output_dir, string(_args_out_prefix, ".significant_qtl_pairs.txt.gz")) 
            for fn in full_summary_files
                _df_sig_qtls_df = DataFrame()
                println_to_file(string("    * Processing ", fn), log_file)
                _cis_qtl_pairs = read_qtl_pairs_file(joinpath(summary_dir, fn))
                if use_detail_thresholds
                    _cis_qtl_pairs = _cis_qtl_pairs[_cis_qtl_pairs[:, ["pval_g1", "pval_g2"][n_tests]].<=nanmaximum(df_tops_detail[:, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]), :]
                else
                    _cis_qtl_pairs = _cis_qtl_pairs[_cis_qtl_pairs[:, ["pval_g1", "pval_g2"][n_tests]].<=nanmaximum(df_tops[:, ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]), :]
                end
                pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                pheno_ids = pheno_ids[findall(pheno_ids .∈ (significant_phenotypes,))]
                for _pid in pheno_ids
                    _pval_mt = _cis_qtl_pairs[_cis_qtl_pairs.pheno_id.==_pid, :]
                    if use_detail_thresholds
                        pval_nominal_threshold = df_tops_detail[findfirst(df_tops_detail.pheno_id .== _pid), ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]
                        if !isnan(pval_nominal_threshold)
                            append!(_df_sig_qtls_df, _pval_mt[_pval_mt[:, ["pval_g1", "pval_g2"][n_tests]].<=pval_nominal_threshold, :])
                        end
                    else
                        if with_group
                            _gid = df_tops_detail.group_id[findfirst(df_tops_detail.pheno_id .== _pid)]
                            pval_nominal_threshold = df_tops[findfirst(df_tops.group_id .== _gid), ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]
                            if !isnan(pval_nominal_threshold)
                                append!(_df_sig_qtls_df, _pval_mt[_pval_mt[:, ["pval_g1", "pval_g2"][n_tests]].<=pval_nominal_threshold, :])
                            end
                        else
                            pval_nominal_threshold = df_tops[findfirst(df_tops.pheno_id .== _pid), ["pval_g1_threshold", "pval_g2_threshold"][n_tests]]
                            if !isnan(pval_nominal_threshold)
                                append!(_df_sig_qtls_df, _pval_mt[_pval_mt[:, ["pval_g1", "pval_g2"][n_tests]].<=pval_nominal_threshold, :])
                            end
                        end
                    end
                end
                n_significant_qtls += size(_df_sig_qtls_df, 1)
                CSV.write(out_file, _df_sig_qtls_df, delim="\t", compress=true, append=fn != full_summary_files[1])
            end
            println_to_file(string(" *** Extracted a total of ", n_significant_qtls, " significant variant-phenotype pairs."), log_file)
        end
        return 0
    end
    if isnothing(_args_multiple_testing)
        error("--multiple-testing must be specified!")
    else
        multiple_testing_method = _args_multiple_testing
    end
    calcu_variant_threshold = true
    if multiple_testing_method == "acat2"
        multiple_testing_method = "acat"
        use_acat2 = true
    else
        use_acat2 = false
    end
    if multiple_testing_method == "acat"
        use_acat2 = true
    end
    storey_lambda = _args_storey_lambda
    out_file_middle_part = multiple_testing_method
    df_tops_mt = copy(df_tops)
    df_tops_mt_detail = !with_group ? nothing : df_tops_detail
    chroms = nothing
    if issubset(["chrom"], names(df_tops_mt))
        chroms = unique(df_tops_mt.chrom)
    end
    mt_method_list = Dict("acat" => "ACAT", "storey" => "Storey", "bh" => "BenjaminiHochberg", "bonf" => "Bonferroni", "by" => "BenjaminiYekutieli", "clipper" => "ClipperQTL", "beta_approx" => "Beta approximation")
    use_orig_mt_method = false
    if multiple_testing_method in ["beta_approx", "clipper"]
        DOF = nothing
        n_perms = nothing
        if !isfile(replace(_args_cis_output[1], r"txt.gz$" => "info"))
            req_colnames = ["dof"]
            if !issubset(req_colnames, names(df_tops))
                error(req_colnames, " are required for `--mode cis_mt`")
            end
            DOF = df_tops.dof[1]
        else
            df_info = CSV.File(replace(_args_cis_output[1], r"txt.gz$" => "info")) |> DataFrame
            DOF = df_info.dof[1]
            n_perms = df_info.num_perm[1]
            orig_mt_method = df_info.mt_method[1]
            use_orig_mt_method = multiple_testing_method == orig_mt_method
        end
        if length(_args_cis_output) == 1
            file_perm_txt = replace(_args_cis_output[1], ".txt.gz" => ".perm.arrow")
        else
            match_suffixes = [match(r".[0-9]+.txt|.[0-9]+.txt.gz", fn).match for fn in _args_cis_output]
            match_prefixes = replace.(_args_cis_output, r"[0-9]+.txt|[0-9]+.txt.gz" => "")
            file_perm_txt = string.(match_prefixes, "perm", replace.(match_suffixes, r".txt$|.txt.gz$" => ".arrow"))
        end
        df_perm = Arrow.Table(file_perm_txt) |> DataFrame
    elseif multiple_testing_method in ["acat"]
        df_info = CSV.File(replace(_args_cis_output[1], r"txt.gz$|[0-9]+.txt$|[0-9]+.txt.gz$" => "info")) |> DataFrame
        orig_mt_method = df_info.mt_method[1]
        use_orig_mt_method = multiple_testing_method == orig_mt_method
    end
    if multiple_testing_method in ["acat", "storey", "bh", "bonf", "by"]
        println_to_file(string("    * Re-perform multiple-testing correction: first-level correction -> ", mt_method_list[multiple_testing_method], ", second-level correction -> ", mt_method_list["bh"]), log_file)
        @runif use_acat2 println_to_file(string("    * Use nested ACAT test for each phenotype group."), log_file)
        @runif calcu_variant_threshold != "acat" println_to_file(string("    * 'pval_g", n_tests, "_threshold' will not be updated for ", multiple_testing_method, "."), log_file)
        if multiple_testing_method == "storey"
            println_to_file(string("    * lambda=", storey_lambda, " was used for Storey."), log_file)
            out_file_middle_part = string(multiple_testing_method, "_", storey_lambda)
        end
        pval_colname = ["pval_g1", "pval_g2"][n_tests]
        rm_colnames_list = pval_colname * "_" .* ["acat", "storey", "bh", "bonf", "by"]
        if !with_group
            pval_adj_colname = ["pval_g1_pheno", "pval_g2_group"][n_tests]
        else
            pval_adj_colname = ["pval_g1_group", "pval_g2_group"][n_tests]
        end
        qval_colname = ["qval_g1", "qval_g2"][n_tests]
        rm_colnames = rm_colnames_list[rm_colnames_list.∈(names(df_tops_mt),)]
        @runif length(rm_colnames) > 0 select!(df_tops_mt, Not(rm_colnames))
        @runif issubset([qval_colname], names(df_tops_mt)) select!(df_tops_mt, Not(qval_colname))
        @runif !use_orig_mt_method df_tops_mt[:, pval_adj_colname] .= NaN
        @runif !use_orig_mt_method df_tops_mt[:, qval_colname] .= NaN
        if with_group
            rm_colnames = rm_colnames_list[rm_colnames_list.∈(names(df_tops_mt_detail),)]
            @runif length(rm_colnames) > 0 select!(df_tops_mt_detail, Not(rm_colnames))
            @runif !use_orig_mt_method df_tops_mt_detail[:, pval_adj_colname] .= NaN
        end
        full_summary_files = _args_cis_pairs_output 
        if length(full_summary_files) == 0
            all_files = readdir(summary_dir)
            full_summary_files = all_files[.!isnothing.(findfirst.(string(file_prefix, ".cis_qtl_pairs"), all_files))]
        else
            summary_dir = dirname(cis_summary_prefix)
        end
        println_to_file(string("    ** Found ", length(full_summary_files), " *.cis_qtl_pairs.* files"), log_file)
        if length(full_summary_files) == 0
            error("All required *.cis_qtl_pairs.*.txt.gz files must be located in the same directory as '", basename(_args_cis_output[1]), "'. Alternatively, you can provide their paths using the '--cis-qtl-pairs' option.")
        else
            if length(unique(replace.(full_summary_files, r".txt.gz|.arrow|.jld|.jld2" => ""))) != length(full_summary_files)
                error("Found duplicated *.cis_qtl_pairs.* files with different suffixes, please specify the '--cis-qtl-pairs' option to prevent error.")
            end
            sorted_indices = sortperm(replace.(full_summary_files, r".txt.gz|.arrow|.jld|.jld2" => "", r"(.)*.cis_qtl_pairs." => ""), lt=compare_chromosomes) 
            full_summary_files = full_summary_files[sorted_indices]
            if isnothing(chroms)
                chroms = replace.(full_summary_files, r".txt.gz|.arrow|.jld|.jld2" => "", r"(.)*.cis_qtl_pairs." => "")
            end
            df_condidate_pval_threshold = DataFrame()
            @runif !use_orig_mt_method for fn in full_summary_files
                println_to_file(string("    *** Processing ", fn), log_file)
                _cis_qtl_pairs = read_qtl_pairs_file(joinpath(summary_dir, fn))
                pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                group_ids = nothing
                @runif with_group group_ids = unique(df_tops_mt_detail.group_id[df_tops_mt_detail.pheno_id.∈(pheno_ids,)])
                println_to_file(string("        Number of phenotypes: ", length(pheno_ids)), log_file)
                if !with_group
                    spinlock = Threads.SpinLock()
                    @threads for _pid in pheno_ids
                        _pval_mt = _cis_qtl_pairs[_cis_qtl_pairs.pheno_id.==_pid, pval_colname]
                        pval_adj = NaN
                        if multiple_testing_method == "acat"
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
                    spinlock = Threads.SpinLock()
                    @threads for _gid in group_ids
                        _pids = df_tops_mt_detail.pheno_id[df_tops_mt_detail.group_id.==_gid]
                        group_size = df_tops_mt.group_size[findfirst(df_tops_mt.group_id .== _gid)]
                        _pval_group = reshape(_cis_qtl_pairs[_cis_qtl_pairs.pheno_id.∈(_pids,), pval_colname], :, length(_pids))
                        _pval_mt = vec(minimum(_pval_group, dims=2))
                        pval_adj = NaN
                        pval_adj_each = repeat([NaN], group_size)
                        if multiple_testing_method == "acat"
                            pval_adj_each .= [ACATest(x; is_check=false) for x in eachcol(_pval_group)]
                            if !use_acat2
                                pval_adj = ACATest(_pval_mt; is_check=false)
                            elseif use_acat2
                                pval_adj = ACATest(pval_adj_each; is_check=false)
                            end
                        elseif multiple_testing_method == "storey"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochbergAdaptive(Storey(storey_lambda)))) 
                            pval_adj_each .= [minimum(MultipleTesting.adjust(x, BenjaminiHochbergAdaptive(Storey(storey_lambda)))) for x in eachcol(_pval_group)]
                        elseif multiple_testing_method == "bh"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochberg()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, BenjaminiHochberg())) for x in eachcol(_pval_group)]
                        elseif multiple_testing_method == "bonf"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, Bonferroni()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, Bonferroni())) for x in eachcol(_pval_group)]
                        elseif multiple_testing_method == "by"
                            pval_adj = minimum(MultipleTesting.adjust(_pval_mt, BenjaminiYekutieli()))
                            pval_adj_each .= [minimum(MultipleTesting.adjust(_pval_mt, BenjaminiYekutieli())) for x in eachcol(_pval_group)]
                        end
                        df_tops_mt[df_tops_mt.group_id.==_gid, pval_adj_colname] .= pval_adj
                        df_tops_mt_detail[df_tops_mt_detail.group_id.==_gid, pval_adj_colname] .= pval_adj_each
                    end
                end
                if _args_debug
                    cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
                    df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
                    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt.txt.gz")), df_tops_mt, delim="\t", compress=true)
                    if with_group
                        cols_float = [df_tops_mt_detail[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt_detail))]
                        df_tops_mt_detail[:, cols_float] = round.(df_tops_mt_detail[:, cols_float], sigdigits=6)
                        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt.detail.txt.gz")), df_tops_mt_detail, delim="\t", compress=true)
                    end
                end
                println_to_file(string("        Elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
            end
            pval_mt = df_tops_mt[:, pval_adj_colname]
            if any(isnan.(pval_mt))
                @warn string("Omit ", df_tops_mt.pheno_id[isnan.(pval_mt)], " due to missing P-values.")
                df_tops_mt = df_tops_mt[.~isnan.(pval_mt),:]
            end
            println_to_file(string("    *** Note: '", qval_colname, "' were re-calculated based on '", pval_adj_colname, "'."), log_file)
            qval_mt = MultipleTesting.adjust(pval_mt, BenjaminiHochberg())
            df_tops_mt[:, qval_colname] .= qval_mt
            mt_col = n_tests > 2 ? "joint" : string("g", n_tests)
            @timeit to "Get final acat threshold" @runif calcu_variant_threshold && multiple_testing_method == "acat" begin
                println_to_file(string("    Calculating the ", string("pval_", mt_col, "_threshold"), " ..."), log_file)
                fdr_closest_i = find_closest_below_fast(qval_mt, _args_fdr)
                @timeit to "Get 1" if fdr_closest_i != -1
                    fdr_closest_i_up = find_closest_up_fast(qval_mt, _args_fdr)
                    if with_group
                        critical_acat_p = sqrt(df_tops_mt[fdr_closest_i, string("pval_", mt_col, "_group")] * df_tops_mt[fdr_closest_i_up, string("pval_", mt_col, "_group")])
                    else
                        critical_acat_p = sqrt(df_tops_mt[fdr_closest_i, string("pval_", mt_col, "_pheno")] * df_tops_mt[fdr_closest_i_up, string("pval_", mt_col, "_pheno")])
                    end
                    for fn in full_summary_files
                        println_to_file(string("    *** Processing ", fn), log_file)
                        _cis_qtl_pairs = read_qtl_pairs_file(joinpath(summary_dir, fn))
                        pheno_ids = unique(_cis_qtl_pairs.pheno_id)
                        group_ids = nothing
                        println_to_file(string("        Number of phenotypes: ", length(pheno_ids)), log_file)
                        if with_group 
                            group_ids = unique(df_tops_mt_detail.group_id[df_tops_mt_detail.pheno_id.∈(pheno_ids,)])
                            println_to_file(string("        Number of groups: ", length(group_ids)), log_file)
                        end
                        @threads for _pid in pheno_ids
                            _pval_mt = _cis_qtl_pairs[_cis_qtl_pairs.pheno_id.==_pid, pval_colname]
                            _pval_mt = _pval_mt[.!isnan.(_pval_mt)]
                            if length(_pval_mt) > 0
                                _, pval_nominal_threshold, _ = calcu_target_pval_threshold_acat!(_pval_mt, critical_acat_p)
                                if !with_group 
                                    df_tops_mt[findfirst(df_tops_mt.pheno_id.==_pid), string("pval_", mt_col, "_threshold")] = pval_nominal_threshold
                                else
                                    df_tops_mt_detail[findfirst(df_tops_mt_detail.pheno_id.==_pid), string("pval_", mt_col, "_threshold")] = pval_nominal_threshold
                                end
                            end
                        end
                        if with_group
                            matchidx = vmatch(df_tops_mt_detail.pheno_id, df_tops_mt.pheno_id)
                            df_tops_mt[.!isnothing.(matchidx), string("pval_", mt_col, "_threshold")] .= df_tops_mt_detail[matchidx[.!isnothing.(matchidx)], string("pval_", mt_col, "_threshold")]
                            cols_float = [df_tops_mt_detail[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt_detail))]
                            df_tops_mt_detail[:, cols_float] = round.(df_tops_mt_detail[:, cols_float], sigdigits=6)
                            CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_mt.detail.txt.gz")), df_tops_mt_detail, delim="\t", compress=true)
                        end
                    end
                end
                println_to_file(string("    Done."), log_file)
            end
            cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
            df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
            CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".cis_mt.txt.gz")), df_tops_mt, delim="\t", compress=true)
            n_significant_phenotypes = sum(qval_mt .< FloatT(_args_fdr))
            logtxt = string("  *** ", n_significant_phenotypes, " out of ", length(qval_mt), " tested phenotypes with qval_", mt_col, "<", _args_fdr, ".")
            println_to_file(logtxt, log_file)
        end
        return 0
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
        println_to_file(string("    * Re-perform multiple-testing correction using ", mt_method_list[multiple_testing_method]), log_file)
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
                if calcu_variant_threshold
                    pval_nominal_threshold = get_approx_p_from_r(sort(Matrix(df_perm_mt[:, 4:end]); dims=2, rev=true)[:, ceil(Int, n_perms * _args_fdr)], DOF)
                    df_tops_mt.pval_g1_threshold .= pval_nominal_threshold
                end
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
            logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(mt_keep_index), " tested phenotypes with qval_g", n_tests, "<", _args_fdr, ".")
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
                    if calcu_variant_threshold
                        pval_nominal_threshold = get_approx_p_from_r(sort(Matrix(df_perm_mt[:, idx_back]); dims=2, rev=true)[:, ceil(Int, n_perms * _args_fdr)], DOF - 1)
                        df_tops_mt.pval_g2_threshold .= pval_nominal_threshold
                    end
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
                logtxt = string("  *** ", n_significant_phenotypes, " out of ", sum(mt_keep_index), " tested phenotypes with qval_", ["g1", "g2", "g_joint"][i], "<", _args_fdr, ".")
                println_to_file(logtxt, log_file)
            end
        end
    end
    cols_float = [df_tops_mt[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt))]
    df_tops_mt[:, cols_float] = round.(df_tops_mt[:, cols_float], sigdigits=6)
    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt.txt.gz")), df_tops_mt, delim="\t", compress=true)
    if with_group
        cols_float = [df_tops_mt_detail[1, x] isa AbstractFloat for x in range(1, ncol(df_tops_mt_detail))]
        df_tops_mt_detail[:, cols_float] = round.(df_tops_mt_detail[:, cols_float], sigdigits=6)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".", out_file_middle_part, ".cis_mt.detail.txt.gz")), df_tops_mt_detail, delim="\t", compress=true)
    end
    println_to_file(string("+++ Total elapsed time (h:m:s:ms): ", format_milliseconds(now() - time_start)), log_file)
    if _args_debug
        println_to_file(string(to), log_file)
        println()
    end
    return 0
end
