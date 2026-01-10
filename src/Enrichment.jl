function getIntervalCollection(annotbed::DataFrame)
        ncol = size(annotbed,2)
        annotbed = rename(annotbed, [:Chromosome, :Start, :End, :Name, Symbol.(string.("Score", 1:(ncol-4)))...])
        annotbed = Tables.rowtable(annotbed)
        result = GenomicIntervalCollection([GenomicInterval(qy.Chromosome, qy.Start, qy.End, '?', qy.Name) for qy in annotbed], true)
        return result
end
function getIntervalCollection(fs::String)
    annotbed = CSV.File(fs, delim="\t", header=0) |> DataFrame
    getIntervalCollection(annotbed)
end
function getInterval(A::GenomicIntervalCollection, B::GenomicIntervalCollection; only_meta::Union{Nothing,Vector{String}}=nothing)
    j = 0
    for i in eachoverlap(A, B)
        j+=1
    end
    if isnothing(only_meta)
        df = DataFrame(groupname_A=repeat([""],j),first_A=1,last_A=1,metadata_A="",groupname_B="",first_B=1,last_B=1,metadata_B="")
        j = 0
        for i in eachoverlap(A, B)
            j+=1
            df[j,:] .= (i[1].groupname, i[1].first, i[1].last, i[1].metadata, i[2].groupname, i[2].first, i[2].last, i[2].metadata)
        end
    else
        df = DataFrame()
        if "A" in only_meta
            df.metadata_A = repeat([""],j)
        end
        if "B" in only_meta
            df.metadata_B = repeat([""],j)
        end
        j = 0
        for i in eachoverlap(A, B)
            j+=1
            if "A" in only_meta
                df.metadata_A[j] = i[1].metadata
            end
            if "B" in only_meta
                df.metadata_B[j] = i[2].metadata
            end
        end
    end
    return df
end
function read_cis_file(cis_file::String; threshold::Union{Nothing,Float64}=nothing)
    df_tops = CSV.read(cis_file, DataFrame, header=true, select=[:variant_id], buffer_in_memory=true)
    if !isnothing(threshold)
        subset!(df_tops, :qval_g1 => x -> x .<= threshold)
    end
    return df_tops
end
function read_bkg_bim(bkg_plink_prefix::String; calcu_maf::Bool=true)
    bkg_bim = CSV.read(string(bkg_plink_prefix, ".bim"), DataFrame, header=false, select=[:1, :2, :4], types=[String, String, Int, Int, String1, String1], buffer_in_memory=!_args_low_mem)
    rename!(bkg_bim, :1 => :chr, :2 => :variant_id, :3 => :pos)
    if calcu_maf
        bkg_bim.maf = get_allele_maf(bkg_plink_prefix)
        bed = DataFrame("chr" => replace.(bkg_bim.chr, r"^chr" => ""), "start" => bkg_bim.pos, "end" => bkg_bim.pos, "variant_id" => bkg_bim.variant_id, "maf" => bkg_bim.maf) 
    else
        bed = DataFrame("chr" => replace.(bkg_bim.chr, r"^chr" => ""), "start" => bkg_bim.pos, "end" => bkg_bim.pos, "variant_id" => bkg_bim.variant_id) 
    end
    return bed
end
function add_chr_pos_to_cis_tops(target_data::DataFrame, bkg_data::DataFrame)
    cis_tops_bed = leftjoin(target_data, bkg_data, on = :variant_id)
    cis_tops_bed.chr = cis_tops_bed.chr |> Vector{String}
    cis_tops_bed.start = cis_tops_bed.start |> Vector{Int}
    cis_tops_bed.end = cis_tops_bed.end |> Vector{Int}
    cis_tops_bed.maf = cis_tops_bed.maf |> Vector{Float64}
    if size(bkg_data,2) == 5
        return cis_tops_bed[:,["chr","start","end","variant_id","maf"]]
    elseif size(bkg_data,2) == 6
        cis_tops_bed.ldscore = cis_tops_bed.ldscore |> Vector{Float64}
        return cis_tops_bed[:,["chr","start","end","variant_id","maf","ldscore"]]
    end
end
function get_bkg_maf_match!(bkg_perm_df::DataFrame, unique_mafs::Vector{Float64}, freq_mafs::Vector{Int}, bkg_data_sorted::DataFrame, maf_match::Float64)
    n_mafs = length(unique_mafs)
    n_target = sum(freq_mafs)
    selected_indices = zeros(Int, n_target)
    j1 = 1
    j2 = 0
    for i in 1:n_mafs
        low = unique_mafs[i] - maf_match
        high = unique_mafs[i] + maf_match
        range_low = searchsortedfirst(bkg_data_sorted.maf, low)
        range_high = searchsortedlast(bkg_data_sorted.maf, high)
        j2 += freq_mafs[i]
        selected_indices[j1:j2] .= sample(range_low:range_high,freq_mafs[i]; replace=true, ordered=false)
        j1 = j2 + 1
    end
    bkg_perm_df .= bkg_data_sorted[sort(selected_indices), 1:4]
end
function get_bkg_maf_ld_match!(bkg_perm_df::DataFrame, unique_mafs::Vector{Float64}, freq_mafs::Vector{Int}, bkg_data_sorted::DataFrame, maf_match::Float64; target_data::Union{Nothing,DataFrame}=nothing,ld_match::Union{Nothing,Float64}=nothing)
    n_mafs = length(unique_mafs)
    n_target = sum(freq_mafs)
    selected_indices = zeros(Int, n_target)
    j1 = 1
    j2 = 0
    @runif !isnothing(ld_match) target_ldscore_std = std(target_data.ldscore) * ld_match
    for i in 1:n_mafs
        low = unique_mafs[i] - maf_match
        high = unique_mafs[i] + maf_match
        range_low = searchsortedfirst(bkg_data_sorted.maf, low)
        range_high = searchsortedlast(bkg_data_sorted.maf, high)
        j2 += freq_mafs[i]
        if isnothing(ld_match)
            selected_indices[j1:j2] .= sample(range_low:range_high, freq_mafs[i]; replace=true, ordered=false)
        else
            target_ldscore = sort(target_data.ldscore[target_data.maf .== unique_mafs[i]])
            ld_match_indices = zeros(freq_mafs[i])
            range_low_high = range_low:range_high
            bkg_data_sorted_ldscore_subset = bkg_data_sorted.ldscore[range_low_high]
            for lds_i in eachindex(target_ldscore)
                lds = target_ldscore[lds_i]
                lds_low = lds - target_ldscore_std
                lds_high = lds + target_ldscore_std
                range_subset_ld_match = range_low_high[lds_low .<= bkg_data_sorted_ldscore_subset .<= lds_high]
                ld_match_indices[lds_i] = sample(range_subset_ld_match) 
            end
            selected_indices[j1:j2] .= ld_match_indices
        end
        j1 = j2 + 1
    end
    bkg_perm_df .= bkg_data_sorted[sort(selected_indices), 1:4]
end
function read_annot_bed(annot_file_list::String)
    chromatin_category_data = CSV.read(annot_file_list, DataFrame, header = false)
    rename!(chromatin_category_data, 1 => :chr, 2 => :start, 3 => :end, 4 => :category)
    chromatin_category_data.chr .= replace.(chromatin_category_data.chr, r"^chr"=>"")
    chromatin_category_data.start .+= 1 
    return chromatin_category_data
end
function enrichment_permutation(
    n_perms::Int64, 
    bkg_data::DataFrame, 
    target_data::DataFrame, 
    annotation_IC::Vector{GenomicIntervalCollection},
    annot_file_list::Union{Tuple{String}, Vector{String}};
    maf_match::Union{Float64, Nothing}=nothing,
    ld_match::Union{Float64, Nothing}=nothing,
    use_region::Bool=false
    )
    if use_region
        maf_match = nothing
    else
        if !isnothing(maf_match)
            target_data_maf_match = sort(target_data, :maf)
            if isnothing(ld_match)
                bkg_data_sorted = sort(bkg_data, :maf)
            else
                bkg_data_sorted = sort(bkg_data, [:maf, :ldscore])
            end
            unique_mafs = unique(target_data_maf_match.maf)
            freq_mafs = [count(x -> x == val, target_data_maf_match.maf) for val in unique_mafs]
        end
        if !isnothing(ld_match) && isnothing(maf_match)
            maf_match = 0.02
            println_to_file(" * The '--ld-match' option was used without specifying '--maf-match'. The '--maf-match' parameter has been set to 0.02.", log_file)
        end
    end
    n_rows_target = size(target_data, 1)
    random_result_permutation = Vector{DataFrame}(undef, n_perms)
    thread_local_storage = [(DataFrame([Vector{eltype(col)}(undef, n_rows_target) for col in eachcol(bkg_data[:,1:4])], names(bkg_data[:,1:4]))) for _ in 1:nthreads()]
    block_size = ceil(Int, n_perms / nthreads())
    iter_collects = collect(Iterators.partition(1:n_perms, block_size))
    @threads for blocki in eachindex(iter_collects) 
        tid = blocki
        bkg_perm_df = thread_local_storage[tid]
        for ra in iter_collects[blocki]
            if !isnothing(maf_match) 
                if isnothing(ld_match)
                    get_bkg_maf_match!(bkg_perm_df, unique_mafs, freq_mafs, bkg_data_sorted, maf_match)
                else
                    get_bkg_maf_ld_match!(bkg_perm_df, unique_mafs, freq_mafs, bkg_data_sorted, maf_match; target_data=target_data, ld_match=ld_match)
                end
            else
                rng = Random.MersenneTwister(ra) 
                random_idx = rand(rng, 1:size(bkg_data, 1), n_rows_target)
                bkg_perm_df .= bkg_data[random_idx, 1:4]
                if use_region
                    bkg_perm_df.end .= bkg_perm_df.end .+ target_data.len .- 1
                end
            end
            bkg_IC = getIntervalCollection(bkg_perm_df)
            bkg_overlap_results = DataFrame()
            for i in eachindex(annot_file_list)
                overlap_res = getInterval(bkg_IC, annotation_IC[i]; only_meta=["B"])
                overlap_res.bed_file .= basename(annot_file_list[i])
                append!(bkg_overlap_results, overlap_res)
            end
            bkg_overlap_results = DataFrames.combine(groupby(bkg_overlap_results, [:metadata_B, :bed_file]), nrow => :count)
            bkg_overlap_results.prop = bkg_overlap_results.count ./ n_rows_target
            random_result_permutation[ra] = bkg_overlap_results
            println_to_file("    Permutated $(ra)/$(n_perms).", log_file)
        end
    end
    return random_result_permutation
end
function calcu_enrichment_p(target_overlap_result::DataFrame, permuted_results::Vector{DataFrame}, n_perms::Int64)
    Fold_Enrichment_DF = DataFrame()
    for ra in 1:n_perms
        merged_ = innerjoin(target_overlap_result, permuted_results[ra], on=[:metadata_B, :bed_file], makeunique=true)
        merged_.Fold_Enrichment .= merged_.prop ./ merged_.prop_1 
        merged_.Random_idx .= ra
        Fold_Enrichment_DF = vcat(Fold_Enrichment_DF, merged_)
    end
    mean_sd_merged_ = DataFrames.combine(
        groupby(Fold_Enrichment_DF, [:metadata_B, :bed_file]),
        :Fold_Enrichment => (x -> Statistics.mean(skipmissing(x))) => :fold,
        :Fold_Enrichment => (x -> std(skipmissing(x))) => :sd
    )
    final_merged_ = innerjoin(target_overlap_result, mean_sd_merged_, on=[:metadata_B, :bed_file])
    Z_score = (final_merged_.fold .- 1) ./ final_merged_.sd
    final_merged_.pval = ccdf.(Ref(Normal()), Z_score)
    sort!(final_merged_, [order(:bed_file), order(:fold, rev=true)])
    select!(final_merged_, [:bed_file, :metadata_B, :count, :prop, :fold, :sd, :pval])
    rename!(final_merged_, Dict(:bed_file => :annot, :metadata_B => :category))
    return final_merged_
end
function runOmiGA_enrich(
    target_file::String,
    annot_file_list::Union{Tuple{String}, Vector{String}};
    bkg_plink_prefix::Union{Nothing, String}=nothing,
    ldscore_file::Union{Nothing, String}=nothing,
    maf_match::Union{Nothing, Float64}=nothing,
    ld_match::Union{Nothing, Float64}=nothing,
    threshold::Union{Nothing, Float64}=nothing,
    n_perms::Int=1000,
    )
    println_to_file("Loading and dealing with target dataset...", log_file)
    use_region = false
    target_file_format = "bed"
    @timeit to "Load target" if endswith(target_file, r"bed|bed.gz")
        target_data = CSV.read(target_file, DataFrame, header = false, buffer_in_memory=true)
        target_data[:,1] .= replace.(target_data[:,1],r"^chr"=>"")
        target_data[:,2] .+= 1 
        target_data.len .= target_data[:,3] .- target_data[:,2] .+ 1
        use_region = any(target_data[:,3] .- target_data[:,2] .> 1)
        if !isnothing(maf_match)
            maf_match = nothing
            @warn "The option --maf-match is invalid when using .bed file."
        end
    else
        target_file_format = "omiga"
        target_data = read_cis_file(target_file; threshold)
    end
    begin
        println_to_file("Loading and dealing with background dataset...", log_file)
        bkg_data = DataFrame()
        @timeit to "Load bkg" if isnothing(ldscore_file) 
            bkg_data = read_bkg_bim(bkg_plink_prefix)
        else
            df_ldscore = CSV.File(ldscore_file, types=Dict(1 => String, 2 => String)) |> DataFrame
            bkg_data = DataFrame(:chr => df_ldscore.chr, :start => df_ldscore.bp, :end => df_ldscore.bp, :variant_id => df_ldscore.SNP, :maf => df_ldscore.MAF, :ldscore => df_ldscore.ldscore)
            println_to_file("[INFO] LD score file loaded, and the 'ldscore' column will be used for LD-match analysis.", log_file)
        end
        if target_file_format == "omiga"
            target_data = add_chr_pos_to_cis_tops(target_data, bkg_data)
            target_data[:,1] .= replace.(target_data[:,1],r"^chr"=>"")
        end
        n_rows_target = size(target_data, 1)
        target_data_IC = getIntervalCollection(target_data)
    end
    begin
        annotation_IC = Vector{GenomicIntervalCollection}(undef, length(annot_file_list))
        @threads for i in eachindex(annot_file_list)
            annotation_IC[i] = getIntervalCollection(read_annot_bed(annot_file_list[i]))
        end
        println_to_file("[INFO] $(length(annot_file_list)) annotation file(s) loaded.", log_file)
    end
    println_to_file("Performing $(n_perms) times permutation-based enrichment analysis, start at: $(now())", log_file)
    target_overlap_result = DataFrame()
    for i in eachindex(annot_file_list)
        chromatin_category_data = annotation_IC[i]
        target_overlap_result_1 = getInterval(target_data_IC, chromatin_category_data)
        target_overlap_result_1.bed_file .= basename(annot_file_list[i])
        target_overlap_result = vcat(target_overlap_result, target_overlap_result_1)
    end
    target_overlap_result = DataFrames.combine(groupby(target_overlap_result, [:metadata_B, :bed_file]), nrow => :count)
    target_overlap_result.prop = target_overlap_result.count ./ n_rows_target
    @timeit to "Permutation" permuted_results = enrichment_permutation(n_perms, bkg_data, target_data, annotation_IC, annot_file_list; maf_match=maf_match, ld_match=ld_match, use_region=use_region)
    enrichment_result = calcu_enrichment_p(target_overlap_result, permuted_results, n_perms)
    println_to_file(" * Enrichment analysis completed.", log_file)
    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".enrich.csv")), enrichment_result, compress=false)
    if _args_debug
        println_to_file(string(to), log_file)
        println()
    end
    return enrichment_result
end
