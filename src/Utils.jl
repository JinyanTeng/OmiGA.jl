function read_qtl_pairs_file(fn::String)
    if endswith(fn, r".txt|.txt.gz")
        return CSV.File(fn) |> DataFrame
    elseif endswith(fn, r".jld2|.jld")
        return JLD2.load(fn)["data"]
    elseif endswith(fn, "arrow")
        return Arrow.Table(fn) |> DataFrame
    else
        error("Unknown file type.")
    end
end
function write_qtl_pairs_file(fn::String, dt::DataFrame, compress::Bool; compress_algo::Symbol=:zstd)
    if endswith(fn, r".txt|.txt.gz")
        CSV.write(fn, dt, delim="\t", compress=compress)
    elseif endswith(fn, r".jld2|.jld")
        JLD2.save(fn, "data", dt)
    elseif endswith(fn, "arrow")
        if compress
            Arrow.write(fn, dt, compress=compress_algo)
        else
            Arrow.write(fn, dt)
        end
    else
        error("Unknown file type.")
    end
end
function write_qtl_pairs_file(fn::String, dt::DataFrame; format::Union{String, Nothing}=nothing)
    if endswith(fn, r".txt")
        CSV.write(fn, dt, delim="\t")
    elseif endswith(fn, r".txt.gz")
        CSV.write(fn, dt, delim="\t", compress=true)
    elseif endswith(fn, r".jld2|.jld")
        if isnothing(format) || format in ["jld","jld2"]
            JLD2.save(fn, "data", dt)
        elseif format in ["jld_compress","jld2_compress"]
            JLD2.save(fn, "data", dt; compress=true)
        end
    elseif endswith(fn, "arrow")
        if isnothing(format) || format == "arrow"
            Arrow.write(fn, dt)
        elseif format == "arrow_zstd"
            Arrow.write(fn, dt, compress=:zstd)
        elseif format == "arrow_lz4"
            Arrow.write(fn, dt, compress=:lz4)
        end
    else
        error("Unknown file type.")
    end
end
function json_write(fn::AbstractString, data)
    open(fn, "w") do f
        JSON.print(f, data, 2)  
    end
end
function nanmaximum(arr)
    filtered = filter(!isnan, arr)
    isempty(filtered) ? NaN : maximum(filtered)
end
function nanminimum(arr)
    filtered = filter(!isnan, arr)
    isempty(filtered) ? NaN : minimum(filtered)
end
function find_closest_below_fast(vec::AbstractVector{T}, target::T) where T
    closest_val = -Inf
    found = false
    i = 0
    final_i = -1
    if issorted(vec)
        for x in vec
            i += 1
            if x <= target && x > closest_val
                closest_val = x
                found = true
                final_i = i
            else
                break
            end
        end
    else
        for x in vec
            i += 1
            if x <= target && x > closest_val
                closest_val = x
                found = true
                final_i = i
            end
        end
    end
    return final_i
end
function find_closest_up_fast(vec::AbstractVector{T}, target::T) where T
    closest_val = Inf
    found = false
    i = 0
    final_i = -1
    if issorted(vec)
        for x in vec
            i += 1
            if x >= target && x < closest_val
                closest_val = x
                found = true
                final_i = i
            else
                break
            end
        end
    else
        for x in vec
            i += 1
            if x >= target && x < closest_val
                closest_val = x
                found = true
                final_i = i
            end
        end
    end
    return final_i
end
function find_closest_below_fast(df::DataFrame, target::T) where T
    pids = unique(df.pheno_id)
    df_pval_threshold_acat = DataFrame(pheno_id = pids, pval_threshold = NaN)
    @threads for i in 1:length(pids)
        useidx = df.pheno_id .== pids[i]
        good_i = find_closest_below_fast(df.cutoff_acatp[useidx], target)
        if good_i == -1
            continue
        else
            if sum(useidx) > good_i
                df_pval_threshold_acat.pval_threshold[i] = sqrt(df.cutoff_pval[useidx][good_i] * df.cutoff_pval[useidx][good_i + 1])
            else
                df_pval_threshold_acat.pval_threshold[i] = df.cutoff_pval[useidx][good_i] + 1e-10
            end
        end
    end
    return df_pval_threshold_acat
end
function runBE(X::AbstractMatrix{T}; B::Int=20, alpha::Float64=0.05, verbose::Bool=false) where T
    @assert alpha >= 0 && alpha <= 1 "alpha must be between 0 and 1."
    n = size(X, 1)  
    p = size(X, 2)  
    d = min(n, p)   
    if verbose
        println("Running PCA on permuted data...")
    end
    testStatsPerm = zeros(d,B)
    @threads for b in 1:B
        if verbose
            println("b=$b out of $B permutations...")
        end
        XPermuted = similar(X)
        for j in 1:p
            XPermuted[:, j] .= shuffle(X[:, j])
        end
        prcompResultPerm = fit(PCA, XPermuted,pratio=1)  
        importanceTablePerm = prcompResultPerm.prinvars 
        PVEsPerm = importanceTablePerm / sum(importanceTablePerm) 
        testStatsPerm[:,b] .= PVEsPerm
    end
    if verbose
        println("Running PCA on the unpermuted data...")
    end
    prcompResult = fit(PCA, X, pratio=1)  
    importanceTable = prcompResult.prinvars
    PVEs = importanceTable / sum(importanceTable)
    pValues = (sum(testStatsPerm .≥ PVEs, dims=2) .+ 1) ./ (B + 1)
    for j in 2:d
        if pValues[j] < pValues[j-1]
            pValues[j] = pValues[j-1]
        end
    end
    numOfPCsChosen = count(p -> p ≤ alpha, pValues)
    return prcompResult, numOfPCsChosen, pValues
end
function runBE(X::AbstractMatrix{T}, prcompResult::PCA{T}; B::Int=20, alpha::Float64=0.05, verbose::Bool=false) where T
    @assert alpha >= 0 && alpha <= 1 "alpha must be between 0 and 1."
    n = size(X, 1)  
    p = size(X, 2)  
    d = length(prcompResult.prinvars)   
    if verbose
        println("Running PCA on permuted data...")
    end
    testStatsPerm = zeros(d, B)
    @threads for b in 1:B
        if verbose
            println("b=$b out of $B permutations...")
        end
        XPermuted = similar(X)
        for j in 1:p
            XPermuted[:, j] .= shuffle(X[:, j])
        end
        prcompResultPerm = fit(PCA, XPermuted,pratio=1)  
        importanceTablePerm = prcompResultPerm.prinvars 
        PVEsPerm = importanceTablePerm / sum(importanceTablePerm) 
        testStatsPerm[1:length(PVEsPerm),b] .= PVEsPerm
    end
    if verbose
        println("Running PCA on the unpermuted data...")
    end
    importanceTable = prcompResult.prinvars
    PVEs = importanceTable / sum(importanceTable)
    pValues = (sum(testStatsPerm .≥ PVEs, dims=2) .+ 1) ./ (B + 1)
    for j in 2:d
        if pValues[j] < pValues[j-1]
            pValues[j] = pValues[j-1]
        end
    end
    numOfPCsChosen = count(p -> p ≤ alpha, pValues)
    return numOfPCsChosen, pValues
end
function runElbow(X::AbstractMatrix{T}) where T
    prcompResult = fit(PCA, X, pratio=1)
    importanceTable = prcompResult.prinvars
    x = 1:length(importanceTable)  
    y = importanceTable / sum(importanceTable)  
    x1 = x[1]  
    y1 = y[1]  
    x2 = x[end]  
    y2 = y[end]  
    x0 = x
    y0 = y
    distancesDenominator = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    distancesNumerator = abs.((x2 - x1) .* (y1 .- y0) .- (x1 .- x0) .* (y2 - y1))
    distances = distancesNumerator / distancesDenominator
    numOfPCsChosen = argmax(distances)  
    return prcompResult, numOfPCsChosen, distances
end
function runElbow(prcompResult::PCA{T}) where T
    importanceTable = prcompResult.prinvars
    x = 1:length(importanceTable)  
    y = importanceTable / sum(importanceTable)  
    x1 = x[1]  
    y1 = y[1]  
    x2 = x[end]  
    y2 = y[end]  
    x0 = x
    y0 = y
    distancesDenominator = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    distancesNumerator = abs.((x2 - x1) .* (y1 .- y0) .- (x1 .- x0) .* (y2 - y1))
    distances = distancesNumerator / distancesDenominator
    numOfPCsChosen = argmax(distances)  
    return numOfPCsChosen, distances
end
function calculate_vif(X::AbstractMatrix; add_intercept::Bool=true)
    n_features = size(X, 2)
    if add_intercept
        X = hcat(ones(size(X, 1)), X)
    end
    vif_values = zeros(n_features)
    for i in 1:n_features
        y = X[:, i+add_intercept]
        X_temp = X[:, setdiff(1:n_features+add_intercept, i+add_intercept)]
        model = lm(X_temp, y)
        r_squared = r²(model)
        vif_values[i] = 1 / (1 - r_squared)
    end
    vif_df = DataFrame(Feature = 1:n_features, VIF = vif_values)
    return sort(vif_df, :VIF, rev=true) 
end
function check_chr_order(chrstr::AbstractVector{T}) where {T}
    is_need_reorder = false
    chroms = unique(chrstr)
    for ch in chroms
        b = findall(chrstr .== ch)
        if (maximum(b) - minimum(b) + 1) > length(b)
            is_need_reorder = true
            break
        end
    end
    return is_need_reorder
end
function range_vector_intersection_set(ranges::AbstractArray{<:AbstractRange}, vec::AbstractVector)
    s = Set(vec)  
    return [any(x -> x in s, r) for r in ranges]
end
function vec_in_ranges(vec::AbstractVector, ranges::AbstractVector{<:AbstractRange})
    result = Vector{Bool}(undef, length(vec))
    for (i, x) in enumerate(vec)
        result[i] = any(r -> x in r, ranges)
    end
    return result
end
function compare_chromosomes(a, b)
    get_num = s -> begin
        num_match = match(r"\d+", s)
        if num_match !== nothing
            parsed = tryparse(Int, num_match.match)
            return parsed === nothing ? typemax(Int) : parsed
        else
            return typemax(Int)
        end
    end
    get_num(a) < get_num(b)
end
function count_max_consecutive_true(vec)  
    count = 0  
    max_count = 0  
    for v in vec  
        if v  
            count = count + 1  
            max_count = max(max_count, count)  
        else  
            count = 0  
        end  
    end  
    return max_count  
end  
function rematch_sample_id!(iid::Vector{T}, id_map::DataFrame) where {T}
    matchidx = repeat([nothing], length(iid)) |> Vector{Union{Nothing, Int64}}
    for i in axes(id_map,2)
        _matchidx = vmatch(id_map[:,i],iid)
        _matchedidx = .!isnothing.(_matchidx)
        erroridx = matchidx[.!isnothing.(matchidx) .& _matchedidx] .!= _matchidx[.!isnothing.(matchidx) .& _matchedidx]
        if any(erroridx)
            error(string(iid[erroridx], " appear multiple times in the id map file!"))
        else
            matchidx[.!isnothing.(_matchidx)] .= _matchidx[.!isnothing.(_matchidx)]
        end
    end
    if any(isnothing.(matchidx))
        @warn string("Cannot found ", iid[isnothing.(matchidx)], " in the id map file!")
        iid[.!isnothing.(matchidx)] .= id_map[matchidx[.!isnothing.(matchidx)],1]
    else
        iid .= id_map[matchidx,1]
    end
    return iid
end
function first_nonnan_argmin(arr::AbstractArray{T}) where {T}
    min_val = Inf
    min_idx = 0
    for (idx, val) in enumerate(arr)
        if !isnan(val) && val < min_val
            min_val = val
            min_idx = idx
        end
    end
    if min_idx == 0
        return -1
    end
    return min_idx
end
function first_nonnan_argmax(arr::AbstractArray{T}) where {T}
    max_val = -Inf
    max_idx = 0
    for (idx, val) in enumerate(arr)
        if !isnan(val) && val > max_val
            max_val = val
            max_idx = idx
        end
    end
    if max_idx == 0
        return -1
    end
    return max_idx
end
function extract_max_by_group(groups, values)
    group_dict = Dict()
    for (group, value) in zip(groups, values)
        if !haskey(group_dict, group) || value > group_dict[group]
            group_dict[group] = value
        end
    end
    max_values = [group_dict[group] for group in groups]
    unique_max_values = unique(max_values)
    return unique_max_values
end
function count_consecutive(arr::Vector; type::String="first")
    if type == "first"
        counts = []
        first_indices = [1]
        count = 1
        for i in 2:length(arr)
            if arr[i] == arr[i-1]
                count += 1
            else
                push!(counts, count)
                count = 1
                push!(first_indices, i)
            end
        end
        push!(counts, count)
        return counts, first_indices
    elseif type == "last"
        isempty(arr) && return [], []  
        counts = Int[]
        last_indices = Int[]
        count = 1
        for i in 2:length(arr)
            if arr[i] == arr[i-1]
                count += 1
            else
                push!(counts, count)
                push!(last_indices, i - 1)
                count = 1  
            end
        end
        push!(counts, count)
        push!(last_indices, length(arr))
        return counts, last_indices
    end
end
function is_collinear(v1, v2, cor_cutoff)
    abscor = abs(cor(v1, v2))  
    return abscor > cor_cutoff  
end
function remove_collinear_columns(matrix::AbstractArray; cor_cutoff::AbstractFloat=0.9, return_index::Bool=false)
    num_cols = size(matrix, 2)  
    cols_to_keep = [1]  
    for i in 2:num_cols
        is_new_col_collinear = any(j -> is_collinear(matrix[:, i], matrix[:, j], cor_cutoff), cols_to_keep)
        if !is_new_col_collinear
            push!(cols_to_keep, i)  
        end
    end
    if return_index
        return cols_to_keep
    else
        return matrix[:, cols_to_keep]  
    end
end
function remove_multicollinearity(matrix::AbstractArray; vif_cutoff::AbstractFloat=10.0, return_index::Bool=false)
    cols_to_keep = 1:size(matrix, 2)
    vif_df = calculate_vif(matrix; add_intercept=false)
    while any(vif_df.VIF .>= vif_cutoff)
        cols_to_keep = cols_to_keep[setdiff(1:length(cols_to_keep), vif_df.Feature[1])]
        num_cols = size(matrix, 2)
        matrix = matrix[:,setdiff(1:num_cols, vif_df.Feature[1])]
        vif_df = calculate_vif(matrix; add_intercept=false)
    end
    if return_index
        return cols_to_keep
    else
        return matrix
    end
end
function create_design_matrix(categories::Vector)
    unique_categories = unique(categories)
    design_matrix = zeros(Int, length(categories), length(unique_categories))
    for (i, category) in enumerate(categories)
        category_index = findfirst(x -> x == category, unique_categories)
        design_matrix[i, category_index] = 1
    end
    design_matrix = design_matrix[:, 1:(end-1)]
    return design_matrix
end
function create_design_matrix(categories_mat::Matrix)
    nc = size(categories_mat, 2)
    c_design_matrix = nothing
    for i in 1:nc
        categories = categories_mat[:, i]
        design_matrix = create_design_matrix(categories)
        if i == 1
            c_design_matrix = design_matrix
        else
            c_design_matrix = [c_design_matrix design_matrix]
        end
    end
    return c_design_matrix
end
function ranking(x::AbstractVector; ties_method="average")
    sorted_indices = sortperm(x)
    ranks = similar(x, Float64)
    current_rank = 1.0
    i = 1
    while i <= length(x)
        start = i
        while i < length(x) && x[sorted_indices[i]] == x[sorted_indices[i+1]]
            i += 1
        end
        if ties_method == "average"
            rank_value = (start + i) / 2.0
        elseif ties_method == "min"
            rank_value = start
        elseif ties_method == "max"
            rank_value = i
        elseif ties_method == "first"
            rank_value = start
        elseif ties_method == "last"
            rank_value = i
        elseif ties_method == "random"
            rank_value = start + rand() * (i - start)
        else
            error("Unsupported ties.method: $ties_method")
        end
        for j in start:i
            ranks[sorted_indices[j]] = rank_value
        end
        current_rank += (i - start + 1)
        i += 1
    end
    return ranks
end
function inverse_normal_transform(x::AbstractVector)
    ranks = ranking(x)
    norm_cdf = (ranks) ./ (length(x) + 1)
    quantile.(Normal(), norm_cdf)
end
function inverse_normal_transform(x::AbstractMatrix)
    x0 = similar(x)
    for i in 1:size(x, 2)
        x0[:, i] .= inverse_normal_transform(x[:, i])
    end
    return x0
end
function inverse_normal_transform!(x::AbstractMatrix)
    for i in 1:size(x, 2)
        x[:, i] .= inverse_normal_transform(x[:, i])
    end
end
function At_mul_B(A, B)
    return permutedims(A) * B
end
function At_mul_B!(C, A, B)
    mul!(C, permutedims(A), B)
end
macro runif(cond, expr)
    quote
        if $(esc(cond))
            $(esc(expr))
        end
    end
end
function df_to_float32!(df)
    cols_float64 = [eltype(df[:, x]) == Float64 for x in range(1, ncol(df))]
    df[!, cols_float64] = Float32.(df[!, cols_float64])
end
function df_to_32bit!(df)
    cols_float64 = [eltype(df[:, x]) == Float64 for x in range(1, ncol(df))]
    cols_int64 = [eltype(df[:, x]) == Int64 for x in range(1, ncol(df))]
    df[!, cols_float64] = Float32.(df[!, cols_float64])
    df[!, cols_int64] = Int32.(df[!, cols_int64])
end
function dict_to_txt(my_dict, fn)
    file = open(fn, "w")
    for (key, value) in my_dict
        write(file, "$key $value\n")
    end
    close(file)
end
function write_pars_to_txt(my_dict, fn)
    file = open(fn, "w")
    write(file, "Parameters for MTGA:\n")
    for (key, value) in my_dict
        write(file, "  --$key $value\n")
    end
    close(file)
end
function println_to_file(message::String, filename::String)
    open(filename, "a") do file
        println(message)
        println(file, message)
        flush(file)
    end
end
function print_to_file(message::String, filename::String)
    open(filename, "a") do file
        print(message)
        print(file, message)
        flush(file)
    end
end
function format_milliseconds(milliseconds::Millisecond)
    total_seconds = milliseconds.value ÷ 1000
    hours = total_seconds ÷ 3600
    minutes = (total_seconds ÷ 60) % 60
    seconds = total_seconds % 60
    milliseconds = milliseconds.value % 1000
    return string(hours, ":", minutes, ":", seconds, ":", milliseconds)
end
function vmatch(S, V; rm_nothing=false)
    idx = indexin(V, S)
    if rm_nothing
        idx = idx[.~isnothing.(idx)] |> Vector{Int}
    end
    return idx
end
function vin(S, V)
    idx = [findfirst(isequal(x), S) for x in V]
    idx = Vector{Union{Nothing,Int}}(idx)
    inn = .~isnothing.(idx)
    return inn
end
function lower_tri(n::Int)
    n_ = Int(n * (n + 1) / 2 - n)
    _idxs = repeat([0], n_)
    _idx = 2:n
    sum_idxs = length(_idx)
    _idxs[1:sum_idxs] = _idx
    for i in 1:(n-1)
        _idx = _idx[2:end] .+ n
        sum_idxs += length(_idx)
        _idxs[sum_idxs-length(_idx)+1:sum_idxs] = _idx
    end
    return _idxs
end
function upper_tri(n::Int)
    n_ = Int(n * (n + 1) / 2 - n)
    _idxs = repeat([0], n_)
    for i in 2:n
        _idx = Int((i - 1) * n + 1):Int((i - 1) * n + (i - 1))
        _idxs[Int((i - 1) * (i - 2) / 2 + 1):Int((i) * (i - 1) / 2)] .= _idx
    end
    return _idxs
end
function writeGRMBin(REL, REL_id, prefix::String, type::Type=Float32; fam_id=nothing)
    function sum_i(i)
        return sum(1:i)
    end
    n = size(REL, 1)
    i = sum_i.(1:n)
    out_REL = Array{type}(repeat([NaN], maximum(i)))
    out_REL[i] = diag(REL)
    idx_up = upper_tri(n)
    out_REL[isnan.(out_REL)] = REL[idx_up]
    open(string(prefix, ".grm.bin"), "w") do f
        write(f, out_REL)
    end
    id = DataFrame(hcat(ifelse(isnothing(fam_id), REL_id, fam_id), REL_id), :auto)
    CSV.write(string(prefix, ".grm.id"), id, delim=" ", header=false)
    return true
end
function readGRMBin(prefix::String; type::Type=Float32, return_type::Type=Float64)
    function sum_i(i)
        return sum(1:i)
    end
    BinFileName = string(prefix, ".grm.bin")
    IDFileName = string(prefix, ".grm.id")
    if !isfile(BinFileName)
        error(string("No found '", BinFileName, "'!"))
    end
    if !isfile(IDFileName)
        error(string("No found '", IDFileName, "'!"))
    end
    id = CSV.read(IDFileName, DataFrame, header=false, types=Dict(:Column1 => String, :Column2 => String)) 
    n = size(id, 1)
    i = sum_i.(1:n)
    i_max = maximum(i)
    grm = zeros(type, i_max)
    open(BinFileName, "r") do f
        grm = Mmap.mmap(f, Vector{type}, i_max, position(f))
    end
    out_REL = Array{return_type}(undef, (n, n))
    out_REL[1:n+1:n*n] = grm[i]
    seq_grm = 1:length(grm)
    tmp = repeat([-1], length(grm))
    tmp[i] .= i
    off = seq_grm[seq_grm.!=tmp]
    out_REL[upper_tri(size(out_REL, 1))] = grm[off]
    out_REL = out_REL'
    out_REL[upper_tri(size(out_REL, 1))] = grm[off]
    return Array(out_REL), id.Column2, id.Column1 
end
function readCSV(file::String, astype::Type; delim::String=" ", rownames::Bool=false, skip_rows::Int=0, read_rows::Int=-1)
    if ~isnothing(match(r".gz", file))
        stream = GZip.open(file)
    else
        stream = open(file, "r")
    end
    io = stream
    len_io = -1
    if skip_rows > 0
        if len_io == -1
            len_io = countlines(io)
            seekstart(io)
        end
        if skip_rows >= len_io
            error("The number of lines skipped beyond the total length of IO.")
        end
        for i_skip in 1:skip_rows
            readline(io)
        end
    end
    if read_rows > 0
        if len_io == -1
            len_io = countlines(io)
            seekstart(io)
        end
        read_rows = (read_rows + skip_rows > len_io) ? (len_io - skip_rows) : read_rows
        strs = Vector{String}(undef, read_rows)
        for i_read in 1:read_rows
            strs[i_read] = readline(io)
        end
    else
        strs = readlines(io)
    end
    close(stream)
    n_rows = length(strs)
    n = n_rows - 1
    jobs = Channel{Int}(n)
    function make_jobs(n)
        for i in 1:n
            put!(jobs, i)
        end
    end
    @async make_jobs(n)
    str_split = split(strs[1], delim)
    line_1 = parse.(astype, str_split[(rownames ? 2 : 1):end])
    n_cols = length(line_1)
    mat = SharedArray{astype}((n_cols, n_rows), init=fill(typemax(astype), (n_cols, n_rows)))
    mat[:, 1] .= line_1
    if rownames
        name_line_1 = str_split[1]
        rownames_Channel = Channel(n_rows)
        put!(rownames_Channel, [1, name_line_1])
    end
    function do_work()
        for job_id in jobs
            if rownames
                put!(rownames_Channel, [job_id + 1, strs[job_id+1][1:Int(match(Regex(delim), strs[job_id+1]).offset)-1]])
            end
            mat[:, job_id+1] .= parse.(astype, split(strs[job_id+1], delim)[(rownames ? 2 : 1):end])
        end
    end
    @threads for i in 1:n_rows 
        @async do_work()
    end
    if rownames
        v_rownames = Vector{String}(undef, n_rows)
        for i in 1:rownames_Channel.sz_max
            row_id, rowname = take!(rownames_Channel)
            v_rownames[row_id] = rowname
        end
        return v_rownames, mat'
    end
    return mat'
end
function readCSV(file::String; delim::String=" ")
    io = open(file, "r")
    strs = readlines(io)
    close(io)
    n_rows = length(strs)
    n = n_rows - 1
    jobs = Channel{Int}(n)
    function make_jobs(n)
        for i in 1:n
            put!(jobs, i)
        end
    end
    @async make_jobs(n)
    line_1 = split(strs[1], delim)
    n_cols = length(line_1)
    mat = Array{String}(undef, (n_rows, n_cols))
    mat[1, :] .= line_1
    function do_work()
        @views for job_id in jobs
            tmp = mat[job_id+1, :]
            tmp[:] .= split(strs[job_id+1], delim)
        end
    end
    @threads for i in 1:n_rows 
        @async do_work()
    end
    return mat
end
function readPlinkBed(prefix::String; type::Type=Int8, model=ADDITIVE_MODEL, center::Bool=false, scale::Bool=false, impute::Bool=false, byrow=true, loadfunc=SNP_matrix)
    bed = loadfunc(string(prefix, ".bed"))
    if !byrow
        bed = bed'
    end
    bim = CSV.read(string(prefix, ".bim"), DataFrame, header=false, buffer_in_memory=!_args_low_mem, types=Dict(1 => String, 2 => String, 3 => Float64, 4 => Int, 5 => String, 6 => String))
    rename!(bim, ["chromosome", "variant", "cM", "position", "a1", "a2"])
    fam = CSV.read(string(prefix, ".fam"), DataFrame, header=false, buffer_in_memory=true, types=Dict(1 => String, 2 => String))
    FID = string.(fam[:, 1])::Vector{String}
    IID = string.(fam[:, 2])::Vector{String}
    if type == UInt8
        return bed, IID, FID, bim
    else
        geno = convert(Matrix{type}, bed)::Matrix{type}
        return geno, IID, FID, bim
    end
end
function readPlinkBed(prefix::String, dominance::Bool, byrow::Bool)
    bed = SnpArray(string(prefix, ".bed"))
    if !byrow
        bed = bed' |> Matrix
        replace!(bed, 0x00 => 0x02, 0x02 => 0x01, 0x03 => 0x00) 
    else
        bed = replace(bed, 0x00 => 0x02, 0x02 => 0x01, 0x03 => 0x00) 
    end
    if dominance
        geno_DOM = copy(bed)
        replace!(geno_DOM, 0x02 => 0x00) 
    end
    bim = CSV.read(string(prefix, ".bim"), DataFrame, header=false, buffer_in_memory=!_args_low_mem, types=Dict(1 => String, 2 => String, 3 => Float64, 4 => Int, 5 => String, 6 => String))
    rename!(bim, ["chromosome", "variant", "cM", "position", "a1", "a2"])
    fam = CSV.read(string(prefix, ".fam"), DataFrame, header=false, buffer_in_memory=true, types=Dict(1 => String, 2 => String))
    FID = string.(fam[:, 1])::Vector{String}
    IID = string.(fam[:, 2])::Vector{String}
    if dominance
        return bed, IID, FID, bim, geno_DOM
    else
        return bed, IID, FID, bim
    end
end
function readPlinkPhen(prefix::String)
    phen = CSV.read(string(prefix, ".phen"), DataFrame, header=false, types=Dict(1 => String, 2 => String))
    FID = string.(phen[:, 1])
    IID = string.(phen[:, 2])
    return phen, IID, FID
end
function readVCF(filename::String; type::Type=Int8, code::String="HT", model=:additive, impute::Bool=false, center::Bool=false, scale::Bool=false, trans::Bool=false)
    fh = openvcf(filename, "r")
    ncount = countlines(fh)
    seekstart(fh)
    nrecord = nrecords(vcffile)
    for l in 1:(ncount-nrecord-1)
        readline(fh)
    end
    header = readline(fh)
    close(fh)
    nsample = nsamples(vcffile)
    header_name = split(header)
    ID = header_name[(length(header_name)-nsample+1):end]
    if code == "GT"
        GENO = convert_gt(type, filename, model=model, impute=impute, center=center, scale=scale, trans=trans, msg="Progress: ")
    elseif code == "HT"
        GENO = convert_ht(type, filename, trans=trans, msg="Progress: ")
        ID = repeat(ID, inner=2)
    elseif code == "DS"
        GENO = convert_ds(type, filename, model=model, impute=impute, center=center, scale=scale, trans=trans, msg="Progress: ")
    else
        error("No this code type.")
    end
    return GENO, ID
end
function get_allele_af(genotype::AbstractMatrix; byrow::Bool=true) 
    dims = ifelse(byrow, 1, 2)
    return vec(sum(genotype, dims=dims)) / (2 * size(genotype, dims))
end
function get_allele_af(genotype::AbstractVector) 
    return sum(genotype) / (2 * length(genotype))
end
function get_allele_maf(genotype::AbstractArray; byrow::Bool=true) 
    maf = get_allele_af(genotype, byrow=byrow)
    maf[maf.>0.5] .= 1 .- maf[maf.>0.5]
    return maf
end
function get_allele_het_rate(genotype::AbstractArray; byrow::Bool=true) 
    dims = ifelse(byrow, 1, 2)
    return vec(sum(genotype .== 1, dims=dims)) / size(genotype, dims)
end
function get_allele_het_rate(genotype::SNP_matrix) 
    return het(genotype)
end
function get_allele_af(genotype::SNP_matrix) 
    return af(genotype)
end
function get_allele_maf(genotype::SNP_matrix) 
    return maf(genotype)
end
function batch_genotype_viewfunc!(v::AbstractVector,genotype::AbstractArray,func::Function,chunk_size::Integer)
    n_snps = size(genotype,2)
    batch_size = n_snps >= 1000000 ? chunk_size : n_snps
    batch_collects = collect(Iterators.partition(1:n_snps,batch_size))
    @threads for j in batch_collects
        v[j] = @views func(genotype[:,j])
    end
end
function getPar(parcard, par_name)
    sp = split(parcard[.!isnothing.(findfirst.(string("@", par_name), parcard))][1], " ")
    r_par = sp[2:end]
    return r_par
end
function parse_commandline(args)
    s = ArgParseSettings()
    add_arg_table!(s,
        ["--DEBUG"],
        Dict(
            :help => "Debug program.",
            :arg_type => Bool,
            :default => false
        ),
        ["--SEED", "-s"],
        Dict(
            :help => "Set random seed.",
            :arg_type => Int,
            :default => 1996
        ),
        ["--WORK_DIR", "-w"],
        Dict(
            :help => "Set up working directory. All cache and output files will present in the directory.",
            :arg_type => String,
            :default => "."
        ),
        ["--DATA_DIR"],
        Dict(
            :help => "Set up data directory. Program find and read the data from the directory.",
            :arg_type => String,
            :default => "."
        ),
        ["--DEBV"],
        Dict(
            :help => "Calculate the de-regressed EBVs based on animal model.",
            :arg_type => Bool,
            :default => false
        ),
        ["--RUN_MODE"],
        Dict(
            :help => "Specify the program running mode [Required]. 1 = Estimate breeding values directly, 2 = Perform Cross-Validation.",
            :arg_type => Int,
            :default => 1
        ),
        ["--CV_MODE"],
        Dict(
            :help => "Specify the Cross-Validation mode. 1 = Random Cross-Validation, 2 = Leave-One-Out.",
            :arg_type => Int,
            :default => 1
        ),
        ["--MODEL"],
        Dict(
            :help => "Specify the linear mixed model of genomic prediction [Required]. 1 = Pedigree-based BLUP, 2 = SNP-based GBLUP, 3 = Single-step GBLUP, 99 = User-defined model which requires to specify '--REL'.",
            :arg_type => Int,
            :default => 1
        ),
        ["--REL"],
        Dict(
            :help => "Specify the prefix of relationship matrix file(s) of PLINK format.",
            :nargs => '*',
            :arg_type => String
        ),
        ["--PEDIGREE"],
        Dict(
            :help => "Specify the filename of pedigree required for BLUP (if MODEL = 1) and ssGBLUP (if MODEL = 3).",
            :arg_type => String,
            :default => ""
        ),
        ["--GENOTYPE"],
        Dict(
            :help => "Specify the filename of genotype required for GBLUP (if MODEL = 2) and ssGBLUP (if MODEL = 3).",
            :arg_type => String,
            :default => ""
        ),
        ["--PHENOTYPE"],
        Dict(
            :help => "Specify the filename of phenotype.",
            :nargs => '+',
            :arg_type => String
        ),
        ["--MISS"],
        Dict(
            :help => "Set missing phenotype value.",
            :arg_type => String,
            :default => "-999"
        ),
        ["--COVAR"],
        Dict(
            :help => "Specify covariate filename.",
            :arg_type => String,
            :default => ""
        ),
        ["--DCOVAR"],
        Dict(
            :help => "Specify index of discrete covariate(s) in file.",
            :nargs => '*',
            :arg_type => Int,
            :default => [0]
        ),
        ["--QCOVAR"],
        Dict(
            :help => "Specify index of quantitative covariate(s) in file.",
            :nargs => '*',
            :arg_type => Int,
            :default => [0]
        ),
        ["--FOLD_REP"],
        Dict(
            :help => "Set the number of fold and replicate in random Cross-Validation (required if CV_MODE = 1).",
            :nargs => 2,
            :arg_type => Int
        ),
        ["--GROUP"],
        Dict(
            :help => "Specify filename of grouping of Cross-validation.",
            :arg_type => String,
            :default => ""
        ),
        ["--REML"],
        Dict(
            :help => "Specify the REML method to estimate the variance components. 1 = JuliaPkg [Default], 2 = GCTA, 3 = LDAK. Combine with PATH of used software and its REML_ALGO if value = 2 or 3.",
            :nargs => '*',
            :arg_type => String,
            :default => ["1"]
        ),
        ["--VARS"],
        Dict(
            :help => "Specify the variance components instead of estimating them from the '--REML'.",
            :nargs => '*',
            :arg_type => Float64
        ),
        ["--H2"],
        Dict(
            :help => "Specify the heritability used in model instead of estimating from software.",
            :nargs => '*',
            :arg_type => Float64
        ),
        ["--WEIGHT_A"],
        Dict(
            :help => "Set the ratio of pedigree in H matirx of ssGBLUP.",
            :arg_type => Float64,
            :default => 0.0
        ),
        ["--SAVE"],
        Dict(
            :help => "Save the intermediate file produced by the program.",
            :arg_type => Bool,
            :default => false
        ),
        ["--THREADS"],
        Dict(
            :help => "Specify the number of threads.",
            :arg_type => Int,
            :default => 2
        ),
        ["--OPEN_PCG"],
        Dict(
            :help => "Open the preconditional conjugate gradiem (PCG) algorithm to solve the MME.",
            :arg_type => Bool,
            :default => false
        ),
        ["--REFERENCE"],
        Dict(
            :help => "Specify the filename of reference individuals.",
            :arg_type => String,
            :default => ""
        ),
        ["--CANDIDATE"],
        Dict(
            :help => "Specify the filename of candidate individuals.",
            :arg_type => String,
            :default => ""
        )
    )
    return parse_args(args, s)
end
function writeParstoDisk(PARs_DICT)
    un_str = .~isnothing.([values(PARs_DICT)...]) .& .~isequal.(typeof.([values(PARs_DICT)...]), String)
    max_len = maximum([length(PARs_DICT[x]) for x in [keys(PARs_DICT)...][un_str]])
    n_key = length(keys(PARs_DICT))
    par_save = reshape(repeat([""], n_key * (max_len + 1)), n_key, max_len + 1)
    par_save[:, 1] = string.("@", [keys(PARs_DICT)...])
    for x in 1:n_key
        tmp = [values(PARs_DICT)...][x]
        if isequal(typeof.(tmp), String)
            par_save[x, 2] = tmp
        elseif isnothing(tmp)
            continue
        elseif length(tmp) > 0
            par_save[x, 2:(1+length(tmp))] .= string.(tmp)
        end
    end
    par_save = par_save[(sum(par_save .== "", dims=2).!=max_len)[:, 1], :]
    CSV.write("PARS.DAT", sort(DataFrame(par_save, :auto)), delim=" ", header=false)
end
function CreateGroup(k_rep::Int, n_fold::Int, n_idvs::Int; seed::Int=0)
    group = Matrix{Int}(undef, n_idvs, k_rep + 1)
    group[:, 1] = 1:n_idvs
    for idx in 1:k_rep
        StatsBase.Random.seed!(idx^3 + 666 + seed)
        sample_id = StatsBase.sample(1:n_idvs, n_idvs, replace=false)
        sample_step = 1:Int(floor(n_idvs / n_fold)):n_idvs+1
        for j in 1:n_fold
            group[sample_id[sample_step[j]:(sample_step[j+1]-1)], idx+1] .= j
        end
        n_rest = n_idvs - (maximum(sample_step) - 1)
        if n_rest > 0
            group[sample_id[end-(n_rest-1):end], idx+1] = StatsBase.sample(1:n_fold, n_rest, replace=false)
        end
    end
    return group
end
function getFixedMatrix(factors)
    ts = term.((1, Symbol.(names(factors))...))
    f = term(:1) ~ foldl(+, ts)
    fixedX = modelmatrix(f, factors)
    return fixedX
end
