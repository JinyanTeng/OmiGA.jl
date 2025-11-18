function get_samples_from_fam(prefix_plink_file::String)
    fam_path = prefix_plink_file * ".fam"
    fam_file = CSV.read(fam_path, DataFrame, header=false)
    samples = size(fam_file, 1)
    return samples
end
function read_bed_to_UInt8(prefix_plink_file::String)
    bed_path = prefix_plink_file * ".bed"
    samples = get_samples_from_fam(prefix_plink_file) 
    rows = cld(samples, 4) 
    io = open(bed_path, "r")
    seek(io, 3) 
    remaining_size = filesize(io) - position(io) 
    gene_coded_num = remaining_size * 4 
    snps::Int64 = gene_coded_num / (rows * 4)
    bed_data = Mmap.mmap(io, Matrix{UInt8}, (rows, snps))
    close(io)
    return bed_data
end
function count_10s_UInt8(A::UInt8)
    sum_10 = 0
    @inbounds for i in [1:2, 3:4, 5:6, 7:8]
        string(A, pad=8, base=2)[i] == "10" ? sum_10 += 1 : sum_10 += 0
    end
    return sum_10
end
function count_10s_het(prefix_plink_file::String, kept_samples::Union{Nothing, Vector{Int64}})
    data = read_bed_to_UInt8(prefix_plink_file)
    rows, cols = size(data)
    _10s_in_bed = zeros(Int, cols)
    COUNT_10_TABLE = [count_10s_UInt8(x) for x in UInt8(0):UInt8(255)]
    if !isnothing(kept_samples)
        MASK_TABLE = let
            masks = zeros(UInt8, 16)  
            for config in 0:15  
                mask = 0b11111111  
                (config & 1 == 0) && (mask &= 0b11111100)  
                (config & 2 == 0) && (mask &= 0b11110011)  
                (config & 4 == 0) && (mask &= 0b11001111)  
                (config & 8 == 0) && (mask &= 0b00111111)  
                masks[config + 1] = mask  
            end
            masks
        end
        row_configs = zeros(UInt8, rows)
        @inbounds for id in kept_samples
            row_idx = (id - 1) ÷ 4 + 1
            frac_index = (id - 1) % 4  
            row_configs[row_idx] |= (0x01 << frac_index)
        end
        involved_rows = findall(!iszero, row_configs)
        data_cols = [@view data[:, col] for col in 1:cols]
        @batch per=core for col in 1:cols
            sum_10 = 0
            col_data = data_cols[col]
            @inbounds @simd for i in involved_rows
                config = row_configs[i]
                mask = col_data[i] & MASK_TABLE[config + 1]
                sum_10 += COUNT_10_TABLE[mask + 1]
            end
            _10s_in_bed[col] = sum_10
        end
    else 
        data_cols = [@view data[:, col] for col in 1:cols]
        @inbounds @simd for col in 1:cols
            tmpcol = data_cols[col]
            sum_10 = 0
            @inbounds @simd for row in 1:rows 
                index_to_find = tmpcol[row]
                @fastmath sum_10 += COUNT_10_TABLE[index_to_find + 1]
            end
            _10s_in_bed[col] = sum_10
        end
    end
    return _10s_in_bed
end
function cal_het(prefix_plink_file::String, counts_10::Vector{Int64}, kept_samples::Union{Nothing, Vector{Int64}})::Vector{Float64}
    data = read_bed_to_UInt8(prefix_plink_file)
    rows, cols = size(data)
    result = zeros(Float64, cols)
    het = 0.0
    if !isnothing(kept_samples)
        rows = length(kept_samples)
        @fastmath indis = rows
        @inbounds @simd for col in 1:cols
            @fastmath het = Float64(counts_10[col]) / indis
            result[col] = het
        end
    else
        samples = get_samples_from_fam(prefix_plink_file) 
        @inbounds @simd for col in 1:cols
            @fastmath het = Float64(counts_10[col]) / samples
            result[col] = het
        end
    end
    return result
end
function get_allele_het_rate(prefix_plink_file::String; kept_samples::Union{Nothing, Vector{Int64}} = nothing)
    counts_10 = count_10s_het(prefix_plink_file, kept_samples)
    het_result = cal_het(prefix_plink_file, counts_10, kept_samples)
    return het_result
end
function count_1s(prefix_plink_file::String, kept_samples::Union{Nothing, Vector{Int64}})
    data = read_bed_to_UInt8(prefix_plink_file)
    rows, cols = size(data)
    _1s_in_bed = zeros(Int, cols)
    COUNT_BITS_TABLE = [count_ones(x) for x in UInt8(0):UInt8(255)]
    if !isnothing(kept_samples)
        MASK_TABLE = let
            masks = zeros(UInt8, 16)  
            for config in 0:15  
                mask = 0b11111111  
                (config & 1 == 0) && (mask &= 0b11111100)  
                (config & 2 == 0) && (mask &= 0b11110011)  
                (config & 4 == 0) && (mask &= 0b11001111)  
                (config & 8 == 0) && (mask &= 0b00111111)  
                masks[config + 1] = mask  
            end
            masks
        end
        row_configs = zeros(UInt8, rows)
        @inbounds for id in kept_samples
            row_idx = (id - 1) ÷ 4 + 1
            frac_index = (id - 1) % 4  
            row_configs[row_idx] |= (0x01 << frac_index)
        end
        involved_rows = findall(!iszero, row_configs)
        data_cols = [@view data[:, col] for col in 1:cols]
        @batch per = core for col in 1:cols
            sum_bits = 0
            col_data = data_cols[col]
            @inbounds @simd for i in involved_rows
                config = row_configs[i]
                mask = col_data[i] & MASK_TABLE[config + 1]
                sum_bits += COUNT_BITS_TABLE[mask + 1]
            end
            _1s_in_bed[col] = sum_bits
        end
    else
        data_cols = [@view data[:, col] for col in 1:cols]
        for col in 1:cols 
            @inbounds begin
                tmpvar = data_cols[col]
                sum_bits = 0
                @simd for row in 1:rows 
                    sum_bits += COUNT_BITS_TABLE[tmpvar[row] + 1]
                end
                _1s_in_bed[col] = sum_bits
            end
        end
    end
    return _1s_in_bed
end
function cal_maf(prefix_plink_file::String, counts_1::Vector{Int64}, kept_samples::Union{Nothing, Vector{Int64}})::Vector{Float64}
    data = read_bed_to_UInt8(prefix_plink_file)    
    cols = size(data, 2)
    result = zeros(Float64, cols)
    maf = 0.0
    if !isnothing(kept_samples)
        div = length(kept_samples) * 2
        for col in 1:cols
            maf = Float64(counts_1[col]) / div
            maf = maf > 0.5 ? 1 - maf : maf
            result[col] = maf
        end
    else
        samples = get_samples_from_fam(prefix_plink_file) 
        div = samples * 2
        for col in 1:cols
            maf = Float64(counts_1[col]) / div
            maf = maf > 0.5 ? 1 - maf : maf
            result[col] = maf
        end
    end
    return result
end
function get_allele_maf(prefix_plink_file::String; kept_samples::Union{Nothing, Vector{Int64}} = nothing)
    counts_1 = count_1s(prefix_plink_file, kept_samples)
    maf_result = cal_maf(prefix_plink_file, counts_1, kept_samples)
    return maf_result
end
function cal_af(prefix_plink_file::String, counts_1::Vector{Int64}, kept_samples::Union{Nothing, Vector{Int64}})::Vector{Float64}
    data = read_bed_to_UInt8(prefix_plink_file)    
    cols = size(data, 2)
    result = zeros(Float64, cols)
    af = 0.0
    if !isnothing(kept_samples)
        div = length(kept_samples) * 2
        for col in 1:cols
            af = Float64(counts_1[col]) / div
            result[col] = 1 - af
        end
    else
        samples = get_samples_from_fam(prefix_plink_file) 
        div = samples * 2
        for col in 1:cols
            af = Float64(counts_1[col]) / div
            result[col] = 1 - af
        end
    end
    return result
end
function get_allele_af(prefix_plink_file::String; kept_samples::Union{Nothing, Vector{Int64}} = nothing)
    counts_1 = count_1s(prefix_plink_file, kept_samples)
    af_result = cal_af(prefix_plink_file, counts_1, kept_samples)
    return af_result
end
