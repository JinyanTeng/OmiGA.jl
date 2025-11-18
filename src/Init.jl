mutable struct Genotype{T}
    genotype::Union{Nothing,AbstractMatrix{T}}
    annotation::DataFrame
    n_samples::Int
    n_snps::Int
    fid::Union{Nothing,Vector{String}}
    iid::Vector{String}
    code_type::String
    centered::Bool
    kept_samples::Union{Nothing,Vector{Int}}
    kept_variants::Union{Nothing,Vector{Int}}
end
function Genotype(geno_file_prefix::String; mac_threshold::T1=1, maf_threshold::T2=T2(0.0), maf_threshold_interaction::T2=T2(0.0), het_threshold_interaction::T2=T2(0.0), het_rate_threshold::T2=T2(1.0), byrow::Bool=false, ITERM::Union{Nothing,AbstractArray{T2}}=nothing, pre_adj_genotype::Bool=true, keep_raw_geno::Bool=false, extract_samples::Union{Vector{String},Nothing}=nothing, id_map::Union{Nothing, DataFrame}=nothing,low_mem::Bool=false,use_dominance::Bool=false) where {T1<:Real,T2<:AbstractFloat}
    genotype=nothing
    loaded_n_samples = 0
    loaded_n_variants = 0
    if sum(isfile.(string.(geno_file_prefix, [".bed", ".bim", ".fam"]))) == 3
        genotype, iid, fid, annotation = readPlinkBed(geno_file_prefix; byrow=byrow,type=UInt8,loadfunc=SNP_matrix)
        loaded_n_samples = length(iid)
        loaded_n_variants = size(annotation,1)
        GC.gc()
    else
        error("Incorrect genotype file, please check!")
    end
    @runif !isnothing(id_map) rematch_sample_id!(iid, id_map)
    filtered = false
    idxlist_kept_geno_samples = 1:loaded_n_samples
    if !isnothing(extract_samples)
        if length(extract_samples) != loaded_n_samples[ifelse(byrow, 1, 2)]
            idxlist_kept_geno_samples = vmatch(iid, extract_samples;rm_nothing=true)
        else
            if any(extract_samples .!= iid)
                idxlist_kept_geno_samples = vmatch(iid, extract_samples;rm_nothing=true)
            end
        end
        if !isnothing(idxlist_kept_geno_samples)
            if any(isnothing.(idxlist_kept_geno_samples))
                nofound_iid = extract_samples[findall(isnothing.(idxlist_kept_geno_samples))]
                error(string("No found samples ", nofound_iid, " in genotype data!"))
            end
            iid = iid[idxlist_kept_geno_samples]
            fid = fid[idxlist_kept_geno_samples]
        end
    end
    n_samples_g, n_snps = length(idxlist_kept_geno_samples), loaded_n_variants
    idxlist_pass_filter_variants = 1:loaded_n_variants
    @runif _args_debug println("compute af") 
    annotation.af = get_allele_af(geno_file_prefix; kept_samples=loaded_n_samples == n_samples_g ? nothing : idxlist_kept_geno_samples)
    @runif USE_Float32 annotation.af = FloatT.(annotation.af)
    GC.gc()
    if !_args_nofilter
        mafs = copy(annotation.af)
        if sum(mafs .> 0.5) == 0
            println_to_file("\n\t* User provided major allele coding genotype or no variants with allele frequency > 0.5", log_file)
        else
            mafs[mafs.>FloatT(0.5)] .= 1 .- mafs[mafs.>FloatT(0.5)]
        end
        if _args_debug
            println_to_file("\n\t* Description for AF of genotype data:", log_file)
            println_to_file(string(summarystats(annotation.af)), log_file)
            println_to_file("\n\t* Description for MAF of genotype data:", log_file)
            println_to_file(string(summarystats(mafs)), log_file)
        end
        macs = round.(Int, mafs .* 2n_samples_g)
        @runif _args_debug println("compute het_rate")
        annotation.het_rate = get_allele_het_rate(geno_file_prefix; kept_samples=loaded_n_samples == n_samples_g ? nothing : idxlist_kept_geno_samples)
        if _args_debug
            println_to_file("\n\t* Description for Het_rate of genotype data:", log_file)
            println_to_file(string(summarystats(annotation.het_rate)), log_file)
        end
        @runif USE_Float32 annotation.het_rate = FloatT.(annotation.het_rate)
        GC.gc()
        if mac_threshold < 1
            mac_threshold = 1
            @info "Require variants with MAC ≥ 1 for analysis."
        end
        mac_eliminate_index = macs .< mac_threshold
        maf_eliminate_index = mafs .< maf_threshold
        if use_dominance
            het_rate_eliminate_index = (annotation.het_rate .>= het_rate_threshold) .| (annotation.het_rate .<= 1 - het_rate_threshold)
        else
            het_rate_eliminate_index = annotation.het_rate .>= het_rate_threshold
        end
        eliminate_index = mac_eliminate_index .| maf_eliminate_index .| het_rate_eliminate_index
        if any(eliminate_index)
            @runif !use_dominance println_to_file(string("\t[INFO] Eliminated ", sum(eliminate_index), " variants (", sum(mac_eliminate_index), " with MAC < ", mac_threshold, ", ", sum(maf_eliminate_index), " with MAF < ", maf_threshold, " and ", sum(het_rate_eliminate_index), " with Heterozygosity rate ≥ ", het_rate_threshold, ")"), log_file)
            @runif use_dominance println_to_file(string("\t[INFO] Eliminated ", sum(eliminate_index), " variants (", sum(mac_eliminate_index), " with MAC < ", mac_threshold, ", ", sum(maf_eliminate_index), " with MAF < ", maf_threshold, " and ", sum(het_rate_eliminate_index), " with Heterozygosity rate ≥ ", het_rate_threshold, " and/or ≤ ", 1 - het_rate_threshold, ")"), log_file)
        end
        snps_pass_maf_threshold = fill(true,size(annotation,1))
        snps_pass_het_threshold = fill(true,size(annotation,1))
        if (maf_threshold_interaction > 0.0) & !isnothing(ITERM)
            num_elements = size(ITERM, 1)
            half_point = round(Int, num_elements * 0.5)
            top_indices = 1:half_point
            bottom_indices = (half_point+1):num_elements
            for int_i in 1:size(ITERM,2)
                sorted_indices = sortperm(ITERM[:, int_i], rev=true)
                top50_percent_indices = sorted_indices[top_indices] 
                bottom50_percent_indices = sorted_indices[bottom_indices] 
                sort!(top50_percent_indices)
                sort!(bottom50_percent_indices)
                maf_top50per = get_allele_maf(geno_file_prefix; kept_samples=idxlist_kept_geno_samples[top50_percent_indices])
                maf_bottom50per = get_allele_maf(geno_file_prefix; kept_samples=idxlist_kept_geno_samples[bottom50_percent_indices])
                het_top50per = get_allele_het_rate(geno_file_prefix; kept_samples=idxlist_kept_geno_samples[top50_percent_indices])
                het_bottom50per = get_allele_het_rate(geno_file_prefix; kept_samples=idxlist_kept_geno_samples[bottom50_percent_indices])
                snps_pass_maf_threshold .&= (maf_top50per .> maf_threshold_interaction) .& (maf_bottom50per .> maf_threshold_interaction)
                snps_pass_het_threshold .&= (het_top50per .< het_threshold_interaction) .& (het_bottom50per .< het_threshold_interaction)
            end
            println_to_file(string("\t[INFO] Eliminated ", sum(.!(snps_pass_maf_threshold .& snps_pass_het_threshold)), " variants due to '--maf-threshold-interaction ", maf_threshold_interaction, "' and '--het-threshold-interaction ", het_threshold_interaction, "'"), log_file)
        end
    end
    idxlist_pass_filter_variants = idxlist_pass_filter_variants[.!eliminate_index .& snps_pass_maf_threshold .& snps_pass_het_threshold]
    annotation = annotation[idxlist_pass_filter_variants,:]
    annotation.index = range(1, size(annotation, 1))    
    code_type = "additive"
    n_snps = size(annotation, 1)
    @runif USE_Float32 df_to_32bit!(annotation)
    if low_mem
        @views genotype = genotype[idxlist_kept_geno_samples,idxlist_pass_filter_variants]
    else
        genotype = genotype[idxlist_kept_geno_samples,idxlist_pass_filter_variants]
    end
    Genotype{eltype(genotype)}(genotype, annotation, n_samples_g, n_snps, fid, iid, code_type, pre_adj_genotype, idxlist_kept_geno_samples, idxlist_pass_filter_variants)
end
mutable struct Dominance{T}
    genotype::Union{Nothing,AbstractMatrix{T}}
    annotation::DataFrame 
    n_samples::Int
    n_snps::Int
    code_type::String
    keep_index::Union{Nothing,Vector{Bool}}
end
function Dominance(geno012::AbstractMatrix{T}, annotation::DataFrame, byrow::Bool=false; geno_file_prefix::String=nothing) where {T<:Number}
    if _args_low_mem
        domGenotype = SNP_matrix(string(geno_file_prefix,".bed"),coding_rule="010")
        @views domGenotype = domGenotype[geno012.indices[1], geno012.indices[2]]
    else
        domGenotype = copy(geno012)
        replace!(domGenotype, 2 => 0)
    end
    domSnpAnnot = copy(annotation)
    n_samples, n_snps = size(domGenotype)[ifelse(byrow, [1, 2], [2, 1])]
    keep_index = nothing
    code_type = "dominance"
    Dominance(domGenotype, domSnpAnnot, n_samples, n_snps, code_type, keep_index)
end
mutable struct Phenotype{T<:AbstractFloat}
    phenotype::Union{Nothing,Matrix{T}}
    annotation::DataFrame
    iid::Union{Nothing,Vector{String}}
    chroms::Vector{String}
    n_samples::Union{Nothing,Int}
    n_phenotypes::Int
    phenotype_pairs::Union{Nothing,DataFrame}
end
function Phenotype(phenotype_bed_file::String; rm_pheno_threshold::T=0.1, chrom::Vector{String}=String[], extract_pheno_file::Union{String,Nothing}=nothing, extract_pheno_name::Vector{String}=String[], exclude_pheno_file::Union{String,Nothing}=nothing, exclude_pheno_name::Vector{String}=String[], id_map::Union{Nothing, DataFrame}=nothing, extract_samples::Union{Nothing, Vector{String}}=nothing) where {T<:AbstractFloat}
    if isfile(phenotype_bed_file)
        phenotype_bed = CSV.File(phenotype_bed_file, header=true, buffer_in_memory=true, types=Dict(1 => String, 2 => Int, 3 => Int, 4 => String)) |> DataFrame
    else
        error("Phenotype file does not exist, please check!")
    end
    @runif USE_Float32 df_to_32bit!(phenotype_bed)
    iid = names(phenotype_bed)[5:end]
    if length(chrom) > 0
        phenotype_bed = phenotype_bed[phenotype_bed[:, 1].∈(chrom,), :]
    end
    phenotype_pairs = nothing
    if !isnothing(extract_pheno_file)
        extract_pheno_df = CSV.File(extract_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(extract_pheno_df, 2) == 2
            rename!(extract_pheno_df, [:phe1, :phe2])
            extract_pheno_list = unique([extract_pheno_df.phe1; extract_pheno_df.phe2])
            phenotype_pairs = extract_pheno_df
        elseif size(extract_pheno_df, 2) == 1
            rename!(extract_pheno_df, [:phe1])
            extract_pheno_list = unique(extract_pheno_df.phe1)
        else
            error("")
        end
        phenotype_bed = phenotype_bed[phenotype_bed[:, 4].∈(extract_pheno_list,), :]
    end
    if !isnothing(exclude_pheno_file)
        exclude_pheno_df = CSV.File(exclude_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(exclude_pheno_df, 2) == 1
            rename!(exclude_pheno_df, [:phe1])
            extract_pheno_list = unique(exclude_pheno_df.phe1)
        else
            error("")
        end
        phenotype_bed = phenotype_bed[.!(phenotype_bed[:, 4] .∈ (extract_pheno_list,)), :]
    end
    if length(extract_pheno_name) > 0
        phenotype_bed = phenotype_bed[phenotype_bed[:, 4].∈(extract_pheno_name,), :]
    end
    if length(exclude_pheno_name) > 0
        phenotype_bed = phenotype_bed[.!(phenotype_bed[:, 4] .∈ (exclude_pheno_name,)), :]
    end
    annotation = phenotype_bed[:, 1:4]      
    rename!(annotation, ["chrom", "start", "end", "pheno_id"])
    annotation.start .+= 1
    annotation.end .+= 1
    phenotype = permutedims(phenotype_bed[:, 5:end]) |> Matrix{FloatT}
    n_samples_e, n_phenotypes = size(phenotype)      
    keep_index = [sum(phenotype[:, i] .>= median(phenotype[:, i])) for i in 1:n_phenotypes] .>= rm_pheno_threshold * n_samples_e
    if sum(keep_index) < n_phenotypes
        println_to_file(string("\n\t* Warning: Eliminated ", n_phenotypes - sum(keep_index), " phenotypes with < ", rm_pheno_threshold * 100, "% samples greater than median value (the expected value is 50% for normally distributed data)"), log_file)
        elp_file = replace(log_file, r".log$" => ".elp")
        println_to_file(string("\t* Description for the eliminated phenotypes can be found in ", elp_file), log_file)
        describe_df = describe(DataFrame(phenotype[:, .!keep_index], :auto), :mean, :std, :min, :q25, :median, :q75, :max, :nuniqueall)
        describe_df.variable .= annotation.pheno_id[.!keep_index]
        CSV.write(elp_file, describe_df, delim="\t")
        phenotype = phenotype[:, keep_index]
        annotation = annotation[keep_index, :]
        n_samples_e, n_phenotypes = size(phenotype)
        if n_phenotypes == 0
            @error "No phenotype was available for analysis as all phenotypes were eliminated by '--rm-pheno-threshold'. Please check the phenotype data!"
        end
    end
    @runif check_chr_order(annotation.chrom) @error "The phenotype data is in unordered! Please sort it by chromosome and position."
    chroms = string.(unique(annotation.chrom))
    @runif USE_Float32 df_to_32bit!(annotation)
    @runif !isnothing(id_map) rematch_sample_id!(iid, id_map)
    if !isnothing(extract_samples)
        order_row_index = nothing
        if length(extract_samples) != length(iid)
            order_row_index = vmatch(iid, extract_samples,rm_nothing=true)
        else
            if any(extract_samples .!= iid)
                order_row_index = vmatch(iid, extract_samples,rm_nothing=true)
            end
        end
        if !isnothing(order_row_index)
            phenotype = phenotype[order_row_index,:]
            iid = iid[order_row_index]
        end
        n_samples_e, n_phenotypes = size(phenotype)
    end
    if _args_normalized_phenotype
        inverse_normal_transform!(phenotype)
    end
    annotation.index = range(1, size(annotation, 1))    
    Phenotype{eltype(phenotype)}(phenotype, annotation, iid, chroms, n_samples_e, n_phenotypes, phenotype_pairs)
end
function Phenotype(phenotype_file::String, pheno_annot_file::Union{String,Nothing}; rm_pheno_threshold::T=0.1, chrom::Vector{String}=String[], extract_pheno_file::Union{String,Nothing}=nothing, extract_pheno_name::Vector{String}=String[], exclude_pheno_file::Union{String,Nothing}=nothing, exclude_pheno_name::Vector{String}=String[], id_map::Union{Nothing, DataFrame}=nothing, extract_samples::Union{Nothing, Vector{String}}=nothing) where {T<:AbstractFloat}
    if isfile(phenotype_file)
        phenotype_bed = CSV.File(phenotype_file, header=true, buffer_in_memory=true, missingstring=["nan", "na", "NA", "NaN", "NAN", ""], types=Dict(1=>String)) |> DataFrame
    else
        error("Phenotype file does not exist, please check!")
    end
    col_iid = findfirst(.!isnothing.(match.(r"IID$", names(phenotype_bed))))
    phenotype_bed[!,col_iid] = string.(phenotype_bed[:,col_iid])
    iid = string.(phenotype_bed[:, col_iid])
    pheno_id = names(phenotype_bed)[(col_iid+1):end]
    phenotype = phenotype_bed[:, (col_iid+1):end] |> Matrix{Union{Missing, FloatT}}
    replace!(phenotype, missing => NaN)
    n_samples_e, n_phenotypes = size(phenotype)      
    annotation = nothing
    @runif !isnothing(pheno_annot_file) if isfile(pheno_annot_file)
        annotation = CSV.File(pheno_annot_file, header=true, buffer_in_memory=true, types=Dict(1 => String, 2 => Int64, 3 => Int64, 4 => String)) |> DataFrame
        rename!(annotation, ["chrom", "start", "end", "pheno_id"])
    else
        error("Annotation file does not exist, please check!")
    end
    if !isnothing(annotation)
        if size(annotation, 1) != n_phenotypes
            error("Mismatched phenotypes number between phenotype and annotation file!")
        end
        if sum(annotation.pheno_id .!= pheno_id) > 0
            error("Mismatched pheno_id between phenotype and annotation!")
        end
        if length(chrom) > 0
            phenotype = phenotype[:, annotation.chrom.∈(chrom,)]
            annotation = annotation[annotation.chrom.∈(chrom,), :]
        end
    else
        annotation = DataFrame(chrom=repeat(["n.s."], n_phenotypes), start_pos=repeat([-1], n_phenotypes), end_pos=repeat([-1], n_phenotypes), pheno_id=pheno_id)
        rename!(annotation, ["chrom", "start", "end", "pheno_id"])
    end
    annotation.start .+= 1
    annotation.end .+= 1
    phenotype_pairs = nothing
    if !isnothing(extract_pheno_file)
        extract_pheno_df = CSV.File(extract_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(extract_pheno_df, 2) == 2
            rename!(extract_pheno_df, [:phe1, :phe2])
            extract_pheno_list = unique([extract_pheno_df.phe1; extract_pheno_df.phe2])
            phenotype_pairs = extract_pheno_df
        elseif size(extract_pheno_df, 2) == 1
            rename!(extract_pheno_df, [:phe1])
            extract_pheno_list = unique(extract_pheno_df.phe1)
        else
            error("")
        end
        extract_index = annotation.pheno_id .∈ (extract_pheno_list,)
        phenotype = phenotype[:, extract_index]
        annotation = annotation[extract_index, :]
    end
    if !isnothing(exclude_pheno_file)
        exclude_pheno_df = CSV.File(exclude_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(exclude_pheno_df, 2) == 1
            rename!(exclude_pheno_df, [:phe1])
            extract_pheno_list = unique(exclude_pheno_df.phe1)
        else
            error("")
        end
        exclude_index = .!(annotation.pheno_id .∈ (extract_pheno_list,))
        phenotype = phenotype[:, exclude_index]
        annotation = annotation[exclude_index, :]
    end
    if length(extract_pheno_name) > 0
        extract_index = annotation.pheno_id .∈ (extract_pheno_name,)
        phenotype = phenotype[:, extract_index]
        annotation = annotation[extract_index, :]
    end
    if length(exclude_pheno_name) > 0
        exclude_index = .!(annotation.pheno_id .∈ (exclude_pheno_name,))
        phenotype = phenotype[:, exclude_index]
        annotation = annotation[exclude_index, :]
    end
    nonan_index = .!isnan.(phenotype[:, 1])
    phenotype = phenotype[nonan_index, :] |> Matrix{FloatT} 
    iid = iid[nonan_index]
    n_samples_e, n_phenotypes = size(phenotype)      
    keep_index = [sum(phenotype[:, i] .> median(phenotype[:, i])) for i in 1:n_phenotypes] .>= rm_pheno_threshold * n_samples_e
    if sum(keep_index) < n_phenotypes
        println_to_file(string("\tWarning: Eliminated ", n_phenotypes - sum(keep_index), " phenotypes with < ", rm_pheno_threshold * 100, "% samples greater than median value (the expected value is 50% for normally distributed data)"), log_file)
        elp_file = replace(log_file, r".log$" => ".elp")
        println_to_file(string("\tDescription for the eliminated phenotypes can be found in ", elp_file), log_file)
        describe_df = describe(DataFrame(phenotype[:, .!keep_index], :auto), :mean, :std, :min, :q25, :median, :q75, :max, :nuniqueall)
        describe_df.variable .= annotation.pheno_id[.!keep_index]
        CSV.write(elp_file, describe_df, delim="\t")
        phenotype = phenotype[:, keep_index]
        annotation = annotation[keep_index, :]
        n_samples_e, n_phenotypes = size(phenotype)
    end
    @runif check_chr_order(annotation.chrom) @error "The phenotype data is in unordered! Please sort it by chromosome and position."
    chroms = string.(unique(annotation.chrom))
    @runif USE_Float32 df_to_32bit!(annotation)
    @runif !isnothing(id_map) rematch_sample_id!(iid, id_map)
    if !isnothing(extract_samples)
        order_row_index = nothing
        if length(extract_samples) != length(iid)
            order_row_index = vmatch(iid, extract_samples,rm_nothing=true)
        else
            if any(extract_samples .!= iid)
                order_row_index = vmatch(iid, extract_samples,rm_nothing=true)
            end
        end
        if !isnothing(order_row_index)
            phenotype = phenotype[order_row_index,:]
            iid = iid[order_row_index]
        end
        n_samples_e, n_phenotypes = size(phenotype)
    end
    if _args_normalized_phenotype
        inverse_normal_transform!(phenotype)
    end
    annotation.index = range(1, size(annotation, 1))    
    Phenotype{eltype(phenotype)}(phenotype, annotation, iid, chroms, n_samples_e, n_phenotypes, phenotype_pairs)
end
function Phenotype_annot(pheno_annot_file::String; chrom::Vector{String}=String[], extract_pheno_file::Union{String,Nothing}=nothing, exclude_pheno_file::Union{String,Nothing}=nothing)
    annotation = nothing
    if isfile(pheno_annot_file)
        annotation = CSV.File(pheno_annot_file, header=true, buffer_in_memory=true, types=Dict(1 => String, 2 => Int64, 3 => Int64, 4 => String), select=[1, 2, 3, 4]) |> DataFrame
        rename!(annotation, ["chrom", "start", "end", "pheno_id"])
    else
        error("Annotation file does not exist, please check!")
    end
    if length(chrom) > 0
        annotation = annotation[annotation.chrom.∈(chrom,), :]
    end
    phenotype_pairs = nothing
    if !isnothing(extract_pheno_file)
        extract_pheno_df = CSV.File(extract_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(extract_pheno_df, 2) == 2
            rename!(extract_pheno_df, [:phe1, :phe2])
            extract_pheno_list = unique([extract_pheno_df.phe1; extract_pheno_df.phe2])
            phenotype_pairs = extract_pheno_df
        elseif size(extract_pheno_df, 2) == 1
            rename!(extract_pheno_df, [:phe1])
            extract_pheno_list = unique(extract_pheno_df.phe1)
        else
            error("")
        end
        extract_index = annotation.pheno_id .∈ (extract_pheno_list,)
        annotation = annotation[extract_index, :]
    end
    if !isnothing(exclude_pheno_file)
        exclude_pheno_df = CSV.File(exclude_pheno_file, header=false, buffer_in_memory=true, types=String) |> DataFrame
        if size(exclude_pheno_df, 2) == 1
            rename!(exclude_pheno_df, [:phe1])
            extract_pheno_list = unique(exclude_pheno_df.phe1)
        else
            error("")
        end
        exclude_index = .!(annotation.pheno_id .∈ (extract_pheno_list,))
        annotation = annotation[exclude_index, :]
    end
    @runif check_chr_order(annotation.chrom) @error "The phenotype data is in unordered! Please sort it by chromosome and position."
    chroms = string.(unique(annotation.chrom))
    n_phenotypes = size(annotation, 1)
    @runif USE_Float32 df_to_32bit!(annotation)
    annotation.index = range(1, size(annotation, 1))    
    Phenotype{FloatT}(nothing, annotation, nothing, chroms, nothing, n_phenotypes, phenotype_pairs)
end
mutable struct Covariates{T<:AbstractFloat}
    X_MME::Matrix{T}
    iid::Union{Vector{String},Nothing}
    n_samples::Int
    n_covars::Int
end
function Covariates(covar_file::Union{String,Nothing}, n_samples::Union{Nothing,Int}=nothing; extract_covar_name::Vector{String}=String[], exclude_covar_name::Vector{String}=String[], dcovar_name::Vector{String}=String[], covar_transpose::Bool=false, rm_collinear_covar::Union{Nothing,T}=nothing, id_map::Union{Nothing, DataFrame}=nothing, extract_samples::Union{Nothing, Vector{String}}=nothing) where {T<:AbstractFloat}
    if !isnothing(covar_file)
        if !isfile(covar_file)
            error(string("Covariates file ", covar_file, " does not exist, please check!"))
        end
        covar = CSV.File(covar_file, transpose=!covar_transpose, buffer_in_memory=true, missingstring=["nan", "na", "NA", "NaN", "NAN", ""], types=Dict(1 => String), typemap=Dict(Int => FloatT)) |> DataFrame
        if covar_transpose
            col_iid = findfirst(.!isnothing.(match.(r"IID$", names(covar))))
            if isnothing(col_iid)
                error("The column 'IID' was not found in the covariates file!")
            end
            if col_iid != 1
                covar = covar[:,col_iid:end]
            end
            covar[!,1] = string.(covar[:,1])
        end
        if length(dcovar_name) == 0
            strcol = [covar[1,x] isa AbstractString for x in 2:size(covar,2)]
            @runif any(strcol) dcovar_name = names(covar[:,2:end])[strcol]
        end
        iid = string.(covar[:, 1])
        @runif !isnothing(id_map) rematch_sample_id!(iid, id_map)
        @runif USE_Float32 df_to_32bit!(covar)
        order_row_index = nothing
        @runif !isnothing(extract_samples) if length(extract_samples) != length(iid)
            order_row_index = vmatch(iid, extract_samples)
        else
            if any(extract_samples .!= iid)
                order_row_index = vmatch(iid, extract_samples)
            end
        end
        if !isnothing(order_row_index)
            if any(isnothing.(order_row_index))
                error(string(extract_samples[findall(isnothing.(order_row_index))], " were not found in the covariates file!"))
            end
            iid = iid[order_row_index]
            covar = covar[order_row_index,:]
        end
        n_samples = length(iid)
        if length(dcovar_name) > 0
            covar[!,dcovar_name] .= string.(covar[:,dcovar_name])
            dcovar = USE_Float32 ? FloatT.(create_design_matrix(Matrix(covar[:, dcovar_name]))) : create_design_matrix(Matrix(covar[:, dcovar_name]))
        else
            dcovar = nothing
        end
        qcovar_index = findall(nonmissingtype.(eltype.(eachcol(covar))) .<: AbstractFloat)
        for qcovar_i in qcovar_index
            covar[ismissing.(covar[:,qcovar_i]),qcovar_i] .= Statistics.mean(skipmissing(covar[:,qcovar_i]))
        end
        if length(extract_covar_name) > 0
            keep_covar_colnames = [names(covar)[1]]
            append!(keep_covar_colnames, names(covar)[names(covar).∈(extract_covar_name,)])
            covar = covar[:, keep_covar_colnames]
        end
        if length(exclude_covar_name) > 0
            keep_covar_colnames = names(covar)[.!(names(covar) .∈ (exclude_covar_name,))]
            covar = covar[:, keep_covar_colnames]
        end
        if !isnothing(dcovar_name)
            qcovar_colnames = names(covar)[.!(names(covar) .∈ (dcovar_name,))][2:end]
        else
            qcovar_colnames = names(covar)[2:end]
        end
        if length(qcovar_colnames) > 0
            covar_var = vec(var(Matrix(covar[:, qcovar_colnames]), dims=1))
            if any(covar_var .< 1e-4)
                println(qcovar_colnames)
                println(covar_var)
                println_to_file(string("The variance of the covariate(s) ", qcovar_colnames[findall(covar_var .< 1e-4)], " < 1e-4, suggesting to exclude it/them using --exclude-covar-name."), log_file) 
            end
            if any(covar_var .> 1)
                println_to_file(string(" [INFO] The covariate(s) ", qcovar_colnames[findall(covar_var .> 1)], " with variance > 1, that were normalized to unit variance."), log_file)
                covar_std = std(Matrix(covar[:, qcovar_colnames[findall(covar_var .> 1)]]), dims=1)
                covar[:, qcovar_colnames[findall(covar_var .> 1)]] ./= covar_std
            end
        end
        if isnothing(dcovar)
            X_MME = covar[:, 2:end] |> Matrix{FloatT}
        else
            if length(qcovar_colnames) > 0 
                X_MME = [dcovar Matrix(covar[:, qcovar_colnames])] |> Matrix{FloatT}
            else
                X_MME = dcovar |> Matrix{FloatT}
            end
        end
        X_MME .-= Statistics.mean(X_MME, dims=1)
        if _args_normalized_covar
            @views inverse_normal_transform!(X_MME)
        end
        if !isnothing(rm_collinear_covar) & (size(X_MME,2) > 1)
            reduced_X_MME = remove_collinear_columns(X_MME, cor_cutoff=rm_collinear_covar)
            if size(reduced_X_MME, 2) != size(X_MME, 2)
                println_to_file(string("    * Eliminated ", size(X_MME, 2) - size(reduced_X_MME, 2), " collinear covariates with absolute correlation > ", rm_collinear_covar, "."), log_file)
                X_MME = reduced_X_MME
            end
        end
        if !_args_no_intercept
            X_MME = [ones(FloatT, n_samples, 1) X_MME]
        end
    else
        if isnothing(n_samples)
            error("Please provide sample size for Covariates")
        end
        iid = nothing
        X_MME = ones(FloatT, n_samples, 1)
    end
    n_covars = size(X_MME, 2)
    Covariates{eltype(X_MME)}(X_MME, iid, n_samples, n_covars)
end
function Covariates(_X_MME::Matrix{T2}, genotype::Union{Matrix{T1},Nothing}, af::Union{Vector{T2},Nothing}, phenotype::Union{Matrix{T2},Nothing}, n_geno_pc::Union{Nothing,Int}, n_phen_pc::Union{Nothing,Int}, dprop::Vector{Float64}, iid::Union{Vector,Nothing}; rm_collinear_covar::Union{Nothing,T3}=nothing, PCAForQTL::Union{Nothing,String}=nothing) where {T1<:Integer,T2<:AbstractFloat,T3<:AbstractFloat}
    n_samples = size(phenotype, 1)
    X_MME = copy(_X_MME)
    _X_c = size(X_MME, 2) - 1
    pc_covar_names = []
    if _X_c > 0
        append!(pc_covar_names, "c" .* string.(Vector(1:_X_c)))
    end
    if length(dprop) == 0
        dprop_gpc = nothing
        dprop_ppc = nothing
    elseif length(dprop) == 1
        dprop_gpc = dprop[1]
        dprop_ppc = dprop[1]
    elseif length(dprop) == 2
        dprop_gpc = dprop[1]
        dprop_ppc = dprop[2]
    else
        error("Incorrect input for `--dprop-pc-covar`.")
    end
    if n_geno_pc > 0
        n_samples, n_snps = size(genotype)
        n_snps = size(genotype, 2)
        if n_snps < 1000000
            pca_snp_index = 1:n_snps |> Vector
        else
            pca_snp_index = floor.(Int, Vector(range(1, n_snps, length=1000000)))
        end
        genotype_adj = genotype[:, pca_snp_index] .- (2 * af[pca_snp_index])'
        geno_pc = fit(PCA, genotype_adj, pratio=1)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".geno.eigvals")), DataFrame(pc=1:size(geno_pc, 2), eigenval=eigvals(geno_pc), varprop=geno_pc.prinvars ./ geno_pc.tprinvar), delim="\t")
        @runif n_geno_pc == n_samples if !isnothing(dprop_gpc)
            prop_explained = eigvals(geno_pc) ./ sum(eigvals(geno_pc))
            d_prop_explained_1 = prop_explained[1:end-1] - prop_explained[2:end]
            d_prop_explained_2 = prop_explained[1:end-2] - prop_explained[3:end]
            d_prop_explained = (d_prop_explained_1[1:end-1] .> dprop_gpc) .| (d_prop_explained_2 .> 2 * dprop_gpc)
            n_geno_pc = minimum([n_geno_pc, count_max_consecutive_true(d_prop_explained)])
        elseif !isnothing(PCAForQTL)
            @info "Performing PCAForQTL for genotype data ..."
            if PCAForQTL == "be"
                n_geno_pc, _ = runBE(genotype_adj, geno_pc; B=20, alpha=0.05, verbose=false)
            elseif PCAForQTL == "elbow"
                n_geno_pc, _ = runElbow(geno_pc)
            else
                error("Incorrect method for PCAForQTL!")
            end
            @info "Done."
        end
        println_to_file(string("    * ", n_geno_pc, " genotype PCs were included as covariates."), log_file)
        append!(pc_covar_names, "geno_pc" .* string.(Vector(1:n_geno_pc)))
        X_MME = [X_MME geno_pc.proj[:, 1:n_geno_pc]]
    end
    if n_phen_pc > 0
        phen_pc = fit(PCA, phenotype, pratio=1)
        CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".phen.eigvals")), DataFrame(pc=1:size(phen_pc, 2), eigenval=eigvals(phen_pc), varprop=phen_pc.prinvars ./ phen_pc.tprinvar), delim="\t")
        @runif n_phen_pc == n_samples if !isnothing(dprop_ppc)
            prop_explained = eigvals(phen_pc) ./ sum(eigvals(phen_pc))
            d_prop_explained_1 = prop_explained[1:end-1] - prop_explained[2:end]
            d_prop_explained_2 = prop_explained[1:end-2] - prop_explained[3:end]
            d_prop_explained = (d_prop_explained_1[1:end-1] .> dprop_ppc) .| (d_prop_explained_2 .> 2 * dprop_ppc)
            n_phen_pc = minimum([n_phen_pc, count_max_consecutive_true(d_prop_explained)])
        elseif !isnothing(PCAForQTL)
            @info "Performing PCAForQTL for phenotype data ..."
            if PCAForQTL == "be"
                n_phen_pc, _ = runBE(phenotype, phen_pc; B=20, alpha=0.05, verbose=false)
            elseif PCAForQTL == "elbow"
                n_phen_pc, _ = runElbow(phen_pc)
            else
                error("Incorrect method for PCAForQTL!")
            end
            @info "Done."
        end
        println_to_file(string("    * ", n_phen_pc, " phenotype PCs were included as covariates."), log_file)
        append!(pc_covar_names, "phen_pc" .* string.(Vector(1:n_phen_pc)))
        X_MME = [X_MME phen_pc.proj[:,1:n_phen_pc]]
    end
    if !isnothing(rm_collinear_covar)
        X_MME_keep_index = remove_collinear_columns(X_MME, cor_cutoff=rm_collinear_covar,return_index=true)
        reduced_X_MME = X_MME[:,X_MME_keep_index]
        pc_covar_names = pc_covar_names[X_MME_keep_index[2:end] .- 1]
        if size(reduced_X_MME, 2) != size(X_MME, 2)
            println_to_file(string("    * Eliminated ", size(X_MME, 2) - size(reduced_X_MME, 2), " collinear covariates with absolute correlation > ", rm_collinear_covar, "."), log_file)
            X_MME = reduced_X_MME
        end
    end
    n_covars = size(X_MME,2)
    println_to_file(string("\t[INFO] Automatically generated covariates were stored in ", joinpath(_args_output_dir, string(_args_out_prefix, ".auto_covar"))), log_file)
    CSV.write(joinpath(_args_output_dir, string(_args_out_prefix, ".auto_covar")), hcat(DataFrame(CovarName=pc_covar_names), DataFrame(X_MME[:,2:end]', iid)), delim="\t")
    Covariates{eltype(X_MME)}(X_MME, iid, n_samples, n_covars)
end
mutable struct InteractionTerm{T}
    ITERM::Union{Nothing,AbstractArray{T}}
    iid::Union{Vector{String},Nothing}
    n_iterms::Int
    n_samples::Int
end
function InteractionTerm(interaction_file::Union{Nothing,String}=nothing, id_map::Union{Nothing,DataFrame}=nothing; extract_samples::Union{Nothing,Vector{String}}=nothing)
    if !isnothing(interaction_file)
        if !isfile(interaction_file)
            error("Interaction file does not exist, please check!")
        end
        interaction_term = CSV.File(interaction_file, header=false, buffer_in_memory=true, types=Dict(1 => String), typemap=Dict(Int => FloatT)) |> DataFrame
        if interaction_term[1, 2] isa AbstractString
            interaction_term = CSV.File(interaction_file, header=true, buffer_in_memory=true, types=Dict(1 => String), typemap=Dict(Int => FloatT)) |> DataFrame
        end
        @runif USE_Float32 df_to_32bit!(interaction_term)
        iid = string.(interaction_term[:, 1])
        @runif !isnothing(id_map) rematch_sample_id!(iid, id_map)
        if !isnothing(_args_interaction_name)
            ITERM = interaction_term[:, _args_interaction_name]
        else
            if _args_experimental
                ITERM = interaction_term[:, 2:end] |> Matrix{FloatT} 
            else
                ITERM = interaction_term[:, 2:2] |> Matrix{FloatT} 
            end
        end
        order_row_index = nothing
        @runif !isnothing(extract_samples) if length(extract_samples) != length(iid)
            order_row_index = vmatch(iid, extract_samples)
        else
            if any(extract_samples .!= iid)
                order_row_index = vmatch(iid, extract_samples)
            end
        end
        if !isnothing(order_row_index)
            iid = iid[order_row_index]
            ITERM = ITERM[order_row_index, :]
        end
        if _args_normalized_interaction
            ITERM .= inverse_normal_transform(ITERM) |> Matrix{FloatT}
        end
        var_INT = [var(ITERM[:, x]) for x in axes(ITERM, 2)]
        if any(var_INT == 0)
            error("The variance of interaction term ", findall(var_INT == 0), " equal to zero!")
        end
        n_iterms = size(ITERM, 2)
        n_samples = length(iid)
        if n_iterms > 1
            if rank(ITERM) < n_iterms
                vif_df = calculate_vif(ITERM; add_intercept=true)
                vif_file = replace(log_file, r".log$" => ".iterm.vif")
                CSV.write(vif_file, vif_df, delim="\t")
                error(string("\tThe interaction term matrix is not full rank (", rank(ITERM), "), please remove the redundant items with high VIF (", vif_file, ")."))
            end
            if rank(hcat(ones(n_samples), ITERM)) <= n_iterms
                @warn string("\tTo avoid multicollinearity, the genetic main effect will be excluded as it is entirely explained by the interactions.")
                global gxe_omit_main_eff = true 
            end
        end
        InteractionTerm{eltype(ITERM)}(ITERM, iid, n_iterms, n_samples)
    else
        ITERM = nothing
        iid = nothing
        n_iterms = 0
        n_samples = 0
        InteractionTerm{Nothing}(ITERM, iid, n_iterms, n_samples)
    end
end
mutable struct PriorHerB
    h2::Union{Nothing,DataFrame}
    B::Union{Nothing,DataFrame}
end
function PriorHerB(prior_herb_file::Union{Nothing,String,Vector{String}}=nothing)
    if prior_herb_file isa String
        her_output = prior_herb_file
        if isfile(her_output)
            prior_h2 = CSV.File(her_output, header=true, buffer_in_memory=true) |> DataFrame
            prior_h2_B_file = replace(her_output, ".txt" => ".B")
            prior_B = CSV.File(prior_h2_B_file, header=true, buffer_in_memory=true) |> DataFrame
        else
            error(string("No found ", her_output))
        end
        @runif USE_Float32 df_to_32bit!(prior_h2)
        @runif USE_Float32 df_to_32bit!(prior_B)
    elseif prior_herb_file isa Vector{String}
        her_output_multi = prior_herb_file
        prior_h2 = DataFrame()
        prior_B = DataFrame()
        for fn in her_output_multi
            _df_hers = CSV.File(fn, header=true, buffer_in_memory=true) |> DataFrame
            prior_h2_B_file = replace(fn, ".txt" => ".B")
            _prior_B = CSV.File(prior_h2_B_file, header=true, buffer_in_memory=true) |> DataFrame
            append!(prior_h2, _df_hers)
            append!(prior_B, _prior_B)
        end
        @runif USE_Float32 df_to_32bit!(prior_h2)
        @runif USE_Float32 df_to_32bit!(prior_B)
    else
        prior_h2 = prior_B = nothing
    end
    PriorHerB(prior_h2, prior_B)
end
struct Kinship{T<:AbstractFloat}
    GRM::Union{Matrix{T},Vector{Matrix{T}}}
    EA::Union{Nothing,Eigen{T,T,Matrix{T},Vector{T}}}
    domGRM::Union{Nothing,Matrix{T}}
    domEA::Union{Nothing,Eigen{T,T,Matrix{T},Vector{T}}}
    n_grm::Int
end
function Kinship(genotype::AbstractMatrix{T1}, adjusted::Bool; afs::Vector{T2}, byrow::Bool=false, return_sparse::Union{Nothing,T3}=nothing, dom_genotype::Union{Nothing,Matrix{UInt8}}=nothing,grm_alpha::Int=0) where {T1<:Real,T2<:AbstractFloat,T3<:AbstractFloat}
    n_snps = size(genotype, byrow ? 2 : 1)
    grm_snps = isnothing(_args_grm_snps) ? n_snps : _args_grm_snps
    @runif n_snps > grm_snps begin
        grm_snp_index = floor.(Int, Vector(range(1, n_snps, length=grm_snps)))
        println_to_file(string("\n\t* ", grm_snps, " variants used for building GRM."), log_file)
    end
    if n_snps <= grm_snps
        GRM = getG(genotype; alpha=grm_alpha, afs=afs, return_eltype=FloatT)
    else
        GRM = getG(genotype[:, grm_snp_index]; alpha=grm_alpha, afs=afs[grm_snp_index], return_eltype=FloatT)
    end
    if _args_debug
        println_to_file("\n\t* Description for diagonal elements of GRM:", log_file)
        println_to_file(string(summarystats(GRM[diagind(GRM)])), log_file)
        println_to_file("\n\t* Description for lower triangular elements of GRM:", log_file)
        println_to_file(string(summarystats(GRM[lower_tri(size(GRM, 1))])), log_file)
    end
    GC.gc()
    if !isnothing(return_sparse)
        GRM[GRM.<return_sparse] .= 0 
    end
    EA = eigen(diagm(repeat([FloatT(1.0)],size(GRM,1))))
    println_to_file("Eigen decomposition ...\t", log_file)
    @timeit to "svd GRM" if _args_eigen_algo == "svd"
        eigen_vec, eigen_val = svd(GRM)
    else 
        eigen_vec, eigen_val = LinearAlgebra.LAPACK.gesdd!('A', copy(GRM))
    end
    EA.vectors .= eigen_vec
    EA.values .= eigen_val
    if !isnothing(dom_genotype)
        if n_snps <= grm_snps
            domGRM = getG_dominance(dom_genotype, code_type="dominance", afs=Float64.(afs), byrow=byrow, return_eltype=FloatT)
        else
            domGRM = getG_dominance(byrow ? dom_genotype[:, grm_snp_index] : dom_genotype[grm_snp_index, :], code_type="dominance", afs=Float64.(afs[grm_snp_index]), byrow=byrow, return_eltype=FloatT)
        end
        if _args_debug
            println_to_file("\n\t* Description for diagonal elements of domGRM:", log_file)
            println_to_file(string(summarystats(domGRM[diagind(domGRM)])), log_file)
            println_to_file("\n\t* Description for lower triangular elements of domGRM:", log_file)
            println_to_file(string(summarystats(domGRM[lower_tri(size(domGRM, 1))])), log_file)
        end
        if !isnothing(return_sparse)
            domGRM[domGRM.<return_sparse] .= 0 
        end
        GC.gc()
        domEA = deepcopy(EA)
        if _args_eigen_algo == "svd"
            eigen_vec, eigen_val = svd(domGRM)
        else
            eigen_vec, eigen_val = LinearAlgebra.LAPACK.gesdd!('A', copy(domGRM))
        end
        domEA.vectors .= eigen_vec
        domEA.values .= eigen_val
    else
        domGRM = nothing
        domEA = nothing
    end
    n_grm = !isnothing(GRM) + !isnothing(domGRM)
    Kinship{eltype(GRM)}(GRM, EA, domGRM, domEA, n_grm)
end
function Kinship(GRM::Matrix{T1}; domGRM::Union{Nothing,Matrix{T1}}=nothing, return_sparse::Union{Nothing,T2}=nothing) where {T1,T2<:AbstractFloat}
    if !isnothing(return_sparse)
        GRM[GRM.<return_sparse] .= 0 
    end
    EA = eigen(diagm(repeat([FloatT(1.0)],size(GRM,1))))
    println_to_file("Eigen decomposition ...\t", log_file)
    @timeit to "svd GRM" if _args_eigen_algo == "svd"
        eigen_vec, eigen_val = svd(GRM)
    else
        eigen_vec, eigen_val = LinearAlgebra.LAPACK.gesdd!('A', copy(GRM))
    end
    EA.vectors .= eigen_vec
    EA.values .= eigen_val
    if !isnothing(domGRM)
        if !isnothing(return_sparse)
            domGRM[domGRM.<return_sparse] .= 0 
        end
        domEA = deepcopy(EA)
        if _args_eigen_algo == "svd"
            eigen_vec, eigen_val = svd(domGRM)
        else
            eigen_vec, eigen_val = LinearAlgebra.LAPACK.gesdd!('A', copy(domGRM))
        end
        domEA.vectors .= eigen_vec
        domEA.values .= eigen_val
    else
        domEA = nothing
    end
    n_grm = !isnothing(GRM) + !isnothing(domGRM)
    Kinship{eltype(GRM)}(GRM, EA, domGRM, domEA, n_grm)
end
function create_OUTDIR_and_LOG(output_dir::String, run_mode::Union{String, Nothing}, chrom::Vector{String}, out_prefix::String)
    mkpath(output_dir)
    if run_mode == "her_est"
        mediate_file_name = "her_est"
    elseif run_mode == "cis"
        mediate_file_name = "cis_qtl"
    elseif run_mode == "cis_independent"
        mediate_file_name = "cis_independent"
    elseif run_mode == "cis_interaction"
        mediate_file_name = "cis_interaction"
    elseif run_mode == "cis_mt"
        mediate_file_name = "cis_mt"
    elseif run_mode == "trans"
        mediate_file_name = "trans_qtl"
    elseif run_mode == "gwas"
        mediate_file_name = "gwas_qtl"
    elseif run_mode == "plot"
        mediate_file_name = "plot"
    elseif run_mode == "birg"
        mediate_file_name = "birg"
    elseif run_mode == "xbirg"
        mediate_file_name = "xbirg"
    elseif run_mode == "simu"
        mediate_file_name = "simu"
    elseif run_mode == "enrich"
        mediate_file_name = "enrich"
    elseif run_mode == "uint_test"
        mediate_file_name = "unit_test"
    else
        mediate_file_name = "omiga"
        @info "No run mode provied, the program will stop after initialization."
    end
    if !isnothing(match(r"/", out_prefix))
        error("The option --prefix cannot contain the special character '/'!")
    end
    if (length(chrom) == 0) || (length(chrom) > 1)
        log_file = string(out_prefix, ".", mediate_file_name, ".log")
    else
        log_file = string(out_prefix, ".", mediate_file_name, ".", chrom[1], ".log")
    end
    log_file = joinpath(output_dir, log_file)
    if isfile(log_file)
        rm(log_file, force=true)
    end
    return log_file
end
function prepare_Data()
    if _args_low_mem
        if _args_run_mode != "cis"
            error("Current version only support 'cis' mode when using '--chunk-size'.")
        end
    end
    use_dominance = false
    if _args_run_mode in ["simu"]
        use_dominance = true
    elseif _args_run_mode in ["her_est"]
        if !isnothing(_args_h2_model)
            use_dominance = !isnothing(match(r"d|D", _args_h2_model))
        end
    else
        use_dominance = !isnothing(match(r"d|D", _args_qtl_map_model))
    end
    @timeit to "Load data from disk" begin
        if length(_args_phenotype_files) > 0
            phenotype_bed_file = _args_phenotype_files[1]
        else
            phenotype_bed_file = nothing
        end
        if length(_args_covar_files) > 0
            covar_file = _args_covar_files[1]
        else
            covar_file = nothing
        end
        if length(_args_extract_pheno_files) > 0
            extract_pheno_file = _args_extract_pheno_files[1]
        else
            extract_pheno_file = nothing
        end
        if length(_args_extract_pheno_files) == 0
            extract_pheno_file = nothing
            extract_pheno_file_2 = nothing
        end
        extract_covar_name = _args_extract_covar_name
        extract_covar_name_2 = _args_extract_covar_name
        exclude_covar_name = _args_exclude_covar_name
        exclude_covar_name_2 = _args_exclude_covar_name
        dcovar_name = _args_dcovar_name
        dcovar_name_2 = _args_dcovar_name
        if _args_run_mode == "xbirg"
            if length(_args_extract_pheno_files) == 2
                extract_pheno_file = _args_extract_pheno_files[1]
                extract_pheno_file_2 = _args_extract_pheno_files[2]
            end
            if length(_args_extract_covar_name) > 0
                extract_covar_name, extract_covar_name_2 = split.(split(_args_extract_covar_name,"|"),",")
                if extract_covar_name == [""]
                    extract_covar_name = String[]
                end
                if extract_covar_name_2 == [""]
                    extract_covar_name_2 = String[]
                end
            end
            if length(_args_exclude_covar_name) > 0
                exclude_covar_name, exclude_covar_name_2 = split.(split(_args_exclude_covar_name,"|"),",")
                if exclude_covar_name == [""]
                    exclude_covar_name = String[]
                end
                if exclude_covar_name_2 == [""]
                    exclude_covar_name_2 = String[]
                end
            end
            if length(_args_dcovar_name) > 0
                dcovar_name, dcovar_name_2 = split.(split(_args_dcovar_name,"|"),",")
                if dcovar_name == [""]
                    dcovar_name = String[]
                end
                if dcovar_name_2 == [""]
                    dcovar_name_2 = String[]
                end
            end
        end
        id_map = nothing
        if !isnothing(_args_id_map_file)
            id_map = CSV.File(_args_id_map_file) |> DataFrame
        end
        geno_sample_ids = nothing
        if !isnothing(_args_geno_file_prefix) 
            if isfile(string(_args_geno_file_prefix,".fam"))
                fam = CSV.read(string(_args_geno_file_prefix, ".fam"), DataFrame, header=false, buffer_in_memory=true, types=Dict(1 => String, 2 => String))
                geno_sample_ids = string.(fam[:, 2])::Vector{String}
                @runif !isnothing(id_map) rematch_sample_id!(geno_sample_ids, id_map)
            else
                error("Not found '", string(_args_geno_file_prefix,".fam'"))
            end
        end
        loaded_pheno = true
        if !isnothing(phenotype_bed_file)
            pheno_file_format = endswith(phenotype_bed_file,r".bed|.bed.gz") ? "bed" : "plink"
        else
            pheno_file_format = "bed"
        end
        if (isnothing(phenotype_bed_file)) & (!isnothing(_args_pheno_annot_file))
            if _args_run_mode in ["simu"]
                pheno_file_format = "only_annot"
                loaded_pheno = false
            else
                error("The option --phenotype must be specified!")
            end
        end
        if isnothing(_args_run_mode) & isnothing(phenotype_bed_file)
            loaded_pheno = false
            PHENO = nothing
            COVAR = nothing
        end
        @runif loaded_pheno println_to_file("Loading phenotype file ...", log_file)
        @runif !loaded_pheno println_to_file("Loading phenotype annotation file ...\t", log_file)
        @runif loaded_pheno elapsed_time = @elapsed if (pheno_file_format == "bed") & (isnothing(_args_pheno_annot_file))
            PHENO = Phenotype(phenotype_bed_file; rm_pheno_threshold=_args_rm_pheno_threshold, chrom=_args_chrom, extract_pheno_file=extract_pheno_file, extract_pheno_name=_args_extract_pheno_name, exclude_pheno_file=_args_exclude_pheno, exclude_pheno_name=_args_exclude_pheno_name, id_map=id_map, extract_samples=geno_sample_ids)
        elseif pheno_file_format == "plink"
            PHENO = Phenotype(phenotype_bed_file, _args_pheno_annot_file; rm_pheno_threshold=_args_rm_pheno_threshold, chrom=_args_chrom, extract_pheno_file=extract_pheno_file, extract_pheno_name=_args_extract_pheno_name, exclude_pheno_file=_args_exclude_pheno, exclude_pheno_name=_args_exclude_pheno_name, id_map=id_map,extract_samples=geno_sample_ids)
        elseif pheno_file_format == "only_annot"
            PHENO = Phenotype_annot(_args_pheno_annot_file; chrom=_args_chrom, extract_pheno_file=extract_pheno_file, exclude_pheno_file=_args_exclude_pheno)
        else
            error("Please check if the option --pheno-format is correct.")
        end
        if !isnothing(_args_pheno_group_file)
            println_to_file("Loading group annotation file ...\t", log_file)
            phenotype_group = CSV.File(_args_pheno_group_file, header=false) |> DataFrame
            rename!(phenotype_group, ["pheno_id", "group_id"])
            if size(phenotype_group, 1) != PHENO.n_phenotypes
                error(string("DimensionMismatch: unmatched length between phenotypes and group annotation!"))
            end
            if sum(PHENO.annotation.pheno_id .!= phenotype_group.pheno_id) > 0
                error("IDMismatch: mismatched phenotype id between phenotype and group files!")
            end
            rep_counts, rep_first_indices = count_consecutive(phenotype_group.group_id)
            _, rep_last_indices = count_consecutive(phenotype_group.group_id; type="last")
            if length(rep_counts) != length(unique(phenotype_group.group_id))
                count_rle_values = countmap(rle(phenotype_group.group_id)[1])
                frequent_elements = Array{eltype(count_rle_values.keys)}(undef, 0)
                for (element, count) in count_rle_values
                    if count > 1
                        push!(frequent_elements, element)  
                    end
                end
                error(string("Inconsecutive group ID: ", frequent_elements, " in phenotype group file, please check!"))
            end
            phenotype_group.first_appear .= false
            phenotype_group.first_appear[rep_first_indices] .= true
            phenotype_group.last_appear .= false
            phenotype_group.last_appear[rep_last_indices] .= true
            PHENO.annotation.group_id = phenotype_group.group_id
            PHENO.annotation.group_id_first = phenotype_group.first_appear
            PHENO.annotation.group_id_last = phenotype_group.last_appear
        end
        @runif loaded_pheno println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
        @runif loaded_pheno @runif PHENO.n_samples < 5 error(string("Only ", PHENO.n_samples, " samples were loaded! Please check whether the sample IDs are matched between genotype and phenotype data. The `--id-map` option can be applied to provide ID conversion."))
        @runif loaded_pheno println_to_file(string("\tLoaded (filtered) ", PHENO.n_samples, " samples, ", PHENO.n_phenotypes, " phenotypes"), log_file)
        if loaded_pheno
            phenotype_sample_id = PHENO.iid
        else
            phenotype_sample_id = nothing
        end
        if length(_args_phenotype_files) > 1
            phenotype_file_2 = _args_phenotype_files[2]
            elapsed_time = @elapsed if (pheno_file_format == "bed") & (isnothing(_args_pheno_annot_file))
                PHENO_2 = Phenotype(phenotype_file_2; rm_pheno_threshold=_args_rm_pheno_threshold, chrom=_args_chrom, extract_pheno_file=extract_pheno_file_2, extract_pheno_name=_args_extract_pheno_name, exclude_pheno_file=_args_exclude_pheno, exclude_pheno_name=_args_exclude_pheno_name, id_map=id_map, extract_samples=phenotype_sample_id)
            elseif pheno_file_format == "plink"
                PHENO_2 = Phenotype(phenotype_file_2, _args_pheno_annot_file; rm_pheno_threshold=_args_rm_pheno_threshold, chrom=_args_chrom, extract_pheno_file=extract_pheno_file_2, extract_pheno_name=_args_extract_pheno_name, exclude_pheno_file=_args_exclude_pheno, exclude_pheno_name=_args_exclude_pheno_name, id_map=id_map, extract_samples=phenotype_sample_id)
            elseif pheno_file_format == "only_annot"
                PHENO_2 = Phenotype_annot(_args_pheno_annot_file; chrom=_args_chrom, extract_pheno_file=extract_pheno_file_2, exclude_pheno_file=_args_exclude_pheno)
            else
                error("Please check if the option --pheno-format is correct.")
            end
        end
        loaded_inte = false
        if _args_run_mode in ["cis_interaction", "simu"]
            if _args_run_mode == "cis_interaction"
                if isnothing(_args_interaction_file)
                    error("The option --interaction must be provided for `--mode cis_interaction`!")
                end
                loaded_inte = true
            end
            println_to_file("Loading interaction file ...", log_file)
            elapsed_time = @elapsed ITERM = InteractionTerm(_args_interaction_file, id_map; extract_samples=phenotype_sample_id)
            println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
            ITERM_vec = ITERM.ITERM
        else
            ITERM = nothing
            ITERM_vec = nothing
        end
        @runif loaded_pheno if !isnothing(covar_file) 
            println_to_file("Loading covariates file ...", log_file)
            elapsed_time = @elapsed COVAR = Covariates(covar_file; extract_covar_name=extract_covar_name, exclude_covar_name=exclude_covar_name, dcovar_name=dcovar_name, covar_transpose=_args_covar_transpose, rm_collinear_covar=_args_rm_collinear_covar, id_map=id_map, extract_samples=phenotype_sample_id)
            println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
            println_to_file(string("\tLoaded ", COVAR.n_samples, " samples, ", COVAR.n_covars - ifelse(_args_no_intercept, 0, 1), " covariates"), log_file)
        end
        loaded_geno = false
        if !isnothing(_args_geno_file_prefix) 
            println_to_file("Loading genotype file ...", log_file)
            elapsed_time = @elapsed GENO = Genotype(_args_geno_file_prefix; mac_threshold=_args_mac_threshold, maf_threshold=FloatT(_args_maf_threshold), maf_threshold_interaction=FloatT(_args_maf_threshold_interaction), het_threshold_interaction=FloatT(_args_het_threshold_interaction), het_rate_threshold=FloatT(_args_het_rate_threshold), byrow=is_sample_byrow, ITERM=ITERM_vec, keep_raw_geno=false, extract_samples=phenotype_sample_id, id_map=id_map, low_mem=_args_low_mem, use_dominance=use_dominance)
            println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
            println_to_file(string("\tLoaded (filtered) ", GENO.n_samples, " samples, ", GENO.n_snps, " variants"), log_file)
            loaded_geno = true
        else
            GENO = nothing
        end
        @runif loaded_pheno if isnothing(covar_file) 
            if !isnothing(PHENO.n_samples)
                COVAR = Covariates(covar_file, PHENO.n_samples)
            else
                COVAR = Covariates(covar_file, GENO.n_samples)
            end
        end
        if length(_args_covar_files) > 1
            covar_file_2 = _args_covar_files[2]
            if !isnothing(covar_file_2)
                println_to_file("Loading covariates file ...\t", log_file)
                elapsed_time = @elapsed COVAR_2 = Covariates(covar_file_2; extract_covar_name=extract_covar_name_2, exclude_covar_name=exclude_covar_name_2, dcovar_name=dcovar_name_2, covar_transpose=_args_covar_transpose, rm_collinear_covar=_args_rm_collinear_covar, id_map=id_map, extract_samples=phenotype_sample_id)
                println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
                println_to_file(string("\tLoaded ", COVAR_2.n_samples, " samples, ", COVAR_2.n_covars - ifelse(_args_no_intercept, 0, 1), " covariates"), log_file)
            else
                if !isnothing(PHENO.n_samples)
                    COVAR_2 = Covariates(covar_file_2, PHENO.n_samples)
                else
                    COVAR_2 = Covariates(covar_file_2, GENO.n_samples)
                end
            end
        end
        if !_args_low_mem
            if !isnothing(_args_geno_pc_covar) || !isnothing(_args_phen_pc_covar) || length(_args_dprop_pc_covar) > 0 || !isnothing(_args_PCAForQTL) 
                n_geno_pc_covar = _args_geno_pc_covar
                n_phen_pc_covar = _args_phen_pc_covar
                if (length(_args_dprop_pc_covar) > 0) || !isnothing(_args_PCAForQTL)
                    @runif isnothing(n_geno_pc_covar) n_geno_pc_covar = GENO.n_samples
                    @runif isnothing(n_phen_pc_covar) n_phen_pc_covar = GENO.n_samples
                else
                    @runif isnothing(n_geno_pc_covar) n_geno_pc_covar = 0
                    @runif isnothing(n_phen_pc_covar) n_phen_pc_covar = 0
                end
                COVAR = Covariates(COVAR.X_MME, isnothing(GENO) ? nothing : GENO.genotype, isnothing(GENO) ? nothing : GENO.annotation.af, isnothing(PHENO) ? nothing : PHENO.phenotype, n_geno_pc_covar, n_phen_pc_covar, isnothing(_args_PCAForQTL) ? _args_dprop_pc_covar : Float64[], isnothing(COVAR.iid) ? PHENO.iid : COVAR.iid; rm_collinear_covar=_args_rm_collinear_covar, PCAForQTL=_args_PCAForQTL)
            end
        else
            @info "--dprop-pc-covar is disable when using --chunk-size."
        end
        @runif !isnothing(ITERM_vec) if !_args_no_add_interaction_to_covar
            COVAR.X_MME = remove_collinear_columns([COVAR.X_MME ITERM_vec], cor_cutoff=_args_rm_collinear_covar, return_index=false)
            COVAR.n_covars = size(COVAR.X_MME, 2)
        end
        @runif loaded_pheno if rank(COVAR.X_MME) < COVAR.n_covars
            @warn string("\tThe covariates matrix is not full rank (", rank(COVAR.X_MME), "), the redundant covariates with high VIF were removed.")
            @suppress_err COVAR.X_MME = remove_multicollinearity(COVAR.X_MME)
            COVAR.n_covars = size(COVAR.X_MME, 2)
        end
        @runif !loaded_inte println_to_file(string("    * ", COVAR.n_covars, " covariates", ifelse(_args_no_intercept, "", " (including intercept)"), " were used."), log_file)
        @runif loaded_inte println_to_file(string("    * ", COVAR.n_covars, " covariates", ifelse(_args_no_intercept, "", " (including intercept and/or interaction terms)"), " were used."), log_file)
        if length(_args_covar_files) > 1
            if rank(COVAR_2.X_MME) < COVAR_2.n_covars
                error(string("\tThe covariates matrix is not full rank (",rank(COVAR_2.X_MME),"), please check!"))
            end
        end
        loaded_domgeno = false
        @runif !loaded_geno DOM = nothing
        @runif loaded_geno if use_dominance
            println_to_file("Processing dominance data ...\t", log_file)
            elapsed_time = @elapsed DOM = Dominance(GENO.genotype, GENO.annotation, is_sample_byrow; geno_file_prefix=_args_geno_file_prefix)
            println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
            loaded_domgeno = true
        else
            DOM = nothing
        end
        if _args_run_mode != "her_est"
            if length(_args_her_output) == 1
                println_to_file("Loading prior HerB file ...", log_file)
                PRIOR_HerB = PriorHerB(_args_her_output[1])
            elseif length(_args_her_output) > 1
                println_to_file("Loading prior HerB files ...", log_file)
                PRIOR_HerB = PriorHerB(_args_her_output)
            else
                PRIOR_HerB = nothing
            end
        else
            PRIOR_HerB = nothing
        end
        println_to_file("Data have been loaded.", log_file)
    end
    @timeit to "Summary and Check data" begin
        @runif loaded_geno & loaded_pheno if sum(GENO.iid .!= PHENO.iid) > 0
            error(string("SampleMismatch: mismatched sample id ", findall(GENO.iid .!= PHENO.iid), " between genotype and phenotype!"))
        end
        @runif !isnothing(covar_file) if sum(COVAR.iid .!= PHENO.iid) > 0
            error(string("SampleMismatch: mismatched sample id ", findall(COVAR.iid .!= PHENO.iid), " between covariate and phenotype!"))
        end
        @runif loaded_inte if sum(ITERM.iid .!= PHENO.iid) > 0
            error(string("SampleMismatch: mismatched sample id ", findall(ITERM.iid .!= PHENO.iid), " between interaction term and phenotype!"))
        end
        @runif loaded_geno & loaded_pheno if GENO.n_samples != PHENO.n_samples
            error(string("DimensionMismatch: genotype has ", GENO.n_samples, " samples but phenotype has ", PHENO.n_samples, " samples!"))
        end
        @runif loaded_pheno if COVAR.n_samples != PHENO.n_samples
            error(string("DimensionMismatch: covariate has ", COVAR.n_samples, " samples but phenotype has ", PHENO.n_samples, " samples!"))
        end
        @runif loaded_pheno if length(unique(PHENO.annotation.pheno_id)) != PHENO.n_phenotypes
            error(string("Duplicate pheno_id exist!"))
        end
    end
    build_grm = true
    if (_args_run_mode in ["simu","plot"]) || _args_linear_model
        build_grm = false
    end
    if build_grm & _args_low_mem
        if length(_args_grm_prefix) == 0
            error("Please specify the '--grm' flag when using '--chunk-size'.")
        end
    end
    @timeit to "Build GRM" begin
        n_grm_specified = length(_args_grm_prefix)
        @runif build_grm println_to_file(string("Number of GRM specified: ", n_grm_specified), log_file)
        elapsed_time = @elapsed if (n_grm_specified > 0) & build_grm
            println_to_file("Loading GRM ...", log_file)
            GRM_list = [zeros(FloatT, PHENO.n_samples, PHENO.n_samples) for _ in 1:n_grm_specified] 
            for grm_i in 1:n_grm_specified
                xgrmmat, xgrmiid, xgrmfid = readGRMBin(_args_grm_prefix[grm_i]; type=Float32, return_type=FloatT) 
                @runif !isnothing(id_map) rematch_sample_id!(xgrmiid, id_map)
                xgrmmat = xgrmmat[vmatch(xgrmiid, PHENO.iid), vmatch(xgrmiid, PHENO.iid)]
                copy!(GRM_list[grm_i], xgrmmat)
            end
            if n_grm_specified == 1
                KIN = Kinship(GRM_list[1]; return_sparse=_args_sparse_grm)
            elseif (n_grm_specified == 2) & use_dominance
                KIN = Kinship(GRM_list[1]; domGRM=GRM_list[2], return_sparse=_args_sparse_grm)
            elseif n_grm_specified > 1
                KIN = Kinship(GRM_list, nothing, nothing, nothing, n_grm_specified)
            end 
        elseif loaded_geno & build_grm
            println_to_file("Building GRM ...", log_file)
            if isnothing(_args_regional_grm)
                if use_dominance
                    KIN = Kinship(GENO.genotype, false; afs=GENO.annotation.af, byrow=is_sample_byrow, dom_genotype=DOM.genotype, return_sparse=_args_sparse_grm, grm_alpha=_args_grm_alpha)
                    if _args_save_grm
                        writeGRMBin(KIN.GRM, GENO.iid, joinpath(_args_output_dir, string(_args_out_prefix, ".a")), Float32; fam_id=GENO.fid)
                        writeGRMBin(KIN.domGRM, GENO.iid, joinpath(_args_output_dir, string(_args_out_prefix, ".d")), Float32; fam_id=GENO.fid)
                    end
                else
                    KIN = Kinship(GENO.genotype, false; afs=GENO.annotation.af, byrow=is_sample_byrow, return_sparse=_args_sparse_grm, grm_alpha=_args_grm_alpha)
                    if _args_save_grm
                        writeGRMBin(KIN.GRM, GENO.iid, joinpath(_args_output_dir, string(_args_out_prefix, ".a")), Float32; fam_id=GENO.fid)
                    end
                end
            else
                variant_regional_map = CSV.File(_args_regional_grm; header=false, types=String) |> DataFrame
                rename!(variant_regional_map, ["variant_id", "region"])
                grm_regions = unique(variant_regional_map.region)
                n_grm_region = length(grm_regions)
                println_to_file(string("Number of regional specified for GRM construction: ", n_grm_region), log_file)
                GRM_list = [zeros(FloatT, PHENO.n_samples, PHENO.n_samples) for _ in 1:n_grm_region] 
                for grm_i in 1:n_grm_region
                    println_to_file(string("    Processing region: ", grm_i, "/", n_grm_region, " <", grm_regions[grm_i],">"), log_file)
                    regional_variants = variant_regional_map.region .== grm_regions[grm_i]
                    b_set = Set(variant_regional_map.variant_id[regional_variants])
                    in_index = [x in b_set for x in GENO.annotation.variant]
                    grm_snp_index = GENO.annotation.index[in_index]
                    GRM = getG(GENO.genotype[:, grm_snp_index]; alpha=_args_grm_alpha, afs=GENO.annotation.af[grm_snp_index], return_eltype=FloatT)
                    copy!(GRM_list[grm_i], GRM)
                    if _args_save_grm
                        writeGRMBin(GRM, GENO.iid, joinpath(_args_output_dir, string(_args_out_prefix, ".", grm_i)), Float32; fam_id=GENO.fid)
                    end
                end
                if n_grm_region == 1
                    KIN = Kinship(GRM_list[1]; return_sparse=_args_sparse_grm)
                elseif (n_grm_region == 2) & use_dominance
                    KIN = Kinship(GRM_list[1]; domGRM=GRM_list[2], return_sparse=_args_sparse_grm)
                elseif n_grm_region > 1
                    KIN = Kinship(GRM_list, nothing, nothing, nothing, n_grm_region)
                end 
            end
        elseif !build_grm
            KIN = nothing
        else
            println_to_file(string("\t* '--genotype' or '--grm' must be specified for analysis!"), log_file)
        end
        println_to_file(string("\tElapsed time: ", elapsed_time, " seconds"), log_file)
    end
    GC.gc()
    if _args_run_mode != "xbirg"
        return PHENO, GENO, KIN, COVAR, ITERM, PRIOR_HerB, DOM
    else
        return PHENO, PHENO_2, GENO, KIN, COVAR, COVAR_2
    end
end
