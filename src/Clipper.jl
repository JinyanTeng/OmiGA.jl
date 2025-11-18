function Clipper(score_exp, score_back; analysis = "enrichment", FDR = 0.05, procedure = nothing, contrast_score = nothing, n_permutation = nothing, seed = 12345)
  if analysis == "differential"
      if isnothing(contrast_score)
          contrast_score = "max"
      end
      if isnothing(procedure)
          procedure = "GZ"
      end
      re = clipper2sided(score_exp = score_exp, score_back = score_back, FDR = FDR,
                        nknockoff = n_permutation,
                        contrastScore_method = contrast_score, importanceScore_method = "diff",
                        FDR_control_method = procedure, ifpowerful = false, seed = seed)
      FDR_nodisc = [length(re_i["discovery"]) == 0 for re_i in re["results"]]
      if any(FDR_nodisc .& (contrast_score == "max"))
          @warn string("At FDR = ", join([FDR][FDR_nodisc], ", "), ", no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold.")
      end
  end
  if isnothing(size(score_exp, 1))
      score_exp = reshape(score_exp, :, 1)
  end
  if isnothing(size(score_back, 1))
      score_back = reshape(score_back, :, 1)
  end
  if analysis == "enrichment"
      if isnothing(procedure)
          procedure = "BC"
      end
      if size(score_exp, 2) != size(score_back, 2)
          procedure = "GZ"
      end
      if isnothing(contrast_score)
          if procedure == "BC"
              contrast_score = "diff"
          elseif procedure == "GZ"
              contrast_score = "max"
          end
      end
      if procedure == "BC"
          re = clipper1sided(score_exp = score_exp, score_back = score_back, FDR = FDR,
                          nknockoff = n_permutation,
                          importanceScore_method = contrast_score,
                          FDR_control_method = procedure, ifpowerful = false, seed = seed)
      elseif procedure == "GZ"
          re = clipper1sided(score_exp, score_back, FDR = FDR,
                          nknockoff = n_permutation,
                          contrastScore_method = contrast_score,
                          FDR_control_method = procedure, ifpowerful = false, seed = seed)
      end
      FDR_nodisc = [length(re_i["discovery"]) == 0 for re_i in re["results"]]
      if any(FDR_nodisc .& (procedure != "aBH"))
          @warn string("At FDR = ", join(FDR[FDR_nodisc], ", "), ", no discovery has been found using current procedure. To make more discoveries, switch to aBH procedure or increase the FDR threshold.")
      end
  end
  contrast_score_value = re["contrastScore"]
  thre = [re_i["thre"] for re_i in re["results"]]
  discoveries = [re_i["discovery"] for re_i in re["results"]]
  q = re["results"][1]["q"]
  re = Dict("contrast_score" => contrast_score,
             "contrast_score_value" => contrast_score_value,
             "FDR" => FDR,
             "contrast_score_thre" => thre,
             "discoveries" => discoveries,
             "q" => q)
  return re
end
function clipper1sided(score_exp, score_back; FDR = 0.05, ifuseknockoff = nothing, nknockoff = nothing, contrastScore_method = nothing, importanceScore_method = "diff", FDR_control_method = nothing, ifpowerful = true, seed = 12345)
    if any(score_exp .< 0) || any(score_back .< 0)
        shift = minimum([minimum(score_exp[.!ismissing.(score_exp)]), minimum(score_exp[.!ismissing.(score_exp)])])
        score_exp .= score_exp .- shift
        score_back .= score_back .- shift
    end
    if ndims(score_exp) == 1
        score_exp = reshape(score_exp, (length(score_exp), 1))
    end
    if ndims(score_back) == 1
        score_back = reshape(score_back, (length(score_back), 1))
    end
    r1 = size(score_exp, 2)
    r2 = size(score_back, 2)
    if size(score_exp, 1) != size(score_back, 1)
        error("score_exp and score_back must have the same number of rows (features)")
    end
    if isnothing(ifuseknockoff)
        if r1 == r2
            ifuseknockoff = false
        else
            ifuseknockoff = true
        end
    end
    if !ifuseknockoff
        if !(r1 == r2)
            println("Caution: no knockoffs are constructed when the numbers of replicates are different; FDR control is not guaranteed")
        end
        if FDR_control_method == "GZ"
            println("FDR_control_method cannot be GZ when no knockoffs are constructed. Switching to BC.")
            FDR_control_method = "BC"
        end
        re = suppressWarnings(clipper_1sided_woknockoff(score_exp = score_exp,
                                    score_back = score_back,
                                    r1 = r1,
                                    r2 = r2,
                                    FDR = FDR,
                                    importanceScore_method = importanceScore_method,
                                    FDR_control_method = FDR_control_method))
        if ifpowerful && FDR_control_method != "BH"
            FDR_nodisc = [length(re_i.discovery) == 0 for re_i in re.results]
            if any(FDR_nodisc)
                println("At FDR = ", join(string.(FDR[FDR_nodisc]), ", "), ", no discovery has been found using FDR control method ", FDR_control_method, "; switching to BH...")
                re_clipperbh = clipper_BH(contrastScore = re.contrastScore, FDR = FDR[FDR_nodisc])
                re.results[FDR_nodisc] = re_clipperbh
            end
        end
    end
    if ifuseknockoff
        if r1 == 1 && r2 == 1
            error("Cannot generate knockoffs when both score_exp and score_back have one column. Please rerun clipper1sided by setting ifuseknockoff = F")
        end
        nknockoff_max = ifelse(r1 == r2, binomial(r1 + r2, r1) ÷ 2 - 1, binomial(r1 + r2, r1) - 1)
        if !isnothing(nknockoff)
            if nknockoff > nknockoff_max || !isa(nknockoff, Int) || nknockoff < 1
                println("nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the exhaustive knockoffs instead.")
            end
            nknockoff = min(nknockoff, nknockoff_max)
        else
            println("nknockoff is not supplied; generate the maximum number of knockoffs: ", nknockoff_max)
            if contrastScore_method == "diff"
                nknockoff = nknockoff_max
            else
                nknockoff = 1
            end
        end
        re = clipper_1sided_wknockoff(score_exp,score_back,r1,r2,FDR,importanceScore_method,contrastScore_method,nknockoff,nknockoff_max,seed)
        if ifpowerful && FDR_control_method != "BH"
            FDR_nodisc = [length(re_i.discovery) == 0 for re_i in re.results]
            if any(FDR_nodisc)
                println("At FDR = ", join(string.(FDR[FDR_nodisc]), ", "), ", no discovery has been found using FDR control method ", FDR_control_method, "; switching to BH...")
                re_clipperbh = clipper_BH(contrastScore = re.contrastScore, nknockoff = nknockoff, FDR = FDR[FDR_nodisc])
                re["results"][FDR_nodisc] = re_clipperbh
            end
        end
    end
    return re
end
function clipper2sided(score_exp, score_back, FDR = 0.05, nknockoff = nothing, contrastScore_method = "max", importanceScore_method = "diff", FDR_control_method = "GZ", ifpowerful = true, seed = 12345)
    if ndims(score_exp) == 1
        score_exp = reshape(score_exp, (length(score_exp), 1))
    end
    if ndims(score_back) == 1
        score_back = reshape(score_back, (length(score_back), 1))
    end
    r1 = size(score_exp, 2)
    r2 = size(score_back, 2)
    if r1 == 1 && r2 == 1
        error("clipper is not yet able to perform two sided identification when either condition has one replicate")
    end
    nknockoff_max = min(ifelse(r1 == r2, binomial(r1 + r2, r1) ÷ 2 - 1, binomial(r1 + r2, r1) - 1), 200)
    if !isnothing(nknockoff)
        if nknockoff > nknockoff_max || !isa(nknockoff, Int) || nknockoff < 1
            println("nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the maximal number of knockoffs instead.")
        end
        nknockoff = min(nknockoff, nknockoff_max)
    else
        if contrastScore_method == "diff"
            nknockoff = min(nknockoff_max, 50)
        else
            nknockoff = 1
        end
    end
    knockoffidx = generate_knockoffidx(r1 = r1, r2 = r2, nknockoff = nknockoff, nknockoff_max = nknockoff_max, seed = seed)
    kappatau_ls = compute_taukappa(score_exp = score_exp, score_back = score_back, r1, r2, true, knockoffidx, importanceScore_method, contrastScore_method)
    re = clipper_GZ(tau = kappatau_ls.tau, kappa = kappatau_ls.kappa, nknockoff = nknockoff, FDR = FDR)
    re = Dict("knockoffidx" => knockoffidx,
              "importanceScore_method" => importanceScore_method,
              "importanceScore" => kappatau_ls.importanceScore,
              "contrastScore_method" => contrastScore_method,
              "contrastScore" => (2 * kappatau_ls.kappa - 1) .* abs.(kappatau_ls.tau),
              "results" => re)
    if ifpowerful && FDR_control_method != "BH"
        FDR_nodisc = [length(re_i["discovery"]) == 0 for re_i in re["results"]]
        if any(FDR_nodisc .& contrastScore_method == "max")
            println("At FDR = ", join(string.(FDR[FDR_nodisc]), ", "), ", no discovery has been found using FDR control method ", FDR_control_method, " ; switching to BH...")
            re_clipperbh = clipper_BH(contrastScore = re["contrastScore"], nknockoff = nknockoff, FDR = FDR[FDR_nodisc])
            re["results"][FDR_nodisc] = re_clipperbh
        end
    end
    return re
end
function clipper_1sided_woknockoff(score_exp, score_back, r1, r2, FDR, aggregation_method, importanceScore_method, FDR_control_method)
    if r1 > 1
        score_exp = aggregate_clipper(score = score_exp, aggregation_method = aggregation_method)
    end
    if r2 > 1
        score_back = aggregate_clipper(score = score_back, aggregation_method = aggregation_method)
    end
    contrastscore = compute_importanceScore_wsinglerep(score_exp = score_exp, score_back = score_back, importanceScore_method = importanceScore_method)
    re = Dict()
    if FDR_control_method == "BC"
        re = suppressWarnings(clipper_BC(contrastScore = contrastscore, FDR = FDR))
    end
    if FDR_control_method == "BH"
        re = clipper_BH(contrastScore = contrastscore, FDR = FDR)
    end
    re = Dict("importanceScore" => contrastscore,
              "importanceScore_method" => importanceScore_method,
              "contrastScore" => contrastscore,
              "contrastScore_method" => importanceScore_method,
              "results" => re)
    return re
end
function clipper_1sided_wknockoff(score_exp, score_back, r1, r2, FDR, importanceScore_method, contrastScore_method, nknockoff, nknockoff_max, seed)
    knockoffidx = generate_knockoffidx(r1, r2, nknockoff, nknockoff_max, seed)
    kappatau_ls = compute_taukappa(score_exp, score_back, r1, r2, false, knockoffidx, importanceScore_method, contrastScore_method)
    re = clipper_GZ(kappatau_ls["tau"], kappatau_ls["kappa"], nknockoff, FDR)
    re = Dict("knockoffidx" => knockoffidx,
              "importanceScore_method" => importanceScore_method,
              "importanceScore" => kappatau_ls["importanceScore"],
              "contrastScore_method" => contrastScore_method,
              "contrastScore" => (2 .* kappatau_ls["kappa"] .- 1) .* abs.(kappatau_ls["tau"]),
              "results" => re)
    return re
end
function aggregate_clipper(score, aggregation_method)
    if aggregation_method == "mean"
        score_single = [Statistics.mean(x[.!ismissing.(x)]) for x in eachrow(score)]
    elseif aggregation_method == "median"
        score_single = [median(x[.!ismissing.(x)]) for x in eachrow(score)]
    end
    return score_single
end
function compute_importanceScore_wsinglerep(score_exp, score_back, importanceScore_method)
    if importanceScore_method == "diff"
        contrastScore = score_exp - score_back
    elseif importanceScore_method == "max"
        contrastScore = max.(score_exp, score_back) .* sign.(score_exp - score_back)
    end
    return vec(contrastScore)
end
function clipper_BC(contrastScore, FDR)
    contrastScore[isna.(contrastScore)] .= 0 
    c_abs = abs.(contrastScore[contrastScore .!= 0])
    c_abs = sort(unique(c_abs))
    i = 1
    emp_fdp = fill!(similar(c_abs), NaN)
    emp_fdp[1] = 1
    while i <= length(c_abs)
        t = c_abs[i]
        emp_fdp[i] = min((1 + sum(contrastScore .<= -t)) / sum(contrastScore .>= t), 1)
        if i >= 2
            emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])
        end
        i += 1
    end
    c_abs = c_abs[.!isna.(emp_fdp)]
    emp_fdp = emp_fdp[.!isna.(emp_fdp)]
    q = emp_fdp[match(contrastScore, c_abs)]
    q[isna.(q)] .= 1
    re = []
    for FDR_i in FDR
        thre = c_abs[minimum(findall(x -> x <= FDR_i, emp_fdp))]
        re_i = Dict("FDR" => FDR_i,
                    "FDR_control" => "BC",
                    "thre" => thre,
                    "q" => q,
                    "discovery" => findall(x -> x >= thre, contrastScore))
        push!(re, re_i)
    end
    return re
end
function generate_knockoffidx(r1, r2, nknockoff, nknockoff_max, seed)
    Random.seed!(seed)
    if nknockoff_max == 200
        knockoffidx = Vector{Vector{Int}}(undef, nknockoff)
        i_knockoff = 1
        while i_knockoff <= nknockoff
            temp = sample(r1 + r2, r1, replace = false)
            if any(!in(1:r1, temp)) && any(in(1:r1, temp))
                knockoffidx[i_knockoff] = temp
                i_knockoff += 1
            else
                continue
            end
        end
    else
        combination_all = collect(combinations(1:r1+r2, r1))
        if r1 == r2
            combination_all = combination_all[1:Int(length(combination_all) * 0.5)][2:end]
        else
            combination_all = combination_all[2:end]
        end
        knockoffidx = sample(combination_all, nknockoff, replace = false)
    end
    return knockoffidx
end
function compute_taukappa(score_exp, score_back, r1, r2, if2sided, knockoffidx, importanceScore_method, contrastScore_method)
    perm_idx = vcat(collect(1:r1), knockoffidx...)
    score_tot = hcat(score_exp, score_back)
    imp_ls = zeros(size(score_tot,1),length(perm_idx))
    i = 0
    for x in perm_idx
        i += 1
        se = score_tot[:, x]
        sb = score_tot[:, setdiff(1:(r1+r2), x)]
        se = aggregate_clipper(se, "mean")
        sb = aggregate_clipper(sb, "mean")
        imp = compute_importanceScore_wsinglerep(se, sb, importanceScore_method)
        if if2sided
            imp = abs.(imp)
        end
        imp_ls[:,i] = imp
    end
    kappatau_ls = []
    for i in range(1,size(imp_ls,1))
        x = imp_ls[i,:]
        kappa = !any(x[2:end] .== maximum(x))
        if isempty(kappa)
            kappa = NaN
        end
        x_sorted = sort(x, rev = true)
        if contrastScore_method == "diff"
            tau = x_sorted[1] - x_sorted[2]
        elseif contrastScore_method == "max"
            tau = x_sorted[1]
        end
        push!(kappatau_ls, Dict("kappa" => kappa, "tau" => tau))
    end
    re = Dict("importanceScore" => imp_ls,
              "kappa" => [x["kappa"] for x in kappatau_ls],
              "tau" => [x["tau"] for x in kappatau_ls])
    return re
end
function clipper_GZ(tau, kappa, nknockoff, FDR)
    n = length(tau)
    contrastScore = (2 .* kappa .- 1) .* abs.(tau)
    contrastScore[ismissing.(contrastScore)] .= 0
    c_abs = sort(unique(abs.(contrastScore[contrastScore .!= 0])))
    emp_fdp = fill(NaN,length(c_abs))
    emp_fdp[1] = 1
    i = 1
    while i <= length(c_abs)
        t = c_abs[i]
        emp_fdp[i] = minimum([(1/nknockoff + 1/nknockoff * sum(contrastScore .<= -t)) / sum(contrastScore .>= t), 1])
        if i >= 2
            emp_fdp[i] = minimum([emp_fdp[i], emp_fdp[i-1]])
        end
        i += 1
    end
    c_abs = c_abs[.!ismissing.(emp_fdp)]
    emp_fdp = emp_fdp[.!ismissing.(emp_fdp)]
    q = fill(1.0,length(contrastScore))
    matchidxs = [findfirst(c_abs .== value) for value in contrastScore] 
    q[.!isnothing.(matchidxs)] .= emp_fdp[matchidxs[.!isnothing.(matchidxs)]]
    re = []
    for FDR_i in FDR
        if !isnothing(findfirst(emp_fdp .<= FDR_i))
            thre = c_abs[findfirst(emp_fdp .<= FDR_i)]
            re_i = Dict("FDR" => FDR_i,
                        "FDR_control" => "BC",
                        "thre" => thre,
                        "discovery" => findall(contrastScore .>= thre),
                        "q" => q)
        else
            re_i = Dict("FDR" => FDR_i,
                        "FDR_control" => "BC",
                        "thre" => NaN,
                        "discovery" => Int[],
                        "q" => q)
        end
        push!(re, re_i)
    end
    return re
end
