struct PartialV{T}  
    vec_ratio::Vector{T}
    eigvals::Vector{T}
    eigvecs::Matrix{T}
end
struct LowRankApproxUDV{T}  
    U::Matrix{T}
    D::Vector{T}
    V::Matrix{T}
    M::Matrix{T}
    method::String
    rank::Int
    _2_norm::Float64
    corr::Float64
end
function LowRankApproxUDV(K2::Matrix{T},approx_method::String,approx_rank::Int,avg_2norm_threshold::AbstractFloat,max_rank::Number) where T
    _2_norm = 1e6
    test_rank = copy(approx_rank)
    iter = 1
    K2U = nothing
    K2D = nothing
    K2V = nothing
    M = nothing
    corr = NAN
    if max_rank isa AbstractFloat
        max_rank_prop = max_rank
        max_rank = Int(round(size(K2,1) * max_rank_prop))
    end
    while true
        if approx_method == "aca"
            K2U, K2V = aca(K2,test_rank)
            K2D = ones(T, test_rank)
            K2V = K2V' |> Matrix
        elseif approx_method == "psvd"
            opts = LRAOptions(rank=test_rank, rtol=1e-20)
            F = psvdfact(K2,opts) 
            K2U = F.U
            K2D = F.S
            K2V = F.Vt
        elseif approx_method == "cur"
            opts = LRAOptions(rank=test_rank, rtol=1e-20)
            U = curfact(K2,opts) 
            F = CUR(K2, U)
            K2U = F.C
            K2D = F.U
            K2V = F.C'
        elseif approx_method == "pheig"
            opts = LRAOptions(rank=test_rank, rtol=1e-20)
            F = pheigfact(K2,opts) 
            K2U = F.vectors
            K2D = F.values
            K2V = F.vectors' |> Matrix
        end
        M = K2U * diagm(K2D) * K2V
        corr = cor(M[:],K2[:])
        _2_norm = norm(M - K2)
        avg_2_norm = _2_norm^2 / length(M)
        println("LRA iter: ",iter, ", approx_rank: ",test_rank,", _2_norm: ", _2_norm, ", avg_2_norm: ", avg_2_norm)
        if (avg_2_norm < avg_2norm_threshold) | (test_rank >= max_rank)
            break
        else
            alpha = _2_norm / sqrt(avg_2norm_threshold * length(M))
            alpha *= 1.5
            test_rank = Int(round(alpha * test_rank))
            if test_rank > max_rank
                test_rank = max_rank
            end
        end
        iter += 1
    end
    test_rank = length(K2D)
    println("Final LRA rank: ", test_rank)
    LowRankApproxUDV{eltype(K2U)}(K2U, K2D, K2V, M, approx_method, test_rank, _2_norm, corr)
end
function reml(VCM_Arr::Vector{Matrix{T}}, y::Vector{T}, X::Matrix{T}, Vi_X::Matrix{T}, Xt_Vi_X_i::Matrix{T}, Hi::Matrix{T}, Py::Vector{T}, P::Matrix{T}, PA::Matrix{T}, Vi::Matrix{T}, cVi::Matrix{T}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvec::Matrix{T}; pred_rand_eff::Bool=false, reml_max_iter::Int=100, EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing, varcmp::Union{Nothing,Vector{T}}=nothing, partialV::Union{Nothing,PartialV{Float32}}=nothing, LRA_factors::Union{Nothing,LowRankApproxUDV{T}}=nothing, em_update::Bool=false, convergence_parameters::T=1e-8, convergence_dlogL::T=NAN) where T
    @timeit to "Prepare for reml_iteration" begin
        _n = length(y) 
        Var_y = var(y)
        _n_cmp = length(VCM_Arr)
        _X_c = size(X, 2) 
        fill!(Vi_X,0)
        fill!(Xt_Vi_X_i,0)
        fill!(Hi,0)
        fill!(Py,0)
        fill!(P,0)
        fill!(PA,0)
        fill!(Vi,0)
        fill!(cVi,0)
        fill!(cVi,0)
        fill!(Di,0)
        fill!(EAvec_Di,0)
        no_constrain = false
        prior_var_flag = !isnothing(varcmp)
        if !prior_var_flag
            varcmp = zeros(T, _n_cmp)
            varcmp .= Var_y / _n_cmp
        else
            if em_update
                prior_var_flag = false
            else
                prior_var_flag = true
            end
        end
    end
    lgL, converged_flag = reml_iteration(VCM_Arr, y, X, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, P, PA, Vi, cVi, Di, EAvec_Di, EAvec; prior_var_flag=prior_var_flag, reml_max_iter=reml_max_iter, EA=EA, partialV=partialV, LRA_factors=LRA_factors, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
    @timeit to "Output b_bat and Hsq" begin
        if pred_rand_eff 
            u = zeros(T, _n, _n_cmp)
            for i in 1:_n_cmp
                u[:, i] .= mat_mul(VCM_Arr[i], Py) * varcmp[i]
            end
        else
            u = nothing
        end
        b_bat = zeros(T, _X_c)
        mat_mul!(b_bat, mat_mul(Xt_Vi_X_i, Vi_X'), y)
        Hsq = zeros(T, _n_cmp - 1)
        VarHsq = zeros(T, _n_cmp - 1)
        Vp, VarVp = calcu_Vp(varcmp, Hi)
        for i in eachindex(Hsq)
            Hsq[i], VarHsq[i] = calcu_hsq(i, Vp, VarVp, varcmp, Hi)
        end
    end
    abc_uni = Dict(
            :ΣG => varcmp[1:(_n_cmp-1)],
            :ΣG_se => nothing,
            :Σe => varcmp[_n_cmp],
            :Σe_se => nothing,
            :Σp => Vp,
            :Σp_se => VarVp,
            :h2 => Hsq,
            :h2_se => VarHsq,
            :logL => lgL,
            :logL0 => nothing,
            :n => _n,
            :Fix_eff => b_bat,
            :Fix_eff_se => nothing,
            :u => u,
            :converged => converged_flag
            )
    return abc_uni
end
function reml(VCM_Arr::Vector{Matrix{T}}, y::Vector{T}, X::Matrix{T}; pred_rand_eff::Bool=false, reml_max_iter::Int=100, EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing, varcmp::Union{Nothing,Vector{T}}=nothing, partialV::Union{Nothing,PartialV{Float32}}=nothing, LRA_factors::Union{Nothing,LowRankApproxUDV{T}}=nothing, em_update::Bool=false, convergence_parameters::T=1e-6, convergence_dlogL::T=NAN) where T
    @timeit to "Prepare for reml_iteration" begin
        _n = length(y) 
        Var_y = var(y)
        _n_cmp = length(VCM_Arr)
        _X_c = size(X, 2) 
        Vi_X = zeros(T, _n, _X_c)
        Xt_Vi_X_i = zeros(T, _X_c, _X_c)
        Hi = zeros(T, _n_cmp, _n_cmp)
        Py = zeros(T, _n)
        no_constrain = false
        prior_var_flag = !isnothing(varcmp)
        if !prior_var_flag
            varcmp = zeros(T, _n_cmp)
            varcmp .= Var_y / _n_cmp
        else
            if em_update
                prior_var_flag = false
            else
                prior_var_flag = true
            end
        end
    end
    lgL, converged_flag = reml_iteration(VCM_Arr, y, X, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, prior_var_flag=prior_var_flag, reml_max_iter=reml_max_iter, EA=EA, partialV=partialV, LRA_factors=LRA_factors, convergence_parameters=convergence_parameters, convergence_dlogL=convergence_dlogL)
    @timeit to "Output b_bat and Hsq" begin
        if pred_rand_eff 
            u = zeros(T, _n, _n_cmp)
            for i in 1:_n_cmp
                u[:, i] .= mat_mul(VCM_Arr[i], Py) * varcmp[i]
            end
        else
            u = nothing
        end
        b_bat = zeros(T, _X_c)
        mat_mul!(b_bat, mat_mul(Xt_Vi_X_i, Vi_X'), y)
        Hsq = zeros(T, _n_cmp - 1)
        VarHsq = zeros(T, _n_cmp - 1)
        Vp, VarVp = calcu_Vp(varcmp, Hi)
        for i in eachindex(Hsq)
            Hsq[i], VarHsq[i] = calcu_hsq(i, Vp, VarVp, varcmp, Hi)
        end
    end
    abc_uni = Dict(
            :ΣG => varcmp[1:(_n_cmp-1)],
            :ΣG_se => nothing,
            :Σe => varcmp[_n_cmp],
            :Σe_se => nothing,
            :Σp => Vp,
            :Σp_se => VarVp,
            :h2 => Hsq,
            :h2_se => VarHsq,
            :logL => lgL,
            :logL0 => nothing,
            :n => _n,
            :Fix_eff => b_bat,
            :Fix_eff_se => nothing,
            :u => u,
            :converged => converged_flag
            )
    return abc_uni
end
function reml_iteration(VCM_Arr::Vector{Matrix{T}}, y::Vector{T}, X::Matrix{T}, Vi_X::Matrix{T}, Xt_Vi_X_i::Matrix{T}, Hi::Matrix{T}, Py::Vector{T}, varcmp::Vector{T}, P::Matrix{T}, PA::Matrix{T}, Vi::Matrix{T}, cVi::Matrix{T}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvec::Matrix{T}; prior_var_flag::Bool=false, reml_max_iter::Int=100, no_constrain::Bool=false, EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing, partialV::Union{Nothing,PartialV{Float32}}=nothing, LRA_factors::Union{Nothing,LowRankApproxUDV{T}}=nothing, convergence_parameters::T=1e-6, convergence_dlogL::T=NAN) where T
    if isnan(convergence_dlogL)
        is_calcu_logdet = false
    else
        is_calcu_logdet = true
    end
    @timeit to "Init reml_iteration" begin
        mtd_str = ["AI-REML", "EM-REML"]
        _n_cmp = length(VCM_Arr)
        _var_name = ["Ve"]
        for i in (_n_cmp-1):-1:1
            pushfirst!(_var_name,string("Vg",i))
        end
        _n = length(y)
        Var_y = var(y)
        if !isnothing(LRA_factors) 
            approx_rank = LRA_factors.rank
            Ai_U = zeros(T, _n, approx_rank)
            V_Ai = zeros(T, approx_rank, _n)
            V_Ai_U = zeros(T, approx_rank, approx_rank)
            Ai_U_Di_V_Ai_U_i = zeros(T, _n, approx_rank)
            Ai_U_Di_V_Ai_U_i_V_Ai = zeros(T, _n, _n)
        end
        if !isnothing(EA)
            EAvect = EA.vectors' |> Matrix
        end
        _reml_mtd = 1
        constrain_num = 0; iter = 0; reml_mtd_tmp = copy(_reml_mtd)
        logdet_V = 0.0; logdet_Xt_Vi_X = 0.0; prev_lgL = -1e20; lgL = -1e20; dlogL = 1000.0; dpars = 100.0; dgra = 100.0
        prev_dlogL = 1000.0; unstable_dlogL = false; extreme_dlogL = false; n_unreliable = 0
        prev_prev_varcmp = copy(varcmp)
        prev_varcmp = copy(varcmp)
        varcomp_init = copy(varcmp)
        converged_flag = false
    end
    while iter <= reml_max_iter
        if iter == 0
            prev_varcmp .= varcomp_init
            if prior_var_flag
                print("Prior values of variance components: ")
                for i in 1:_n_cmp
                    print(varcmp[i], "\t")
                end
                println()
            else
                _reml_mtd = 2
                println("Calculating prior values of variance components by EM-REML ...")
            end
        end
        if iter == 1
            _reml_mtd = copy(reml_mtd_tmp)
            print("Prior values of variance components: ")
            for i in 1:_n_cmp
                print(varcmp[i], "\t")
            end
            println()
            if is_calcu_logdet
                print(string("Running ", mtd_str[_reml_mtd], " algorithm ...", "\nIter.\tlogL\tdpars\t"))
            else
                print(string("Running ", mtd_str[_reml_mtd], " algorithm ...", "\nIter.\tdpars\t"))
            end
            for i in 1:_n_cmp
                print(_var_name[i], "\t")
            end
            println()
        end
        if (_n_cmp == 2) && !isnothing(EA)
            @timeit to "calcu_Vi! EAvec_Di" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, EA, Di, EAvec_Di, EAvect, cVi, is_calcu_logdet=is_calcu_logdet)
        elseif (_n_cmp == 3) && !isnothing(LRA_factors) 
            @timeit to "calcu_Vi! Woodbury" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, EA, Di, EAvec_Di, EAvect, LRA_factors, Ai_U, V_Ai, V_Ai_U, Ai_U_Di_V_Ai_U_i, Ai_U_Di_V_Ai_U_i_V_Ai, is_calcu_logdet=is_calcu_logdet)
        elseif (_n_cmp == 3) && !isnothing(partialV)
            @timeit to "calcu_Vi! Pre-calculate" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, partialV, Di, EAvec_Di, EAvec, is_calcu_logdet=is_calcu_logdet)
        else
            @timeit to "calcu_Vi!" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, is_calcu_logdet=is_calcu_logdet)
        end
        if isinf(logdet_V)
            break
        end
        @timeit to "calcu_P!" logdet_Xt_Vi_X = calcu_P!(Vi, Vi_X, Xt_Vi_X_i, P, X, is_calcu_logdet=is_calcu_logdet)
        if isinf(logdet_Xt_Vi_X)
            break
        end
        if _reml_mtd == 1
            @timeit to "ai_reml!" delta, R = ai_reml!(PA, P, Hi, Py, prev_varcmp, varcmp, dlogL, VCM_Arr, y, is_calcu_logdet=is_calcu_logdet, prev_prev_varcmp=prev_prev_varcmp)
            if isnan(delta[1])
                break
            end
        elseif _reml_mtd == 2
            @timeit to "em_reml!" R = em_reml!(PA, P, Py, prev_varcmp, varcmp, VCM_Arr, y)  
        end
        @timeit to "lgL" lgL = -0.5 * (logdet_Xt_Vi_X + logdet_V + dot(y, Py))
        if !no_constrain 
            @timeit to "constrain_varcmp!" constrain_num = constrain_varcmp!(varcmp, y)
        end
        if iter >= 1
            print(iter,"\t")
            @runif is_calcu_logdet print(lgL,"\t")
            print(dpars,"\t")
            for i in 1:_n_cmp
                print(varcmp[i], "\t")
            end
            println()
        end
        dlogL = lgL - prev_lgL
        dpars = norm(varcmp - prev_varcmp)^2 / norm(varcmp)^2
        dgra = norm(R)
        con_1 = dpars < convergence_parameters
        con_2 = abs(dlogL) < convergence_dlogL || (abs(dlogL) < 1e-2 && dlogL < 0)
        con_3 = dgra < 1e-6
        if con_1 || con_2 || con_3
            converged_flag = true
            if con_1
                println("* Parameters change converged.")
            end
            if con_2
                println("* Log-likelihood ratio converged.")
            end
            if con_3
                println("* Partial derivative converged.")
            end
            break
        end
        if iter > ifelse(is_calcu_logdet, 2, 0)
            if sum(varcmp) > 2 * Var_y
                varcmp .= prev_varcmp
                @runif _args_debug println("* sum(varcmp) > 2 * Var_y")
                break
            end
            if abs(lgL) > 1e4
                varcmp .= prev_varcmp
                @runif _args_debug println("* abs(lgL) > 1e4")
                break
            end
            prev_r_t = (prev_varcmp ./ sum(prev_varcmp)) .< 0.0001
            r_t = (varcmp ./ sum(varcmp)) .< 0.0001 
            if sum(r_t .& prev_r_t) > 0
                varcmp[r_t .& prev_r_t] .= 0
            end
                extreme_dlogL = abs(dlogL) > abs(prev_dlogL)
                unstable_dlogL = (dlogL * prev_dlogL) < 0
            if extreme_dlogL | unstable_dlogL
                n_unreliable += 1
                println("Warning: Unstable iteration (n = ",n_unreliable,").")
            end
            unstable_dlogL = false
            extreme_dlogL = false
            if n_unreliable == 3
                @runif _args_debug println("* n_unreliable == 3")
                break
            end
        end
        prev_prev_varcmp .= prev_varcmp
        prev_varcmp .= varcmp
        prev_lgL = copy(lgL)
        prev_dlogL = copy(dlogL)
        iter += 1
    end
    if converged_flag
        println("* AIREML converged.")
    else
        println("* AIREML not converged.")
    end
    return lgL, converged_flag
end
function reml_iteration(VCM_Arr::Vector{Matrix{T}}, y::Vector{T}, X::Matrix{T}, Vi_X::Matrix{T}, Xt_Vi_X_i::Matrix{T}, Hi::Matrix{T}, Py::Vector{T}, varcmp::Vector{T}; prior_var_flag::Bool=false, reml_max_iter::Int=100, no_constrain::Bool=false, EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing, partialV::Union{Nothing,PartialV{Float32}}=nothing, LRA_factors::Union{Nothing,LowRankApproxUDV{T}}=nothing, convergence_parameters::T=1e-8, convergence_dlogL::T=NAN) where T
    if isnan(convergence_dlogL)
        is_calcu_logdet = false
    else
        is_calcu_logdet = true
    end
    @timeit to "Init reml_iteration" begin
        mtd_str = ["AI-REML", "EM-REML"]
        _n_cmp = length(VCM_Arr)
        _var_name = ["Ve"]
        for i in (_n_cmp-1):-1:1
            pushfirst!(_var_name,string("Vg",i))
        end
        _n = length(y)
        Var_y = var(y)
        _X_c = size(X, 2)
        P = zeros(T, _n, _n)
        PA = zeros(T, _n, _n)
        Vi = zeros(T, _n, _n)
        cVi = zeros(T, _n, _n)
        if !isnothing(LRA_factors) 
            approx_rank = LRA_factors.rank
            Ai_U = zeros(T, _n, approx_rank)
            V_Ai = zeros(T, approx_rank, _n)
            V_Ai_U = zeros(T, approx_rank, approx_rank)
            Ai_U_Di_V_Ai_U_i = zeros(T, _n, approx_rank)
            Ai_U_Di_V_Ai_U_i_V_Ai = zeros(T, _n, _n)
        end
        if !isnothing(EA)
            EAvect = EA.vectors' |> Matrix
        end
            Di = Diagonal(zeros(T, _n, _n))
            EAvec_Di = zeros(T, _n, _n)
            EAvec = zeros(T, _n, _n)
        _reml_mtd = 1
        constrain_num = 0; iter = 0; reml_mtd_tmp = copy(_reml_mtd)
        logdet_V = 0.0; logdet_Xt_Vi_X = 0.0; prev_lgL = -1e20; lgL = -1e20; dlogL = 1000.0; dpars = 100.0; dgra = 100.0
        prev_dlogL = 1000.0; unstable_dlogL = false; extreme_dlogL = false; n_unreliable = 0
        prev_prev_varcmp = copy(varcmp)
        prev_varcmp = copy(varcmp)
        varcomp_init = copy(varcmp)
        converged_flag = false
    end
    while iter <= reml_max_iter
        if iter == 0
            prev_varcmp .= varcomp_init
            if prior_var_flag
                print("Prior values of variance components: ")
                for i in 1:_n_cmp
                    print(varcmp[i], "\t")
                end
                println()
            else
                _reml_mtd = 2
                println("Calculating prior values of variance components by EM-REML ...")
            end
        end
        if iter == 1
            _reml_mtd = copy(reml_mtd_tmp)
            print("Prior values of variance components: ")
            for i in 1:_n_cmp
                print(varcmp[i], "\t")
            end
            println()
            if is_calcu_logdet
                print(string("Running ", mtd_str[_reml_mtd], " algorithm ...", "\nIter.\tlogL\tdpars\t"))
            else
                print(string("Running ", mtd_str[_reml_mtd], " algorithm ...", "\nIter.\tdpars\t"))
            end
            for i in 1:_n_cmp
                print(_var_name[i], "\t")
            end
            println()
        end
        if (_n_cmp == 2) && !isnothing(EA)
            @timeit to "calcu_Vi! EAvec_Di" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, EA, Di, EAvec_Di, EAvect, cVi, is_calcu_logdet=is_calcu_logdet)
        elseif (_n_cmp == 3) && !isnothing(LRA_factors) 
            @timeit to "calcu_Vi! Woodbury" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, EA, Di, EAvec_Di, EAvect, LRA_factors, Ai_U, V_Ai, V_Ai_U, Ai_U_Di_V_Ai_U_i, Ai_U_Di_V_Ai_U_i_V_Ai, is_calcu_logdet=is_calcu_logdet)
        elseif (_n_cmp == 3) && !isnothing(partialV)
            @timeit to "calcu_Vi! Pre-calculate" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, partialV, Di, EAvec_Di, EAvec, is_calcu_logdet=is_calcu_logdet)
        else
            @timeit to "calcu_Vi!" logdet_V = calcu_Vi!(Vi, VCM_Arr, prev_varcmp, is_calcu_logdet=is_calcu_logdet)
        end
        if isinf(logdet_V)
            break
        end
        @timeit to "calcu_P!" logdet_Xt_Vi_X = calcu_P!(Vi, Vi_X, Xt_Vi_X_i, P, X, is_calcu_logdet=is_calcu_logdet)
        if isinf(logdet_Xt_Vi_X)
            break
        end
        if _reml_mtd == 1
            @timeit to "ai_reml!" delta, R = ai_reml!(PA, P, Hi, Py, prev_varcmp, varcmp, dlogL, VCM_Arr, y, is_calcu_logdet=is_calcu_logdet, prev_prev_varcmp=prev_prev_varcmp)
            if isnan(delta[1])
                break
            end
        elseif _reml_mtd == 2
            @timeit to "em_reml!" R = em_reml!(PA, P, Py, prev_varcmp, varcmp, VCM_Arr, y)  
        end
        @timeit to "lgL" lgL = -0.5 * (logdet_Xt_Vi_X + logdet_V + dot(y, Py))
        if !no_constrain 
            @timeit to "constrain_varcmp!" constrain_num = constrain_varcmp!(varcmp, y)
        end
        if iter >= 1
            print(iter,"\t")
            @runif is_calcu_logdet print(lgL,"\t")
            print(dpars,"\t")
            for i in 1:_n_cmp
                print(varcmp[i], "\t")
            end
            println()
        end
        dlogL = lgL - prev_lgL
        dpars = norm(varcmp - prev_varcmp)^2 / norm(varcmp)^2
        dgra = norm(R)
        con_1 = dpars < convergence_parameters
        con_2 = abs(dlogL) < convergence_dlogL || (abs(dlogL) < 1e-2 && dlogL < 0)
        con_3 = dgra < 1e-6
        if con_1 || con_2 || con_3
            converged_flag = true
            if con_1
                println("* Parameters change converged.")
            end
            if con_2
                println("* Log-likelihood ratio converged.")
            end
            if con_3
                println("* Partial derivative converged.")
            end
            break
        end
        if iter > ifelse(is_calcu_logdet, 2, 0)
            if sum(varcmp) > 2 * Var_y
                varcmp .= prev_varcmp
                @runif _args_debug println("* sum(varcmp) > 2 * Var_y")
                break
            end
            if abs(lgL) > 1e4
                varcmp .= prev_varcmp
                @runif _args_debug println("* abs(lgL) > 1e4")
                break
            end
            prev_r_t = (prev_varcmp ./ sum(prev_varcmp)) .< 0.0001
            r_t = (varcmp ./ sum(varcmp)) .< 0.0001 
            if sum(r_t .& prev_r_t) > 0
                varcmp[r_t .& prev_r_t] .= 0
            end
                extreme_dlogL = abs(dlogL) > abs(prev_dlogL)
                unstable_dlogL = (dlogL * prev_dlogL) < 0
            if extreme_dlogL | unstable_dlogL
                n_unreliable += 1
                println("Warning: Unstable iteration (n = ",n_unreliable,").")
            end
            unstable_dlogL = false
            extreme_dlogL = false
            if n_unreliable == 3
                @runif _args_debug println("* n_unreliable == 3")
                break
            end
        end
        prev_prev_varcmp .= prev_varcmp
        prev_varcmp .= varcmp
        prev_lgL = copy(lgL)
        prev_dlogL = copy(dlogL)
        iter += 1
    end
    if converged_flag
        println("* AIREML converged.")
    else
        println("* AIREML not converged.")
    end
    return lgL, converged_flag
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}, partialV::Union{Nothing,PartialV{Float32}}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvec::Matrix{T}; is_calcu_logdet::Bool=true, verbose::Bool=false) where T
    ratio = maximum([prev_varcmp[1] / prev_varcmp[2], prev_varcmp[2] / prev_varcmp[1]])
    @runif verbose println("Ratio between vg1 and vg2: ", ratio)
    if ratio > maximum(partialV.vec_ratio)
        @timeit to "Vi-Excessive ratio" logdet_V = calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T})
    else 
        if is_calcu_logdet
            _n_cmp = length(VCM_Arr)
            @timeit to "Vi3 1.0" fill!(Vi, 0)
            Vi[diagind(Vi)] .= prev_varcmp[_n_cmp]
            @timeit to "Vi3 1.1" for i in 1:(_n_cmp-1)
                Vi .+= VCM_Arr[i] .* prev_varcmp[i]
            end
            @timeit to "Vi3 1.2" begin
                cholesky!(Vi)
                logdet_V = logdet(UpperTriangular(Vi)) * 2
            end
        else
            logdet_V = NAN
        end
        ratio_index = argmin(abs.(ratio .- partialV.vec_ratio)) * 2
        if prev_varcmp[1] >= prev_varcmp[2]
            ratio_index -= 1
        end
        j = ratio_index
        n_samples = size(Vi,1)
        range_inds = n_samples*(j-1)+1:n_samples*j
        vals_inds = CartesianIndices((range_inds,))
        vecs_inds = CartesianIndices((range_inds,1:n_samples))
        @timeit to "Vi3 EAvec" copyto!(EAvec, CartesianIndices((1:n_samples,1:n_samples)), partialV.eigvecs, vecs_inds)
        @timeit to "Vi3 1.3" begin
            Di[diagind(Di)] .= partialV.eigvals[vals_inds]
            Di[diagind(Di)] .*= sum(prev_varcmp[1:2])
            Di[diagind(Di)] .+= prev_varcmp[3]
            Di[diagind(Di)] .\= 1.0
        end
        @timeit to "Vi3 1.4" mul!(EAvec_Di, Di, EAvec')
        @timeit to "Vi3 1.5" mat_mul!(Vi, EAvec, EAvec_Di)
    end
    return logdet_V
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}; is_calcu_logdet::Bool=true) where T
    _n_cmp = length(VCM_Arr)
    @timeit to "Vi 1.0" fill!(Vi, 0)
    Vi[diagind(Vi)] .= prev_varcmp[_n_cmp]
    @timeit to "Vi 1.1" for i in 1:(_n_cmp-1)
        Vi .+= VCM_Arr[i] .* prev_varcmp[i]
    end
    logdet_V = NAN
    if isposdef(Vi)
        @timeit to "Vi 1.2" begin
            cVi = cholesky!(Vi)
            @runif is_calcu_logdet logdet_V, _ = logabsdet(cVi)
        end
        @timeit to "Vi 1.3" Vi .= inv(cVi)
    else
        @timeit to "Vi 1.2" begin
            @runif is_calcu_logdet logdet_V, _ = logabsdet(Vi)
        end
        try
            @timeit to "Vi 1.3" Vi .= inv(Vi)
        catch e 
            return Inf
        end
    end
    return logdet_V
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}, EA::Eigen{T, T, Matrix{T}, Vector{T}}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvect::Matrix{T}, LRA_factors::LowRankApproxUDV{T}; is_calcu_logdet::Bool=true) where T
    if is_calcu_logdet
        _n_cmp = length(VCM_Arr)
        @timeit to "Vi 1.0" begin
            fill!(Vi, 0)
            Vi[diagind(Vi)] .= prev_varcmp[_n_cmp]
        end
        @timeit to "Vi 1.1" for i in 1:(_n_cmp-1)
            Vi .+= VCM_Arr[i] .* prev_varcmp[i]
        end
        @timeit to "Vi 1.2" logdet_V, _ = logabsdet(Vi)
    else
        logdet_V = NAN
    end
    @timeit to "Vi 1.3" Di[diagind(Di)] .= 1.0 ./ (prev_varcmp[1] .* EA.values .+ prev_varcmp[3])
    @timeit to "Vi 1.4" mul!(EAvec_Di, Di, EAvect)
    @timeit to "Vi 1.5" mat_mul!(Vi, EA.vectors, EAvec_Di) 
    if prev_varcmp[2] > 0
        @timeit to "Vi 1.6" lDi = Diagonal(1.0 ./ (prev_varcmp[2] .* LRA_factors.D))
        @timeit to "Vi 1.7" Ai_U = mat_mul(Vi, LRA_factors.U) 
        @timeit to "Vi 1.8" V_Ai = mat_mul(LRA_factors.V, Vi) 
        @timeit to "Vi 1.9" Di_V_Ai_U_i = inv(lDi + mat_mul(LRA_factors.V, Ai_U))
        @timeit to "Vi 1.10" Ai_U_DDinv = mat_mul(Ai_U, Di_V_Ai_U_i)
        @timeit to "Vi 1.11" Vi .-= mat_mul(Ai_U_DDinv, V_Ai) 
    end
    return logdet_V
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}, EA::Eigen{T, T, Matrix{T}, Vector{T}}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvect::Matrix{T}, LRA_factors::LowRankApproxUDV{T}, Ai_U::Matrix{T}, V_Ai::Matrix{T}, V_Ai_U::Matrix{T}, Ai_U_Di_V_Ai_U_i::Matrix{T}, Ai_U_Di_V_Ai_U_i_V_Ai::Matrix{T}; is_calcu_logdet::Bool=true) where T
    if is_calcu_logdet
        _n_cmp = length(VCM_Arr)
        @timeit to "Vi 1.0" begin
            fill!(Vi, 0)
            Vi[diagind(Vi)] .= prev_varcmp[_n_cmp]
        end
        @timeit to "Vi 1.1" for i in 1:(_n_cmp-1)
            Vi .+= VCM_Arr[i] .* prev_varcmp[i]
        end
        @timeit to "Vi 1.2" logdet_V, _ = logabsdet(Vi)
    else
        logdet_V = NAN
    end
    @timeit to "Vi 1.3" Di[diagind(Di)] .= 1.0 ./ (prev_varcmp[1] .* EA.values .+ prev_varcmp[3])
    @timeit to "Vi 1.4" mul!(EAvec_Di, Di, EAvect)
    @timeit to "Vi 1.5" mat_mul!(Vi, EAvect', EAvec_Di) 
    if prev_varcmp[2] > 0
        @timeit to "Vi 1.6" mat_mul!(Ai_U, Vi, LRA_factors.U) 
        @timeit to "Vi 1.7" mat_mul!(V_Ai, LRA_factors.V, Vi) 
        @timeit to "Vi 1.8" mat_mul!(V_Ai_U, LRA_factors.V, Ai_U)
        @timeit to "Vi 1.9" V_Ai_U[diagind(V_Ai_U)] .+= 1.0 ./ (prev_varcmp[2] .* LRA_factors.D)
        @timeit to "Vi 1.10" try
            V_Ai_U .= inv(V_Ai_U) 
        catch e
            return Inf
        end
        @timeit to "Vi 1.11" mat_mul!(Ai_U_Di_V_Ai_U_i, Ai_U, V_Ai_U) 
        @timeit to "Vi 1.12" mat_mul!(Ai_U_Di_V_Ai_U_i_V_Ai, Ai_U_Di_V_Ai_U_i, V_Ai)
        @timeit to "Vi 1.13" Vi .-= Ai_U_Di_V_Ai_U_i_V_Ai 
    end
    return logdet_V
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}, EA::Eigen{T, T, Matrix{T}, Vector{T}}, cVi::Matrix{T}; is_calcu_logdet::Bool=true) where T
    if is_calcu_logdet
        @timeit to "Vi 1.0" Vi .= VCM_Arr[2]
        n = size(Vi, 1)
        @timeit to "Vi 1.1" Di = I(n) |> Matrix{Float64}
        @timeit to "Vi 1.2" mat_mul!(Vi, VCM_Arr[1], Di, prev_varcmp[1], prev_varcmp[2])
        @timeit to "Vi 1.3" begin
            cVi .= Vi
            if isposdef!(cVi)
                logdet_V = logdet(UpperTriangular(cVi)) * 2
            else
                logdet_V = logdet(Vi)
            end
        end
    else
        logdet_V = NAN
    end
    @timeit to "Vi 1.4" Di[diagind(Di)] .= 1.0 ./ (prev_varcmp[1] .* EA.values .+ prev_varcmp[2])
    @timeit to "Vi 1.5" mat_mul!(Vi, mat_mul(EA.vectors, Di), EA.vectors')
    return logdet_V
end
function calcu_Vi!(Vi::Matrix{T}, VCM_Arr::Vector{Matrix{T}}, prev_varcmp::Vector{T}, EA::Eigen{T, T, Matrix{T}, Vector{T}}, Di::Diagonal{T, Vector{T}}, EAvec_Di::Matrix{T}, EAvect::Matrix{T}, cVi::Matrix{T}; is_calcu_logdet::Bool=true) where T
    if is_calcu_logdet
        @timeit to "Vi 1.0" Vi .= VCM_Arr[1]
        n = size(Vi, 1)
        @timeit to "Vi 1.1" begin
            rmul!(Vi, prev_varcmp[1])
            Vi[diagind(Vi)] .+= prev_varcmp[2]
        end
        @timeit to "Vi 1.2" begin
            cVi .= Vi
            if isposdef!(cVi)
                logdet_V = logdet(UpperTriangular(cVi)) * 2
            else
                logdet_V, _ = logabsdet(Vi)
            end
        end
    else
        logdet_V = NAN
    end
    @timeit to "Vi 1.3" Di[diagind(Di)] .= 1.0 ./ (prev_varcmp[1] .* EA.values .+ prev_varcmp[2])
    @timeit to "Vi 1.4" mul!(EAvec_Di, Di, EAvect)
    @timeit to "Vi 1.5" mat_mul!(Vi, EA.vectors, EAvec_Di)
    return logdet_V
end
function calcu_P!(Vi::Matrix{T}, Vi_X::Matrix{T}, Xt_Vi_X_i::Matrix{T}, P::Matrix{T}, X::Matrix{T}; is_calcu_logdet::Bool=true) where T
    mat_mul!(Vi_X, Vi, X)
    mat_mul!(Xt_Vi_X_i, X', Vi_X)
    if is_calcu_logdet
        logdet_Xt_Vi_X, _ = logabsdet(Xt_Vi_X_i)
    else
        logdet_Xt_Vi_X = NAN
    end
    try
        Xt_Vi_X_i .= inv(Xt_Vi_X_i)
    catch e
        return Inf
    end
    mat_mul!(P, mat_mul(Vi_X, Xt_Vi_X_i), Vi_X')
    P .-= Vi
    lmul!(-1.0, P)
    return logdet_Xt_Vi_X
end
function calcu_P!(Vi::Matrix{T}, Vi_X::Matrix{T}, Xt_Vi_X_i::Matrix{T}, P::Matrix{T}, X::Matrix{T}, Vi_X_Xt_Vi_X_i::Matrix{T}; is_calcu_logdet::Bool=true) where T
    mat_mul!(Vi_X, Vi, X) 
    mat_mul!(Xt_Vi_X_i, X', Vi_X) 
    if is_calcu_logdet
        logdet_Xt_Vi_X, _ = logabsdet(Xt_Vi_X_i)
    else
        logdet_Xt_Vi_X = NAN
    end
    try
        Xt_Vi_X_i .= inv(Xt_Vi_X_i)
    catch e
        return Inf
    end
    mat_mul!(Vi_X_Xt_Vi_X_i, Vi_X, Xt_Vi_X_i) 
    mat_mul!(P, Vi_X_Xt_Vi_X_i, Vi_X')
    P .-= Vi
    lmul!(-1.0, P)
    return logdet_Xt_Vi_X
end
function calcu_tr_PA!(PA::Matrix{T}, P::Matrix{T}, tr_PA::Vector{T}, VCM_Arr::Vector{Matrix{T}}) where T
    _n_cmp = length(VCM_Arr)
    for i in 1:_n_cmp
        PA .= P .* VCM_Arr[i]
        tr_PA[i] = sum(PA)
    end
end
function calcu_tr_PA!(P::Matrix{T}, tr_PA::Vector{T}, VCM_Arr::Vector{Matrix{T}}) where T
    _n_cmp = length(VCM_Arr)
    PA = similar(P)
    for i in 1:_n_cmp
        PA .= P .* VCM_Arr[i]
        tr_PA[i] = sum(PA)
    end
end
function em_reml!(PA::Matrix{T}, P::Matrix{T}, Py::Vector{T}, prev_varcmp::Vector{T}, varcmp::Vector{T}, VCM_Arr::Vector{Matrix{T}}, y::Vector{T}) where T
    _n_cmp = length(varcmp)
    tr_PA = zeros(T, _n_cmp)
    _n = length(Py)
    calcu_tr_PA!(PA, P, tr_PA, VCM_Arr)  
    mat_mul!(Py, P, y)
    R = zeros(T, _n_cmp)
    @inbounds for i in 1:_n_cmp
        R[i] = dot(mat_mul(VCM_Arr[i],Py), Py)
        varcmp[i] = prev_varcmp[i] - prev_varcmp[i] * prev_varcmp[i] * (tr_PA[i] - R[i]) / _n
    end
    return R
end
function ai_reml!(PA::Matrix{T}, P::Matrix{T}, Hi::Matrix{T}, Py::Vector{T}, prev_varcmp::Vector{T}, varcmp::Vector{T}, dlogL::Float64, VCM_Arr::Vector{Matrix{T}}, y::Vector{T}; is_calcu_logdet::Bool=true, prev_prev_varcmp::Union{Nothing,Vector{T}}=nothing) where T
    _n = length(y)
    _n_cmp = length(varcmp)
    mat_mul!(Py, P, y)
    APy = zeros(T, _n, _n_cmp)
    APy[:,_n_cmp] .= Py 
    @timeit to "ai_reml Apy1" @inbounds for i in 1:(_n_cmp-1)
        APy[:, i] .= mat_mul(VCM_Arr[i], Py)
    end
    @timeit to "ai_reml Apy2"  begin
        R = zeros(T, _n_cmp)
        P_APy = zeros(T, _n, _n_cmp)
        mat_mul!(R, APy', Py)
        mat_mul!(P_APy, P, APy)
        mat_mul!(Hi, APy', P_APy)
    end
    Hi .*= 0.5
    tr_PA = zeros(T, _n_cmp)
    @timeit to "calcu_tr_PA" calcu_tr_PA!(PA, P, tr_PA, VCM_Arr)
    R .= -0.5 * (tr_PA - R)
    delta = zeros(T, _n_cmp)
    try 
        Hi .= inv(Hi)
    catch e
        return fill!(delta, NAN), fill!(R, NAN)
    end
    mul!(delta,Hi,R)
    if is_calcu_logdet
        if dlogL > 1.0 
            varcmp .= prev_varcmp + 0.316 * delta
        else
            varcmp .= prev_varcmp + delta
        end
    else
        if isnothing(prev_prev_varcmp)
            error("'prev_prev_varcmp' must be provided!")
        end
        if norm(varcmp - prev_prev_varcmp)^2 / norm(varcmp)^2 > 0.01
            varcmp .= prev_varcmp + 0.316 * delta
        else
            varcmp .= prev_varcmp + delta
        end
    end
    return delta, R
end
function constrain_varcmp!(varcmp::Vector{T},y::Vector{T}) where T
    _n_cmp = length(varcmp)
    delta = 0.0
    constr_scale = 1e-6
    num = 0
    constrain = zeros(Int, _n_cmp)
    for i in 1:_n_cmp
        if varcmp[i] < 0
            Var_y = var(y)
            delta += Var_y * constr_scale - varcmp[i]
            varcmp[i] = Var_y * constr_scale
            constrain[i] = 1
            num += 1
        end
    end
    delta /= _n_cmp - num
    for i in 1:_n_cmp
        if constrain[i] < 1 && varcmp[i] > delta
            varcmp[i] -= delta
        end
    end
    return num
end
function cg_inv(M::Matrix{T}) where T
    n = size(M,1)
    @timeit to "cg_inv 1.0" Mi = zeros(T, n,n)
    @timeit to "cg_inv 1.1" b = I(n) |> Matrix{Float64}
    @timeit to "cg_inv 1.2" @inbounds Threads.@threads for i in 1:n
        Mi[:, i] .= cg(M, b[:, i])
    end
    return Mi
end
function cg_inv!(M::Matrix{T}, Mi::Matrix{T}) where T
    n = size(M,1)
    b = I(n) |> Matrix{Float64}
    @inbounds Threads.@threads for i in 1:n
        Mi[:, i] .= cg(M, b[:, i])
    end
    return true
end
function calcu_Vp(varcmp::Vector{T}, Hi::Matrix{T}) where T
    _n_cmp = length(varcmp)
    Vp = 0.0
    VarVp = 0.0
    for i in 1:_n_cmp
        Vp += varcmp[i]
        for j in 1:_n_cmp
            VarVp += Hi[i, j]
        end
    end
    return Vp, VarVp
end
function calcu_hsq(i::Int, Vp::Float64, VarVp::Float64, varcmp::Vector{T}, Hi::Matrix{T}) where T
    _n_cmp = length(varcmp)
    V1 = varcmp[i]
    VarV1 = Hi[i, i]
    Cov12 = 0.0
    for j in 1:_n_cmp
        Cov12 += Hi[i, j]
    end
    hsq = V1 / Vp
    var_hsq = (V1 / Vp) * (V1 / Vp) * (VarV1 / (V1 * V1) + VarVp / (Vp * Vp) - (2 * Cov12) / (V1 * Vp))
    return hsq, var_hsq
end
function getOpenMendelVC(G::Matrix{T}, Y::Vector{T}; EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing, X::Union{Nothing,Matrix{T}}=nothing, method::String=["MLE","REML"][2], algo::Symbol=[:MM,:FS][1]) where T
    n = length(Y)
    if isnothing(X)
        X = ones(T, n, 1)
    end
    vcdata = VarianceComponentVariate(Y, X, (G, spdiagm(ones(T, n))))
    vcmodel = VarianceComponentModel(vcdata)
    if !isnothing(EA)
        if method == "MLE"
            logl, vcmodel, Σse, Σcov, Bse, Bcov = fit_mle!(vcmodel, vcdata, EA; algo = algo)
        else
            logl, vcmodel, Σse, Σcov, Bse, Bcov = fit_reml!(vcmodel, vcdata, EA; algo = algo)
        end
    else
        if method == "MLE"
            logl, vcmodel, Σse, Σcov, Bse, Bcov = fit_mle!(vcmodel, vcdata; algo = algo)
        else
            logl, vcmodel, Σse, Σcov, Bse, Bcov = fit_reml!(vcmodel, vcdata; algo = algo)
        end
    end
    h, hse = VarianceComponentModels.heritability(vcmodel.Σ,Σcov)
    abc_uni = Dict(
            :ΣG => vec(vcmodel.Σ[1]),
            :ΣG_se => vec(Σse[1]),
            :Σe => vcmodel.Σ[2][1],
            :Σe_se => Σse[2][1],
            :Σp => nothing,
            :Σp_se => nothing,
            :h2 => h,
            :h2_se => hse,
            :logL => logl,
            :logL0 => nothing,
            :n => n,
            :Fix_eff => vec(vcmodel.B),
            :Fix_eff_se => vec(Bse),
            :u => nothing
            )
    return abc_uni
end
function getMRVCModel(multi_KSs::Vector{Matrix{T}}, Y::Vector{T}; X::Union{Nothing,Matrix{T}}=nothing, se::Bool=true, reml::Bool=true, maxiter::Int=100,verbose::Bool=false) where T
    n = length(Y)
    if isnothing(X)
        X = ones(T, n, 1)
    end
    n_grm = length(multi_KSs) - 1
    if n_grm == 1
        vcmodel = MRTVCModel(Y, X, multi_KSs; se=se, reml=reml)
    else
        vcmodel = MRVCModel(Y, X, multi_KSs; se=se, reml=reml)
    end
    mrvc_model_status = MultiResponseVarianceComponentModels.fit!(vcmodel; maxiter=maxiter,verbose=verbose) 
    h2s, ses, Σs, Σses = MultiResponseVarianceComponentModels.h2(vcmodel)
    Bse = sqrt.(diag(vcmodel.Bcov))
    abc_uni = Dict(
            :ΣG => Σs[1:n_grm],
            :ΣG_se => Σses[1:n_grm],
            :Σe => Σs[end],
            :Σe_se => Σses[end],
            :Σp => nothing,
            :Σp_se => nothing,
            :h2 => h2s[1:n_grm],
            :h2_se => ses[1:n_grm],
            :logL => vcmodel.logl,
            :logL0 => nothing,
            :n => n,
            :Fix_eff => vec(vcmodel.B),
            :Fix_eff_se => vec(Bse),
            :u => nothing,
            :converged => mrvc_model_status.isconverged
            )
    return abc_uni
end
function AIREML(A::AbstractMatrix,y::AbstractVector;X::Union{Nothing,Matrix{T}}=nothing,EA::Union{Nothing,Eigen{T, T, Matrix{T}, Vector{T}}}=nothing,max_inter=100,conv_threshold=1e-4) where T
    _time_strat = time()
    n = length(y)
    if isnothing(X)
        X = ones(n,1)
    end
    y_t = y'
    X_t = X' |> Matrix
    α = 1
    σ_y = var(y)
    σ_1 = σ_y / 2
    σ_2 = σ_y / 2
    Σ_i = [σ_1,σ_2]
    II = spdiagm(ones(n))
    iter = 0
    logL_i = 0
    AI = zeros(T, 2,2)
    DL = zeros(T, 2,1)
    println("Running AI-REML algorithm ...")
    println(string("Iter.","\t","logL","\t","V(G)","\t","V(e)"))
    if isnothing(EA)
        EA = eigen(A)
    end
    while true
        iter += 1
        V = Σ_i[1] .* A .+ Σ_i[2] .* II
        Dinv = spdiagm(1.0 ./ (Σ_i[1] .* EA.values .+ Σ_i[2]))
        Vinv = EA.vectors * Dinv * EA.vectors'
        X_t_Vinv = X_t * Vinv
        X_t_Vinv_X = X_t_Vinv * X
        X_t_Vinv_X_inv = inv(X_t_Vinv_X)
        P = Vinv .- Vinv*X*X_t_Vinv_X_inv*X_t_Vinv
        P_y = P*y
        logL_i_1 = -0.5*(n*log(2*π)+y'*Vinv*y+logdet(V))
        go = logL_i_1 - logL_i
        logL_i = logL_i_1
        P_P_y = P*P_y
        P_A=P*A
        y_t_P_A = y_t*P_A
        y_t_P_A_P_P_y = y_t_P_A*P_P_y
        AI[1,1] = y_t_P_A*P_A*P_y
        AI[1,2] = y_t_P_A_P_P_y
        AI[2,1] = y_t_P_A_P_P_y
        AI[2,2] = y_t*P*P_P_y
        DL[1,1] = tr(P_A) .- y_t_P_A*P_y
        DL[2,1] = tr(P) .- y_t*P_P_y
        DL = -0.5 * DL
        AI = 0.5 * AI
        AIinv_DL = inv(AI)*DL
        Σ_i_1 = Σ_i .+ α .* AIinv_DL
        while sum(Σ_i_1 .< 0) > 0
            α = 0.5 * α
            Σ_i_1 = Σ_i .+ α .* AIinv_DL
        end
        α = 1
        Σ_0 = (Σ_i_1 ./ sum(Σ_i_1) .< 1e-4) .& (Σ_i ./ sum(Σ_i) .< 1e-4)
        if sum(Σ_0) > 0
            Σ_i_1[Σ_0] .= 0
            go = 0
        end
        Σ_i = Σ_i_1
        println(string(iter,"\t",logL_i_1,"\t",Σ_i[1],"\t",Σ_i[2]))
        if (abs(go) <= conv_threshold) | (iter >= max_inter)
            println("REML converged.\n")
            var_cov = inv(AI)*2 
            se_σ_1 = sqrt(var_cov[1,1])
            se_σ_2 = sqrt(var_cov[2,2])
            println("Summary result of REML analysis:")
            println(string("Source","\t","Variance","\t","SE"))
            println(string("V(G)","\t",Σ_i[1],"\t",se_σ_1))
            println(string("V(e)","\t",Σ_i[2],"\t",se_σ_2))
            println(string("Vp","\t",sum(Σ_i),"\t","-"))
            println(string("V(G)/Vp","\t",Σ_i[1]/sum(Σ_i),"\t","-"))
            println("\nSampling variance/covariance of the estimates of variance components:")
            println(string(var_cov[1,1],"\t",var_cov[1,2]))
            println(string(var_cov[2,1],"\t",var_cov[2,2]))
            Dinv = spdiagm(1.0 ./ (Σ_i[1] .* EA.values .+ Σ_i[2]))
            Vinv = EA.vectors * Dinv *EA.vectors'
            b_bat = inv(X' * Vinv * X) * X' * Vinv * y
            u_bat = A * Σ_i[1] * Vinv * (y .- X * b_bat)
            println("\nEstimatesof fixed effects:")
            println(string("Source","\t","Estimate","\t","SE"))
            println(string("mean","\t",b_bat[1],"\t","-"))
            if length(b_bat) > 1
                for i_b_bat in 2:length(b_bat)
                    println(string("X",i_b_bat,"\t",b_bat[i_b_bat],"\t","-"))
                end
            end
            println(string("\nAnalysis finished at ",Dates.now()))
            println(string("Computational time: ",time()-_time_strat," second(s)."))
            abc_uni = Dict(
            :ΣG => Σ_i[1],
            :ΣG_se => nothing,
            :Σe => Σ_i[2],
            :Σe_se => nothing,
            :Σp => sum(Σ_i),
            :Σp_se => nothing,
            :ΣG╱Σp => Σ_i[1]/sum(Σ_i),
            :ΣG╱Σp_se => nothing,
            :logL => logL_i,
            :logL0 => nothing,
            :n => size(A,1),
            :Fix_eff => b_bat,
            :Fix_eff_se => nothing,
            :u_bat => u_bat
            )
            return abc_uni
        end
    end
end
function AIREML(RELs::Vector{<:AbstractArray},y::AbstractVector;X=missing,max_inter=100,conv_threshold=1e-4,is_rep=false)
    _time_strat = time()
    n_REL = length(RELs)
    n = length(y)
    y_t = y'
    if ismissing(X)
        X = ones(length(y))
    end
    X_t = X'
    α = 1
    iter = 0
    logL_i = 0
    n_Σ = n_REL+1
    Σ_i = repeat([std(y)^2/2/(n_REL+1)],n_Σ)
    AI = zeros(T, n_Σ,n_Σ)
    DL = zeros(T, n_Σ,1)
    II = Matrix{Float64}(I(n))
    println("Running AI-REML algorithm ...")
    println(string("Iter.","\t","logL","\t",string([[string("V(G",i,")\t") for i in 1:n_REL]...]...),"V(e)"))
    while true
        iter = iter + 1
        V = II .* Σ_i[n_REL+1]
        @inbounds for i in 1:n_REL
            BLAS.axpby!(Σ_i[i],RELs[i],1,V)
        end
        Vinv = inv(V)
        X_t_Vinv = X_t * Vinv
        X_t_Vinv_X = X_t_Vinv * X
        X_t_Vinv_X_inv = zeros(T, size(X_t,1),size(X_t,1))
        if !is_rep
            try
                X_t_Vinv_X_inv = inv(X_t_Vinv_X)
            catch e
                X_t_Vinv_X_inv = pinv(X_t_Vinv_X)
            end
        else
            X_t_Vinv_X_inv = pinv(X_t_Vinv_X)
        end
        P = Vinv .- Vinv*X*X_t_Vinv_X_inv*X_t_Vinv
        P_y = P*y
        logL_i_1 = -0.5*(n*log(2*π)+y'*Vinv*y+logdet(V))::Float64
        go = logL_i_1 - logL_i
        logL_i = logL_i_1
        y_t_P = y_t*P  
        @inbounds y_t_P_ = [[y_t_P*RELs[i] for i in 1:n_REL]...,y_t_P]  
        @inbounds _P_y = [[RELs[i]*P_y for i in 1:n_REL]...,P_y]  
        @inbounds for AI_row_i in 1:n_Σ
            @inbounds for AI_col_i in AI_row_i:n_Σ
                AI[AI_row_i,AI_col_i] = (y_t_P_[AI_col_i]*P*_P_y[AI_row_i])[1,1]
                AI[AI_col_i,AI_row_i] = AI[AI_row_i,AI_col_i]
            end
        end
        @inbounds P_ = [P*RELs[i] for i in 1:n_REL]
        @inbounds tr_s = [[tr(P_[i]) for i in 1:n_REL]...,tr(P)]
        @inbounds DL_r = [(y_t_P_[i] * P_y)[1,1] for i in 1:n_Σ]
        DL .= -(tr_s .- DL_r)
        AIinv_DL = pinv(AI)*DL
        Σ_i_1 = Σ_i .+ α .* AIinv_DL
        while sum(Σ_i_1 .< 0) > 0
            α = 0.5 * α
            Σ_i_1 = Σ_i .+ α .* AIinv_DL
        end
        α = 1
        Σ_0 = (Σ_i_1 ./ sum(Σ_i_1) .< 1e-4) .& (Σ_i ./ sum(Σ_i) .< 1e-4)
        if sum(Σ_0) > 0
            Σ_i_1[Σ_0] .= 0
        end
        Σ_i = Σ_i_1
        println(string(iter,"\t",logL_i_1,string([[string("\t",Σ_i[i]) for i in 1:n_REL]...]...,"\t",Σ_i[end])))
        if (abs(go) <= conv_threshold) | (iter >= max_inter)
            println("REML converged.\n")
            var_cov = pinv(AI)*2 
            se_σ_1 = sqrt(var_cov[1,1])
            se_σ_2 = sqrt(var_cov[2,2])
            println("Summary result of REML analysis:")
            println(string("Source","\t","Variance","\t","SE"))
            print([string("V(G",i,")","\t",Σ_i[i],"\t","-","\n") for i in 1:n_REL]...)
            println(string("V(e)","\t",Σ_i[end],"\t","-"))
            println(string("Vp","\t",sum(Σ_i),"\t","-"))
            println(string("V(G)/Vp","\t",sum(Σ_i[1:n_REL])/sum(Σ_i),"\t","-"))
            println("\nSampling variance/covariance of the estimates of variance components:")
            b_bat = [nothing]
            println("\nEstimatesof fixed effects:")
            println(string("Source","\t","Estimate","\t","SE"))
            println(string("mean","\t",b_bat[1],"\t","-"))
            if length(b_bat) > 1
                for i_b_bat in 2:length(b_bat)
                    println(string("X",i_b_bat,"\t",b_bat[i_b_bat],"\t","-"))
                end
            end
            println(string("\nAnalysis finished at ",Dates.now()))
            println(string("Computational time: ",time()-_time_strat," second(s)."))
            abc_uni = Dict(
            :ΣG => Σ_i[1:n_REL],
            :ΣG_se => nothing,
            :Σe => Σ_i[end],
            :Σe_se => nothing,
            :Σp => sum(Σ_i),
            :Σp_se => nothing,
            :ΣG╱Σp => Σ_i[1:n_REL] ./ sum(Σ_i),
            :ΣG╱Σp_se => nothing,
            :logL => logL_i,
            :logL0 => nothing,
            :n => n,
            :Fix_eff => b_bat,
            :Fix_eff_se => nothing
            )
            return abc_uni
        end
    end
end
function reml_GCTA(REL1::AbstractMatrix, y0::AbstractVector; work_dir::String=pwd(), gcta_path::String="gcta64", Xfixed=missing, reml_alg::Int=0, n_threads::Int=22, is_debug::Bool=false)
    phen = DataFrame(family_ID = 1:length(y0), individual_ID = 1:length(y0), phenotypes = y0)
    wk_dir_tmp = mktempdir(work_dir;prefix="gcta_",cleanup=false)
    cd(wk_dir_tmp)
    CSV.write("gcta.phen",phen,header=false,delim=' ')
    idx_ = string.(1:size(REL1,1))
    writeGRMBin(REL1, idx_, "REL_tmp")
    write("multi_grm.txt",string("REL_tmp"))
    if ~ismissing(Xfixed)
        Xfixed = Xfixed[:,1] == ones(size(Xfixed,1)) ? Xfixed[:,2:end] : Xfixed
        covar = [DataFrame(fam_id=idx_, ind_id=idx_) DataFrame(Xfixed,:auto)]
        CSV.write("gcta.covar", covar, delim=" ", header=false)
        code = `$gcta_path --reml --mgrm multi_grm.txt --qcovar gcta.covar --reml-alg $reml_alg --pheno gcta.phen --thread-num $n_threads --reml-est-fix --out gcta`
    else
        code = `$gcta_path --reml --mgrm multi_grm.txt --reml-alg $reml_alg --pheno gcta.phen --thread-num $n_threads --reml-est-fix --out gcta`
    end
    run(code)
    hsq = CSV.read("gcta.hsq",DataFrame,silencewarnings=true)
    cd("..")
    if !is_debug
        rm(wk_dir_tmp, recursive=true, force=true)
    end
    return hsq
end
function reml_GCTA(RELs::Vector{<:AbstractArray}, y0::AbstractVector; work_dir::String=pwd(), gcta_path::String="gcta64", Xfixed=missing, reml_alg::Int=0, n_threads::Int=22, is_debug::Bool=false)
    n_RELs = length(RELs)
    phen = DataFrame(family_ID = 1:length(y0), individual_ID = 1:length(y0), phenotypes = y0)
    wk_dir_tmp = mktempdir(work_dir;prefix="gcta_",cleanup=false)
    cd(wk_dir_tmp)
    CSV.write("gcta.phen",phen,header=false,delim=' ')
    idx_ = string.(1:size(RELs[1],1))
    for i in 1:n_RELs
        writeGRMBin(RELs[i], idx_, string("REL",i,"_tmp"))
    end
    multi_grm_s = [string("REL",i,"_tmp\n") for i in 1:n_RELs]
    multi_grm_txt = ""
    for i in 1:length(multi_grm_s)
        multi_grm_txt = string(multi_grm_txt, multi_grm_s[i])
    end
    write("multi_grm.txt",multi_grm_txt)
    if ~ismissing(Xfixed)
        Xfixed = Xfixed[:,1] == ones(size(Xfixed,1)) ? Xfixed[:,2:end] : Xfixed
        covar = DataFrame(hcat(idx_, idx_, Xfixed),:auto)
        CSV.write("gcta.covar", covar, delim=" ", header=false)
        code = `$gcta_path --reml --mgrm multi_grm.txt --covar gcta.covar --reml-alg $reml_alg --pheno gcta.phen --thread-num $n_threads --reml-est-fix --out gcta`
    else
        code = `$gcta_path --reml --mgrm multi_grm.txt --reml-alg $reml_alg --pheno gcta.phen --thread-num $n_threads --reml-est-fix --out gcta`
    end
    run(code)
    hsq = CSV.read("gcta.hsq",DataFrame,silencewarnings=true)
    cd("..")
    if !is_debug
        rm(wk_dir_tmp, recursive=true, force=true)
    end
    return hsq
end
function reml_LDAK(REL1::AbstractMatrix, y0::AbstractVector; work_dir::String=pwd(), ldak_path::String="ldak5.linux.fast", Xfixed=missing, is_debug::Bool=false)
    phen = DataFrame(family_ID = 1:length(y0), individual_ID = 1:length(y0), phenotypes = y0)
    wk_dir_tmp = mktempdir(work_dir;prefix="ldak_",cleanup=false)
    cd(wk_dir_tmp)
    CSV.write("ldak.phen",phen,header=false,delim=' ')
    idx_ = string.(Array(1:size(REL1,1)))
    writeGRMBin(REL1, idx_, "REL_tmp")
    write("multi_grm.txt",string("REL_tmp"))
    if ~ismissing(Xfixed)
        Xfixed = Xfixed[:,1] == ones(size(Xfixed,1)) ? Xfixed[:,2:end] : Xfixed
        covar = DataFrame(hcat(idx_, idx_, Xfixed),:auto)
        CSV.write("ldak.covar", covar, delim=" ", header=false)
        code = `$ldak_path --reml multiKernel --pheno ldak.phen --mgrm multi_grm.txt --covar ldak.covar --permute NO --kinship-details NO --constrain YES`
    else
        code = `$ldak_path --reml multiKernel --pheno ldak.phen --mgrm multi_grm.txt --permute NO --kinship-details NO --constrain YES`
    end
    run(code)
    vars = CSV.read("multiKernel.vars",DataFrame,silencewarnings=true)
    cd("..")
    if !is_debug
        rm(wk_dir_tmp, recursive=true, force=true)
    end
    return vars
end
function reml_LDAK(RELs::Vector{<:AbstractMatrix}, y0::AbstractVector; work_dir::String=pwd(), ldak_path::String="ldak5.linux.fast", Xfixed=missing, is_debug::Bool=false)
    n_RELs = length(RELs)
    phen = DataFrame(family_ID = 1:length(y0), individual_ID = 1:length(y0), phenotypes = y0)
    wk_dir_tmp = mktempdir(work_dir;prefix="ldak_",cleanup=false)
    cd(wk_dir_tmp)
    CSV.write("ldak.phen",phen,header=false,delim=' ')
    idx_ = string.(1:size(RELs[1],1))
    for i in 1:n_RELs
        writeGRMBin(RELs[i], idx_, string("REL",i,"_tmp"))
    end
    multi_grm_s = [string("REL",i,"_tmp\n") for i in 1:n_RELs]
    multi_grm_txt = ""
    for i in 1:length(multi_grm_s)
        multi_grm_txt = string(multi_grm_txt, multi_grm_s[i])
    end
    write("multi_grm.txt",multi_grm_txt)
    if ~ismissing(Xfixed)
        Xfixed = Xfixed[:,1] == ones(size(Xfixed,1)) ? Xfixed[:,2:end] : Xfixed
        covar = DataFrame(hcat(idx_, idx_, Xfixed),:auto)
        CSV.write("ldak.covar", covar, delim=" ", header=false)
        code = `$ldak_path --reml multiKernel --pheno ldak.phen --mgrm multi_grm.txt --covar ldak.covar --permute NO --kinship-details NO --constrain YES`
    else
        code = `$ldak_path --reml multiKernel --pheno ldak.phen --mgrm multi_grm.txt --permute NO --kinship-details NO --constrain YES`
    end
    run(code)
    vars = CSV.read("multiKernel.vars",DataFrame,silencewarnings=true)
    cd("..")
    if !is_debug
        rm(wk_dir_tmp, recursive=true, force=true)
    end
    return vars
end
function getVarianceComponents(REL1::AbstractMatrix, y0::AbstractVector; work_dir=pwd(), eVC_method::String=["JuliaPkg","GCTA","LDAK","AIREML"][1], gcta_path::String="", ldak_path::String="", Xfixed=missing, gcta_reml_alg::Int=0, n_threads::Int=2, is_debug::Bool=false)
    refer_idx = .~isnan.(y0)
    REL1_ = REL1[refer_idx,refer_idx]
    y0_ = Vector{Float64}(y0[refer_idx])
    if ~ismissing(Xfixed)
        Xfixed_ = Xfixed[refer_idx,:]
    else
        Xfixed_ = Xfixed
    end
    if eVC_method == "JuliaPkg"
        abc = getVarianceComponent(REL1_, y0_; X=Xfixed_)
        abc_uni = Dict(
        :ΣG => abc[:vcmodel].Σ[1][1],
        :ΣG_se => abc[:Σse][1][1],
        :Σe => abc[:vcmodel].Σ[2][1],
        :Σe_se => abc[:Σse][2][1],
        :Σp => abc[:vcmodel].Σ[1][1]+abc[:vcmodel].Σ[2][1],
        :Σp_se => NaN,
        :ΣG╱Σp => abc[:h][1],
        :ΣG╱Σp_se => abc[:hse][1],
        :logL => abc[:logl],
        :logL0 => NaN,
        :n => length(y0_),
        :Fix_eff => abc[:vcmodel].B,
        :Fix_eff_se => abc[:Bse]
        )
    elseif eVC_method == "GCTA"
        abc = reml_GCTA(REL1_, y0_; Xfixed=Xfixed_, gcta_path=gcta_path, reml_alg=gcta_reml_alg, work_dir=work_dir, n_threads=n_threads, is_debug=is_debug)
        abc_uni = Dict(
        :ΣG => parse(Float64,abc[1,2]),
        :ΣG_se => abc[1,3],
        :Σe => parse(Float64,abc[2,2]),
        :Σe_se => abc[2,3],
        :Σp => parse(Float64,abc[3,2]),
        :Σp_se => abc[3,3],
        :ΣG╱Σp => parse(Float64,abc[4,2]),
        :ΣG╱Σp_se => abc[4,3],
        :logL => parse(Float64,abc[5,2]),
        :logL0 => parse(Float64,abc[6,2]),
        :n => parse(Int,abc[10,2]),
        :Fix_eff => parse.(Float64,abc[13:end,1]),
        :Fix_eff_se => parse.(Float64,abc[13:end,2])
        )
    elseif eVC_method == "LDAK"
        abc = reml_LDAK(REL1_, y0_; Xfixed= Xfixed_, work_dir=work_dir, ldak_path=ldak_path, is_debug=is_debug)
        abc_uni = Dict(
        :ΣG => abc[1,2],
        :ΣG_se => nothing,
        :Σe => abc[2,2],
        :Σe_se => nothing,
        :Σp => nothing,
        :Σp_se => nothing,
        :ΣG╱Σp => nothing,
        :ΣG╱Σp_se => nothing,
        :logL => nothing,
        :logL0 => nothing,
        :n => nothing,
        :Fix_eff => nothing,
        :Fix_eff_se => nothing
        )
    elseif eVC_method == "AIREML"
        abc_uni = AIREML(REL1_,y0_;X=Xfixed_)
    end
    return abc_uni
end
function getVarianceComponents(RELs::Vector{<:AbstractArray}, y0::AbstractVector; work_dir=pwd(), eVC_method::String=["JuliaPkg","GCTA","LDAK","AIREML"][1], gcta_path::String="", ldak_path::String="", Xfixed=missing, gcta_reml_alg::Int=0, n_threads::Int=2, is_debug::Bool=false)
    refer_idx = .~isnan.(y0)
    RELs_ = @views [RELs[i][refer_idx,refer_idx] for i in 1:length(RELs)]
    n_RELs = length(RELs_)
    y0_ = @view y0[refer_idx]
    if ~ismissing(Xfixed)
        Xfixed_ = @view Xfixed[refer_idx,:]
    else
        Xfixed_ = Xfixed
    end
    if eVC_method == "GCTA"
        abc = reml_GCTA(RELs_, y0_; Xfixed= Xfixed_, gcta_path=gcta_path, reml_alg=gcta_reml_alg, work_dir=work_dir, n_threads=n_threads, is_debug=is_debug)
        abc_uni = Dict(
        :ΣG => parse.(Float64,abc[1:n_RELs,2]),
        :ΣG_se => Float64.(abc[1:n_RELs,3]),
        :Σe => parse(Float64,abc[n_RELs+1,2]),
        :Σe_se => Float64.(abc[n_RELs+1,3]),
        :Σp => parse(Float64,abc[n_RELs+2,2]),
        :Σp_se => Float64.(abc[n_RELs+2,3]),
        :ΣG╱Σp => parse.(Float64,abc[(n_RELs+3):(n_RELs+n_RELs+2),2]),
        :ΣG╱Σp_se => Float64.(abc[(n_RELs+3):(n_RELs+n_RELs+2),3]),
        :logL => parse(Float64,abc[n_RELs+n_RELs+5,2]),
        :logL0 => parse(Float64,abc[n_RELs+n_RELs+6,2]),
        :n => parse(Int,abc[n_RELs+n_RELs+10,2]),
        :Fix_eff => parse.(Float64,abc[(n_RELs+n_RELs+13):end,1]),
        :Fix_eff_se => parse.(Float64,abc[(n_RELs+n_RELs+13):end,2])
        )
    elseif eVC_method == "LDAK"
        abc = reml_LDAK(RELs_, y0_; Xfixed= Xfixed_, work_dir=work_dir, ldak_path=ldak_path, is_debug=is_debug)
        abc_uni = Dict(
        :ΣG => abc[1:n_RELs,2],
        :ΣG_se => nothing,
        :Σe => abc[end,2],
        :Σe_se => nothing,
        :Σp => nothing,
        :Σp_se => nothing,
        :ΣG╱Σp => nothing,
        :ΣG╱Σp_se => nothing,
        :logL => nothing,
        :logL0 => nothing,
        :n => nothing,
        :Fix_eff => nothing,
        :Fix_eff_se => nothing
        )
    elseif eVC_method == "AIREML"
        abc_uni = AIREML(RELs_,y0_;X=Xfixed_)
    end
    return abc_uni
end
