function idul_data(eG, exppheno, X_MME)
    D = copy(eG.values)
    D[D.<0] .= 0
    _1_D = [ones(length(D)) D]
    xQ = eG.vectors' * X_MME 
    yQ = eG.vectors' * exppheno 
    return xQ, yQ, D, _1_D
end
function get_lm_residuals(X_MME::Matrix{T}, y::Vector{T}, weights::Vector{T}) where {T}
    W = spdiagm(weights)
    beta0 = (X_MME' * W * X_MME) \ (X_MME' * W * y)
    yresid = y - X_MME * beta0
    return yresid, beta0
end
function get_lm_residuals(X_MME::Matrix{T}, y::Matrix{T}, weights::Vector{T}) where {T}
    W = spdiagm(weights)
    beta0 = (X_MME' * W * X_MME) \ (X_MME' * W * y)
    yresid = y - X_MME * beta0
    return yresid, beta0
end
function get_lm_coef(X_MME::Matrix{T}, y::Vector{T}, weights::Vector{T}) where {T}
    W = spdiagm(weights)
    beta0 = (X_MME' * W * X_MME) \ (X_MME' * W * y)
    return beta0
end
function get_lm_coef(X_MME::Matrix{T}, y::Matrix{T}, weights::Vector{T}) where {T}
    W = spdiagm(weights)
    beta0 = (X_MME' * W * X_MME) \ (X_MME' * W * y)
    return beta0
end
function get_IDUL_VarianceComponent(EA::Eigen{T,T,Matrix{T},Vector{T}}, exppheno::Vector{T}, X_MME::Matrix{T}; xQ::Union{Nothing,Matrix{T}}=nothing, max_iter::Int=100, thre::T=T(1e-6), algo::String="reml") where {T}
    if isnothing(xQ)
        xQ = EA.vectors' * X_MME
    end
    yQ = EA.vectors' * exppheno
    _, h2, vc_B = idul_h(xQ, yQ, EA.values; max_iter=max_iter, thre=thre, algo=algo)
    yadj = exppheno - X_MME * vc_B
    Var_y = var(yadj)
    NaN_T = T.(NaN)
    Σ_i = [NaN_T, NaN_T]
    Σ_i[1] = h2 * Var_y
    Σ_i[2] = Var_y - Σ_i[1]
    abc_uni = Dict(
        :ΣG => [Σ_i[1]],
        :ΣG_se => [NaN_T],
        :Σe => Σ_i[2],
        :Σe_se => NaN_T,
        :Σp => Var_y,
        :Σp_se => NaN_T,
        :h2 => [h2],
        :h2_se => [NaN_T],
        :logL => NaN_T,
        :logL0 => NaN_T,
        :n => NaN_T,
        :Fix_eff => vc_B,
        :Fix_eff_se => nothing
    )
    return abc_uni
end
function idul_h(xQ::Matrix{T}, yQ::Vector{T}, evals::Vector{T}; max_iter::Int=100, thre=1e-6, algo::String="mle") where {T}
    h = T.(0.1) 
    eta0 = h / (1 - h)
    eta1 = copy(eta0)
    steps = 0
    diff1 = 1
    n_samples = length(yQ)
    _1_D = [ones(T, n_samples) evals]
    coefs = zeros(size(xQ, 2))
    while (diff1 > thre) & (steps < max_iter)
        steps += 1
        H = eta0 .* evals .+ 1
        w1 = 1 ./ H
        w2 = abs2.(w1)
        resid, coefs = get_lm_residuals(xQ, yQ, w1)
        r2 = abs2.(resid)
        fit_coef = get_lm_coef(_1_D, r2, w2)
        if algo == "mle"
            eta1 = fit_coef[2] / fit_coef[1]
        elseif algo == "reml"
            ta1 = fit_coef[2] / fit_coef[1]
            ta2 = fit_coef[2] / Statistics.mean(r2 ./ H)
            eta1 = (ta1 + ta2) / 2
        end
        if eta1 >= 1e5
            eta1 = T.(1e5)
        elseif eta1 < 1e-5
            eta1 = T.(1e-5)
        end
        diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
        eta0 = copy(eta1)
    end 
    h = eta1 / (1 + eta1)
    return eta1, h, coefs
end
function idul_plus_assoc_test!(XS0::Matrix{T}, y::Vector{T}, X_MME::Matrix{T}, evals::Vector{T}; test_method::String="wald", is_calcu_pv::Bool=true) where {T}
    test_method = "wald"
    max_iter = 100
    thre = 1e-6
    y = _struct_KIN.EA.vectors' * _struct_PHENO.expression[1, :] |> Vector{FloatT}
    _2p = 2 .* _struct_GENO.snp_annotation.af
    _2pq = _2p .* (1 .- _struct_GENO.snp_annotation.af)
    _2pq = FloatT.(_2pq)
    _2pq = _2pq' |> Matrix
    XS0 = _struct_GENO.genotype .- _2pq 
    XS0 .= _struct_KIN.EA.vectors' * XS0
    evals = FloatT.(_struct_KIN.EA.values)
    evals[evals.<0] .= 0
    X_MME = _struct_COVAR.X_MME
    xQ = _struct_KIN.EA.vectors' * X_MME
    xQ = FloatT.(xQ)
    _ni, _ns = size(XS0)
    _c = size(xQ, 2)
    _c1 = _c + 1
    xW0 = [ones(FloatT, _ni) xQ]
    d = WaldTest(1, _ni - _c1)
    df_test = zeros(FloatT, _ns, 5)
    @notimeit to @inbounds @threads for xi in 1:_ns
        if xi % 1000 == 0
            println(xi)
        end
        xW = copy(xW0)
        xW[:, 1] .= XS0[:, xi] 
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            xtxi = zeros(FloatT, _c1, _c1)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
            mat = zeros(FloatT, _c1, _c1)
        end
        steps = 0
        eta0 = 1
        eta1 = 0.1f0
        rh1r = 0.0f0
        begin
            eta0 = copy(eta1)
            hh = eta0 .* evals .+ 1
            h2inv .= hh 
            h2inv .\= 1 
            h1inv .= sqrt.(h2inv)
            X .= xW .* h1inv
            mat_mul!(xtx, X', X)
            xtx[diagind(xtx)] .+= 1.0f-10
            y1 .= y .* h1inv
            mat_mul!(xty, X', y1)
            begin 
                try
                    xtxi .= inv(cholesky!(Symmetric(xtx)))
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtxi, xty)
            end
            r2 .= y1 .- X * beta1 
            r2 .*= r2 
            tauinv = Statistics.mean(r2) 
        end
        tauinv1 = 0
        like_eta0 = 0
        like_eta1 = 0
        if test_method == "wald" 
            like_eta0 = -log(tauinv) * (_ni - _c)
            like_eta0 -= sum(log.(hh))
            like_eta0 -= logdet(UpperTriangular(xtx)) * 2
            while true
                rh1r = sum(r2)
                rh2r = dot(r2, h2inv)
                trh1 = sum(h2inv)
                trh2 = dot(h2inv, h2inv)
                mh = trh1 / _ni
                vh = trh2 / _ni - mh * mh
                X .= xW .* h2inv
                mat_mul!(xtx, X', X) 
                mul!(mat, xtxi, xtx) 
                fp1 = trh1 - tr(mat) - (_ni - _c1) * rh2r / rh1r
                eta1 = eta0 + fp1 * eta0 / _ni / vh * (0.5f0 + 0.5f0 * mh)
                eta1 = clamp(eta1, 1.0f-5, 1.0f5)
                while true
                    if abs(eta1 / (1 + eta1) - eta0 / (1 + eta0)) < 1.0f-6
                        break
                    end
                    steps += 1
                    hh .= eta1 .* evals .+ 1
                    h2inv .= 1 ./ hh
                    h1inv .= sqrt.(h2inv)
                    X .= xW .* h1inv
                    mul!(xtx, X', X)
                    xtx[diagind(xtx)] .+= 1.0f-10
                    y1 .= y .* h1inv
                    mul!(xty, X', y1)
                    beta1 .= xtx \ xty
                    e1 = y1 - X * beta1
                    r2 .= abs2.(e1)
                    tauinv1 = Statistics.mean(r2)
                    like_eta1 = -log(tauinv1) * (_ni - _c1)
                    like_eta1 -= sum(log.(hh))
                    like_eta1 -= logdet(xtx)
                    if like_eta1 >= like_eta0
                        break
                    end
                    eta1 = (eta1 + eta0) / 2
                end
                if abs(eta1 / (1 + eta1) - eta0 / (1 + eta0)) < 1.0f-6
                    break
                end
                eta0 = copy(eta1)
                tauinv = copy(tauinv1)
                like_eta0 = copy(like_eta1)
            end
        end
        if test_method == "lrt"
            like_eta0 = -log(tauinv) * _ni
            like_eta0 -= sum(log.(hh))
            while true
                v1 = ones(_ni, 1)
                X2 = zeros(_ni, 2)
                X2[:, 1] = v1 .* h2inv
                X2[:, 2] = D .* h2inv
                xtx2 = X2' * X2
                xty2 = X2' * r2
                beta2 = xtx2 \ xty2
                mu = beta2[1]
                ga = beta2[2]
                eta1 = (ga / tauinv) + (1 - mu / tauinv) * eta0 / 2.0
                eta1 = clamp(eta1, 1e-5, 1e5)
                steps += 1
                while true
                    if abs(eta1 / (1 + eta1) - eta0 / (1 + eta0)) < 1e-6
                        break
                    end
                    hh = eta1 * D + v1
                    h2inv = 1 ./ hh
                    h1inv = sqrt.(h2inv)
                    X = [xs .* h1inv ww2 .* h1inv]
                    xtx = X' * X + diagm(t1)
                    y1 = y .* h1inv
                    xty = X' * y1
                    beta1 = xtx \ xty
                    e1 = y1 - X * beta1
                    r2 = abs2.(e1)
                    tauinv1 = Statistics.mean(r2)
                    like_eta1 = -log(tauinv1) * _ni
                    like_eta1 -= sum(log.(hh))
                    if like_eta1 >= like_eta0
                        break
                    end
                    eta1 = (eta1 + eta0) / 2.0
                end
                if abs(eta1 / (1 + eta1) - eta0 / (1 + eta0)) < 1e-6
                    break
                end
                eta0 = copy(eta1)
                tauinv = copy(tauinv1)
                like_eta0 = copy(like_eta1)
            end
        end
        begin
            tauinv = rh1r / (_ni - _c1) 
            bvar = tauinv * inv(xtx)[1, 1]
            if test_method == "wald"
                tstats = beta1[1] / bvar * beta1[1]
                @runif is_calcu_pv pval = ccdf(d, Float64(tstats)) 
            end
            if test_method == "lrt"
                loglike1 = 0.0
                loglike1 += sum(log.(h1inv))
                loglike1 -= log(sum(r2)) * _ni / 2
                if lrts < 0
                    lrts = 0
                end
                @runif is_calcu_pv pval = ccdf(Chisq(1), Float64(lrts))
            end
            η = eta1 / (1 + eta1)
            β = beta1[1]
            σ = sqrt(abs(bvar))
            if is_calcu_pv
                df_test[xi, :] .= [steps, η, β, σ, pval]
            else
                df_test[xi, :] .= [steps, η, β, σ, tstats]
            end
        end
    end
    return eta1, h, coefs
end
if false
    evals = FloatT.(_struct_KIN.EA.values)
    evals[evals.<0] .= 0
    evecs = _struct_KIN.EA.vectors
    _2p = 2 .* _struct_GENO.snp_annotation.af
    _2pq = _2p .* (1 .- _struct_GENO.snp_annotation.af)
    _2pq = FloatT.(_2pq)
    _2pq = _2pq' |> Matrix
    SQ = _struct_GENO.genotype .- _2p 
    SQ .= evecs' * SQ 
    yQ = evecs' * _struct_PHENO.expression[1, :] |> Vector{FloatT}
    X_MME = _struct_COVAR.X_MME
    xQ = evecs' * X_MME
    xQ = FloatT.(xQ)
end
function idul_assoc_test(SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, is_calcu_pv::Bool=true) where {T}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    _ni, _ns = size(SQ)
    _c = size(xQ, 2)
    _c1 = _c + 1
    xW0 = [ones(FloatT, _ni) xQ]
    d = WaldTest(1, _ni - _c1)
    df_test = zeros(FloatT, _ns, 5)
    @notimeit to @inbounds @threads for xi in 1:_ns
        xW = copy(xW0)
        xW[:, 1] .= SQ[:, xi] 
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            xtxi = zeros(FloatT, _c1, _c1)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
            mat = zeros(FloatT, _c1, _c1)
        end
        eta0 = 1
        eta1 = init_eta
        diff1 = 1
        rh1r = FloatT(0.0)
        for steps in 1:max_iter
            begin
                h2inv .= eta1 .* evals .+ 1 
                h2inv .\= 1 
                h1inv .= sqrt.(h2inv)
                broadcast!(*, X, xW, h1inv)
                mat_mul!(xtx, X', X)
                xtx[diagind(xtx)] .+= FloatT(1e-10)
                broadcast!(*, y1, yQ, h1inv)
                mat_mul!(xty, X', y1)
                begin 
                    try
                        xtxi .= inv(cholesky!(Symmetric(xtx)))
                    catch e
                        @debug string(xi, ": ", e)
                        df_test[xi, :] .= NAN
                        continue
                    end
                    mul!(beta1, xtxi, xty)
                end
                r2 .= y1 .- X * beta1 
                r2 .*= r2 
            end
            diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
            if (diff1 <= thre) | (steps == max_iter)
                begin
                    tauinv = rh1r / (_ni - _c1) 
                    bvar = tauinv * xtxi[1, 1]
                    if test_method == "wald"
                        tstats = beta1[1] / bvar * beta1[1]
                        @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                    end
                    β = beta1[1]
                    σ = sqrt(abs(bvar))
                    if is_calcu_pv
                        df_test[xi, :] .= [steps, eta1, β, σ, pv]
                    else
                        df_test[xi, :] .= [steps, eta1, β, σ, tstats]
                    end
                end
                break
            end
            eta0 = copy(eta1)
            if test_method == "wald"
                rh1r = sum(r2)
                rh2r = dot(r2, h2inv)
                trh1 = sum(h2inv)
                trh2 = dot(h2inv, h2inv)
                mh = trh1 / _ni
                vh = trh2 / _ni - mh * mh
                broadcast!(*, X, xW, h2inv)
                mat_mul!(xtx, X', X) 
                mul!(mat, xtxi, xtx) 
                fp1 = trh1 - tr(mat) - (_ni - _c1) * rh2r / rh1r
                begin
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                end
                eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
            end
            if test_method == "lrt"
            end
        end
    end
    return df_test
end
if false
    evals = FloatT.(_struct_KIN.EA.values)
    evals[evals.<0] .= 0
    evecs = _struct_KIN.EA.vectors
    _2p = 2 .* _struct_GENO.snp_annotation.af
    _2pq = _2p .* (1 .- _struct_GENO.snp_annotation.af)
    _2pq = FloatT.(_2pq)
    _2pq = _2pq' |> Matrix
    SQ = _struct_GENO.genotype .- _2p' 
    SQ .= evecs' * SQ 
    yQ = evecs' * _struct_PHENO.expression[1, :] |> Vector{FloatT}
    X_MME = _struct_COVAR.X_MME
    xQ = evecs' * X_MME
    xQ = FloatT.(xQ)
    df_test = copy(df_subs)
    SQ = copy(cis_genotype[:, idul_prior_index])
    yQ = copy(testpheno)
    xQ
    evals = copy(glo_EA.values)
    max_iter=_args_idul_max_iter
    thre=FloatT.(_args_idul_converge)
    approx=false
    init_eta=FloatT(1.5)
    return_detail=true
    is_calcu_pv=false
    test_method="wald"
    FloatT=Float32
    testdata = [df_test,SQ,yQ,xQ,evals,max_iter,thre,approx,init_eta,return_detail,is_calcu_pv,test_method]
    save("../testdata_idul.jld2","data",testdata)
    df_test,SQ,yQ,xQ,evals,max_iter,thre,approx,init_eta,return_detail,is_calcu_pv,test_method = load("../testdata_idul.jld2")["data"]
    WaldTest = FDist
end
function idul_assoc_test_plus!(df_test::AbstractMatrix{T1}, SQ::AbstractMatrix{T}, yQ::AbstractVector{T}, xQ::AbstractMatrix{T}, evals::AbstractVector{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true,add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c1 = _c + 1
        xW0 = [ones(FloatT, _ni) xQ]
        d = WaldTest(1, _ni - _c1)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c1)
        xtx_init = zeros(FloatT, _c1, _c1)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c1)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        y1_init .= yQ .* h1inv_init
        mat_mul!(xty_init, X_init', y1_init)
        if max_iter > 1
            X2_init = similar(X_init)
            xtx2_init = similar(xtx_init)
            broadcast!(*, X2_init, xW0, h2inv_init)
            mat_mul!(xtx2_init, X2_init', X2_init) 
        end
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c1) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
            xtx2 = zeros(FloatT, _c1, _c1)
            mat = zeros(FloatT, _c1, _c1)
        end
        @views for xi in iter_collects[blocki]
            xW[:, 1] .= SQ[:, xi] 
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                if steps == 1
                    h2inv .= h2inv_init
                    X .= X_init
                    X[:, 1] .= xW[:, 1] .* h1inv_init
                    xtx .= xtx_init
                    xtx[1, :] .= xtx[:, 1] .= vec(X[:, 1]' * X)
                    xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                    xty .= xty_init
                    y1 .= y1_init
                    xty[1] = X[:, 1]' * y1
                else
                    h2inv .= eta1 .* evals .+ 1 
                    h2inv .\= 1 
                    h1inv .= sqrt.(h2inv)
                    broadcast!(*, X, xW, h1inv)
                    mat_mul!(xtx, X', X) 
                    xtx[diagind(xtx)] .+= FloatT(1e-10)
                    broadcast!(*, y1, yQ, h1inv)
                    mat_mul!(xty, X', y1)
                end
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtx, xty) 
                mul!(r2, X, beta1)
                r2 .= y1 .- r2
                r2 .*= r2 
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    tauinv = sum(r2) / (_ni - _c1)
                    bvar = tauinv * xtx[1, 1] 
                    if test_method == "wald"
                        tstats = beta1[1] / bvar * beta1[1]
                        @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                    end
                    β = beta1[1]
                    σ = sqrt(abs(bvar))
                    if return_detail
                        if is_calcu_pv
                            df_test[xi, :] .= [steps, eta1, β, σ, pval]
                        else
                            df_test[xi, :] .= [steps, eta1, β, σ, tstats]
                        end
                    else
                        if is_calcu_pv
                            df_test[xi, :] .= [β, σ, pval]
                        else
                            df_test[xi, :] .= [β, σ, tstats]
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    if steps == 1
                        X .= X2_init
                        X[:, 1] .= xW[:, 1] .* h2inv_init
                        xtx2 .= xtx2_init
                        xtx2[1, :] .= xtx2[:, 1] .= vec(X[:, 1]' * X)
                    else
                        broadcast!(*, X, xW, h2inv)
                        mat_mul!(xtx2, X', X) 
                    end
                    mul!(mat, xtx, xtx2) 
                    fp1 = trh1 - tr(mat) - (_ni - _c1) * rh2r / rh1r
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_assoc_test_approx!(df_test::AbstractMatrix{T1}, SQ::AbstractMatrix{T}, yQ::AbstractVector{T}, xQ::AbstractMatrix{T}, evals::AbstractVector{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), init_eta::Union{T,Nothing}=nothing, is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c1 = _c + 1
        xW0 = [ones(FloatT, _ni) xQ]
        d = WaldTest(1, _ni - _c1)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c1)
        xtx_init = zeros(FloatT, _c1, _c1)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c1)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    broadcast!(*, SQ, SQ, h1inv_init)
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @timeit to "threads" @inbounds @threads for blocki in 1:length(iter_collects)
        begin
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
        end
        X .= X_init
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            begin
                X[:, 1] .= SQ[:, xi]
                xtx .= xtx_init
                xtx[1, :] .= xtx[:, 1] .= vec(X[:, 1]' * X)
                xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                xty[1] = dot(X[:, 1], y1_init)
            end
            begin 
                try
                    if FloatT2 == Float32
                        LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                    else
                        xtx .= inv(Symmetric(xtx)) 
                    end
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtx, xty) 
            end
            mul!(r2, X, beta1)
            r2 .= y1_init .- r2
            r2 .*= r2 
            tauinv = sum(r2) / (_ni - _c1)
            bvar = tauinv * xtx[1, 1] 
            β = beta1[1]
            σ = sqrt(abs(bvar))
            if test_method == "wald"
                tstats = β / bvar * β
                @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                if is_calcu_pv
                    df_test[xi, :] .= [β, σ, pval]
                else
                    df_test[xi, :] .= [β, σ, tstats]
                end
            end
        end
    end
    return df_test
end
function idul_global_test_approx!(df_test::AbstractMatrix{T}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ_h1inv_init_t_y1_init::Vector{T}, X_init_SQ_h1inv_init::Matrix{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), init_eta::Union{T,Nothing}=nothing, is_calcu_pv::Bool=true, pval_threshold::Union{Nothing,T}=nothing) where {T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c1 = _c + 1
        xW0 = [ones(FloatT, _ni) xQ]
        d = WaldTest(1, _ni - _c1)
    end
    if !isnothing(pval_threshold)
        tstat_threshold = quantile(d, 1 - pval_threshold)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c1)
        xtx_init = zeros(FloatT, _c1, _c1)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c1)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        xtx_init[diagind(xtx_init)] .+= FloatT(1e-10)
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    broadcast!(*, SQ, SQ, h1inv_init)
    @timeit to "1.2" begin
        mat_mul!(SQ_h1inv_init_t_y1_init, SQ', y1_init)
        mat_mul!(X_init_SQ_h1inv_init, X_init', SQ)
        X_init_SQ_h1inv_init[1, :] .= vec(sum(abs2, SQ, dims=1)) .+ FloatT(1e-10)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @timeit to "threads" @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        @timeit to "threads-prealloc" begin
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
        end
        @timeit to "threads-init" begin
            X .= X_init
            xty .= xty_init
        end
        @timeit to "for" @views for xi in iter_collects[blocki]
            @timeit to "for 1" begin
                X[:, 1] .= SQ[:, xi]
                xtx .= xtx_init
                xtx[1, :] .= xtx[:, 1] .= X_init_SQ_h1inv_init[:, xi]
                xty[1] = SQ_h1inv_init_t_y1_init[xi]
            end
            @timeit to "for 2" begin 
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtx, xty) 
            end
            @timeit to "for 3" begin
                mul!(r2, X, beta1)
                broadcast!(-, r2, y1_init, r2)
                broadcast!(abs2, r2, r2)
                tauinv = sum(r2) / (_ni - _c1)
                bvar = tauinv * xtx[1, 1] 
                β = beta1[1]
                σ = sqrt(abs(bvar))
            end
            if test_method == "wald"
                tstats = β / bvar * β
                @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                if isnothing(pval_threshold)
                    if is_calcu_pv
                        df_test[xi, :] .= [β, σ, pval]
                    else
                        df_test[xi, :] .= [β, σ, tstats]
                    end
                else
                    if tstats >= tstat_threshold
                        if is_calcu_pv
                            df_test[xi, :] .= [β, σ, pval]
                        else
                            df_test[xi, :] .= [β, σ, tstats]
                        end
                    end
                end
            end
        end
    end
end
function idul_global_test_approx_no_covar!(df_test::AbstractMatrix{T}, SQ::Matrix{T}, yQ::Vector{T}, evals::Vector{T}, SQ_h1inv_init_t_y1_init::Vector{T}, X_init_SQ_h1inv_init::Matrix{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), init_eta::Union{T,Nothing}=nothing, is_calcu_pv::Bool=true, pval_threshold::Union{Nothing,T}=nothing) where {T<:AbstractFloat}
    FloatT = T
    begin
        _ni, _ns = size(SQ)
        _c = 0
        _c1 = _c + 1
        xW0 = ones(FloatT, _ni)
        d = WaldTest(1, _ni - _c1)
    end
    if !isnothing(pval_threshold)
        tstat_threshold = quantile(d, 1 - pval_threshold)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c1)
        y1_init = zeros(FloatT, _ni)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        xtx_init = dot(X_init, X_init) 
        xtx_init += FloatT(1e-10)
        broadcast!(*, y1_init, yQ, h1inv_init)
        xty_init = dot(X_init, y1_init)
    end
    broadcast!(*, SQ, SQ, h1inv_init)
    @timeit to "1.2" begin
        SQ_h1inv_init_t_y1_init = mat_mul(SQ', y1_init)
        X_init_SQ_h1inv_init = vec(sum(abs2, SQ, dims=1)) .+ FloatT(1e-10)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @timeit to "threads" @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        @timeit to "threads-prealloc" begin
            X = zeros(FloatT, _ni, _c1)
            r2 = zeros(FloatT, _ni)
        end
        @timeit to "threads-init" begin
            X .= X_init
        end
        @timeit to "for" @views for xi in iter_collects[blocki]
            @timeit to "for 1" begin
                X[:, 1] .= SQ[:, xi]
                xtx = X_init_SQ_h1inv_init[xi]
                xty = SQ_h1inv_init_t_y1_init[xi]
            end
            @timeit to "for 2" begin 
                    xtxi = 1 / xtx
                beta1 = xtxi * xty
            end
            @timeit to "for 3" begin
                mul!(r2, X, beta1)
                broadcast!(-, r2, y1_init, r2)
                broadcast!(abs2, r2, r2)
                tauinv = sum(r2) / (_ni - _c1)
                bvar = tauinv * xtxi
                β = beta1
                σ = sqrt(abs(bvar))
            end
            if test_method == "wald"
                tstats = β / bvar * β
                @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                if isnothing(pval_threshold)
                    if is_calcu_pv
                        df_test[xi, :] .= [β, σ, pval]
                    else
                        df_test[xi, :] .= [β, σ, tstats]
                    end
                else
                    if tstats >= tstat_threshold
                        if is_calcu_pv
                            df_test[xi, :] .= [β, σ, pval]
                        else
                            df_test[xi, :] .= [β, σ, tstats]
                        end
                    end
                end
            end
        end
    end
end
function idul_assoc_test!(df_test::Matrix{T}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true) where {T}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c1 = _c + 1
        xW0 = [ones(FloatT, _ni) xQ]
        d = WaldTest(1, _ni - _c1)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c1) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c1)
            xtx = zeros(FloatT2, _c1, _c1)
            xtxi = zeros(FloatT, _c1, _c1)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c1)
            beta1 = zeros(FloatT, _c1)
            r2 = zeros(FloatT, _ni)
            mat = zeros(FloatT, _c1, _c1)
        end
        @views for xi in iter_collects[blocki]
            xW[:, 1] .= SQ[:, xi] 
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                begin
                    h2inv .= eta1 .* evals .+ 1 
                    h2inv .\= 1 
                    h1inv .= sqrt.(h2inv)
                    broadcast!(*, X, xW, h1inv)
                    mat_mul!(xtx, X', X) 
                    xtx[diagind(xtx)] .+= FloatT(1e-10)
                    broadcast!(*, y1, yQ, h1inv)
                    mat_mul!(xty, X', y1)
                    begin 
                        try
                            xtxi .= inv(cholesky!(Symmetric(xtx)))
                        catch e
                            @debug string(xi, ": ", e)
                            df_test[xi, :] .= NAN
                            continue
                        end
                        mul!(beta1, xtxi, xty)
                    end
                    r2 .= y1 .- X * beta1 
                    r2 .*= r2 
                end
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    begin
                        tauinv = sum(r2) / (_ni - _c1)
                        bvar = tauinv * xtxi[1, 1]
                        if test_method == "wald"
                            tstats = beta1[1] / bvar * beta1[1]
                            @runif is_calcu_pv pval = ccdf(d, Float64(tstats))
                        end
                        begin
                            β = beta1[1]
                            σ = sqrt(abs(bvar))
                        end
                        if return_detail
                            if is_calcu_pv
                                df_test[xi, :] .= [steps, eta1, β, σ, pval]
                            else
                                df_test[xi, :] .= [steps, eta1, β, σ, tstats]
                            end
                        else
                            if is_calcu_pv
                                df_test[xi, :] .= [β, σ, pval]
                            else
                                df_test[xi, :] .= [β, σ, tstats]
                            end
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    broadcast!(*, X, xW, h2inv)
                    mat_mul!(xtx, X', X) 
                    mul!(mat, xtxi, xtx) 
                    fp1 = trh1 - tr(mat) - (_ni - _c1) * rh2r / rh1r
                    begin
                        delta = fp1 * eta0 / _ni / vh
                        t1 = eta0 + delta
                        t2 = eta0 + delta * mh
                        eta1 = (t1 + t2) / 2
                    end
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_assoc_test!(df_test::Matrix{T}, xs0::Matrix{Int8}, norm_factors::AbstractVecOrMat{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, evecs::Matrix{T}, batch_collects::Vector{UnitRange{Int64}}, S0::Matrix{T}, S1::Matrix{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false) where {T}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    @inbounds @views for batchi in 1:length(batch_collects)
        if batchi < length(batch_collects)
            S0 .= xs0[:, batch_collects[batchi]]
            @timeit to "S0_broadcast!" broadcast!(-, S0, S0, norm_factors[:, batch_collects[batchi]])
            @timeit to "S0.=mat_mul" S0 .= mat_mul(evecs', S0)
            @timeit to "S0 do test" if approx
                idul_assoc_test_approx!(df_test[batch_collects[batchi], :], S0, yQ, xQ, evals; init_eta=init_eta)
            else
                idul_assoc_test_plus!(df_test[batch_collects[batchi], :], S0, yQ, xQ, evals; max_iter=max_iter, init_eta=init_eta)
            end
        else
            S1 .= xs0[:, batch_collects[batchi]]
            @timeit to "S1_broadcast!" broadcast!(-, S1, S1, norm_factors[:, batch_collects[batchi]])
            @timeit to "S1.=mat_mul" S1 .= mat_mul(evecs', S1)
            @timeit to "S1 do test" if approx
                idul_assoc_test_approx!(df_test[batch_collects[batchi], :], S1, yQ, xQ, evals; init_eta=init_eta)
            else
                idul_assoc_test_plus!(df_test[batch_collects[batchi], :], S1, yQ, xQ, evals; max_iter=max_iter, init_eta=init_eta)
            end
        end
    end
end
function idul_two_assoc_test!(df_test::Matrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Matrix{T}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + 2
        xW0 = [ones(FloatT, _ni, 2) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c2) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
            xtx2 = zeros(FloatT, _c2, _c2)
            mat = zeros(FloatT, _c2, _c2)
        end
        @views for xi in iter_collects[blocki]
            xW[:, 1] .= SQ[:, xi] 
            xW[:, 2] .= SQ2[:, xi] 
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                h2inv .= eta1 .* evals .+ 1 
                h2inv .\= 1 
                h1inv .= sqrt.(h2inv)
                broadcast!(*, X, xW, h1inv)
                mat_mul!(xtx, X', X) 
                xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                broadcast!(*, y1, yQ, h1inv)
                mat_mul!(xty, X', y1)
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                        df_test[xi, :] .= NAN
                        continue
                end
                mul!(beta1, xtx, xty) 
                r2 .= y1 .- X * beta1 
                r2 .*= r2 
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    begin
                        tauinv = sum(r2) / (_ni - _c2)
                        bvar1 = tauinv * xtx[1, 1] 
                        bvar2 = tauinv * xtx[2, 2] 
                        if test_method == "wald"
                            tstats1 = beta1[1] / bvar1 * beta1[1]
                            @runif is_calcu_pv pval1 = ccdf(d, Float64(tstats1))
                            tstats2 = beta1[2] / bvar2 * beta1[2]
                            @runif is_calcu_pv pval2 = ccdf(d, Float64(tstats2))
                        end
                        β1 = beta1[1]
                        σ1 = sqrt(bvar1)
                        β2 = beta1[2]
                        σ2 = sqrt(bvar2)
                        if return_detail
                            if is_calcu_pv
                                df_test[xi, :] .= [steps, eta1, β1, σ1, pval1, β2, σ2, pval2]
                            else
                                df_test[xi, :] .= [steps, eta1, β1, σ1, tstats1, β2, σ2, tstats2]
                            end
                        else
                            if is_calcu_pv
                                df_test[xi, :] .= [β1, σ1, pval1, β2, σ2, pval2, NaN]
                            else
                                df_test[xi, :] .= [β1, σ1, tstats1, β2, σ2, tstats2, NaN]
                            end
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    broadcast!(*, X, xW, h2inv)
                    mat_mul!(xtx2, X', X) 
                    mul!(mat, xtx, xtx2) 
                    fp1 = trh1 - tr(mat) - (_ni - _c2) * rh2r / rh1r
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test!(df_test::Matrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    n_tests = length(SQ2) + 1
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c2) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
            xtx2 = zeros(FloatT, _c2, _c2)
            mat = zeros(FloatT, _c2, _c2)
        end
        @views for xi in iter_collects[blocki]
            xW[:, 1] .= SQ[:, xi] 
            [xW[:, ii+1] .= SQ2[ii][:, xi] for ii in 1:n_tests-1]
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                h2inv .= eta1 .* evals .+ 1 
                h2inv .\= 1 
                h1inv .= sqrt.(h2inv)
                broadcast!(*, X, xW, h1inv)
                mat_mul!(xtx, X', X) 
                xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                broadcast!(*, y1, yQ, h1inv)
                mat_mul!(xty, X', y1)
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                        df_test[xi, :] .= NAN
                        continue
                end
                mul!(beta1, xtx, xty) 
                r2 .= y1 .- X * beta1 
                r2 .*= r2 
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    begin
                        tauinv = sum(r2) / (_ni - _c2)
                        bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
                        β_ = beta1[1:n_tests]
                        σ_ = sqrt.(bvar_)
                        if test_method == "wald"
                            tstats_ = β_ .^ 2 ./ bvar_
                            @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
                        end
                        if return_detail
                            if is_calcu_pv
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= pval_
                            else
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= tstats_
                            end
                        else
                            if is_calcu_pv
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= pval_
                            else
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= tstats_
                            end
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    broadcast!(*, X, xW, h2inv)
                    mat_mul!(xtx2, X', X) 
                    mul!(mat, xtx, xtx2) 
                    fp1 = trh1 - tr(mat) - (_ni - _c2) * rh2r / rh1r
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test!(df_test::Matrix{T1}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    n_tests = length(SQ2)
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ2[1])
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c2) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
            xtx2 = zeros(FloatT, _c2, _c2)
            mat = zeros(FloatT, _c2, _c2)
        end
        @views for xi in iter_collects[blocki]
            [xW[:, ii] .= SQ2[ii][:, xi] for ii in 1:n_tests]
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                h2inv .= eta1 .* evals .+ 1 
                h2inv .\= 1 
                h1inv .= sqrt.(h2inv)
                broadcast!(*, X, xW, h1inv)
                mat_mul!(xtx, X', X) 
                xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                broadcast!(*, y1, yQ, h1inv)
                mat_mul!(xty, X', y1)
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                        df_test[xi, :] .= NAN
                        continue
                end
                mul!(beta1, xtx, xty) 
                r2 .= y1 .- X * beta1 
                r2 .*= r2 
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    begin
                        tauinv = sum(r2) / (_ni - _c2)
                        bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
                        β_ = beta1[1:n_tests]
                        σ_ = sqrt.(bvar_)
                        if test_method == "wald"
                            tstats_ = β_ .^ 2 ./ bvar_
                            @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
                        end
                        if return_detail
                            if is_calcu_pv
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= pval_
                            else
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= tstats_
                            end
                        else
                            if is_calcu_pv
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= pval_
                            else
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= tstats_
                            end
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    broadcast!(*, X, xW, h2inv)
                    mat_mul!(xtx2, X', X) 
                    mul!(mat, xtx, xtx2) 
                    fp1 = trh1 - tr(mat) - (_ni - _c2) * rh2r / rh1r
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test!(df_test::Matrix{T1}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}, thread_local_storage::Vector{Dict{String, Array{T}}}; test_method::String="wald", max_iter::Int=100, thre::T=T.(1e-6), approx::Bool=false, init_eta::Union{T,Nothing}=nothing, return_detail::Bool=false, is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    n_tests = length(SQ2)
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    if ((test_method == "lrt") | approx) & (isnothing(init_eta))
        init_eta, _, _ = idul_h(xQ, yQ, evals; max_iter=max_iter, thre=thre, algo="reml")
    end
    if approx
        max_iter = 1
    end
    begin
        _ni, _ns = size(SQ2[1])
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        xW = zeros(FloatT, _ni, _c2) 
        copy!(xW, xW0)
        begin
            h2inv = zeros(FloatT, _ni)
            h1inv = zeros(FloatT, _ni)
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            y1 = zeros(FloatT, _ni)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
            xtx2 = zeros(FloatT, _c2, _c2)
            mat = zeros(FloatT, _c2, _c2)
        end
        @views for xi in iter_collects[blocki]
            [xW[:, ii] .= SQ2[ii][:, xi] for ii in 1:n_tests]
            eta0 = 1
            eta1 = init_eta
            diff1 = 1
            rh1r = FloatT(0.0)
            for steps in 1:max_iter
                h2inv .= eta1 .* evals .+ 1 
                h2inv .\= 1 
                h1inv .= sqrt.(h2inv)
                broadcast!(*, X, xW, h1inv)
                mat_mul!(xtx, X', X) 
                xtx[diagind(xtx)] .+= FloatT(add_to_diag)
                broadcast!(*, y1, yQ, h1inv)
                mat_mul!(xty, X', y1)
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                        df_test[xi, :] .= NAN
                        continue
                end
                mul!(beta1, xtx, xty) 
                r2 .= y1 .- X * beta1 
                r2 .*= r2 
                diff1 = abs(eta1 / (1 + eta1) - eta0 / (1 + eta0))
                if (diff1 <= thre) | (steps == max_iter)
                    begin
                        tauinv = sum(r2) / (_ni - _c2)
                        bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
                        β_ = beta1[1:n_tests]
                        σ_ = sqrt.(bvar_)
                        if test_method == "wald"
                            tstats_ = β_ .^ 2 ./ bvar_
                            @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
                        end
                        if return_detail
                            if is_calcu_pv
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= pval_
                            else
                                df_test[xi, 1:2] .= [steps, eta1]
                                df_test[xi, 3:3:(n_tests*3+2)] .= β_
                                df_test[xi, 4:3:(n_tests*3+2)] .= σ_
                                df_test[xi, 5:3:(n_tests*3+2)] .= tstats_
                            end
                        else
                            if is_calcu_pv
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= pval_
                            else
                                df_test[xi, 1:3:n_tests*3] .= β_
                                df_test[xi, 2:3:n_tests*3] .= σ_
                                df_test[xi, 3:3:n_tests*3] .= tstats_
                            end
                        end
                    end
                    break
                end
                eta0 = copy(eta1)
                if test_method == "wald"
                    rh1r = sum(r2)
                    rh2r = dot(r2, h2inv)
                    trh1 = sum(h2inv)
                    trh2 = dot(h2inv, h2inv)
                    mh = trh1 / _ni
                    vh = trh2 / _ni - mh * mh
                    broadcast!(*, X, xW, h2inv)
                    mat_mul!(xtx2, X', X) 
                    mul!(mat, xtx, xtx2) 
                    fp1 = trh1 - tr(mat) - (_ni - _c2) * rh2r / rh1r
                    delta = fp1 * eta0 / _ni / vh
                    t1 = eta0 + delta
                    t2 = eta0 + delta * mh
                    eta1 = (t1 + t2) / 2
                    eta1 = FloatT(clamp(eta1, 1e-5, 1e5))
                end
                if test_method == "lrt"
                end
            end
        end
    end
    return df_test
end
function idul_two_assoc_test_approx!(df_test::AbstractMatrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Matrix{T}; test_method::String="wald", init_eta::Union{T,Nothing}=nothing, max_iter::Int=100, thre::T=T.(1e-6), is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + 2
        xW0 = [ones(FloatT, _ni, 2) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c2)
        xtx_init = zeros(FloatT, _c2, _c2)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c2)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    broadcast!(*, SQ, SQ, h1inv_init)
    broadcast!(*, SQ2, SQ2, h1inv_init)
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        begin
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
        end
        X .= X_init
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            X[:, 1] .= SQ[:, xi]
            X[:, 2] .= SQ2[:, xi]
            xtx .= xtx_init
            xtx[1:2, :] .= X[:, 1:2]' * X
            xtx[:, 1:2] .= xtx[1:2, :]'
            xtx[diagind(xtx)] .+= FloatT(add_to_diag)
            xty[1:2] .= X[:, 1:2]' * y1_init
            try
                LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
            catch e
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            mul!(beta1, xtx, xty) 
            mul!(r2, X, beta1)
            r2 .= y1_init .- r2
            r2 .*= r2 
            tauinv = sum(r2) / (_ni - _c2)
            bvar1 = tauinv * xtx[1, 1]  
            bvar2 = tauinv * xtx[2, 2] 
            if test_method == "wald"
                tstats1 = beta1[1] / bvar1 * beta1[1]
                @runif is_calcu_pv pval1 = ccdf(d, Float64(tstats1))
                tstats2 = beta1[2] / bvar2 * beta1[2]
                @runif is_calcu_pv pval2 = ccdf(d, Float64(tstats2))
            end
            β1 = beta1[1]
            σ1 = sqrt(bvar1)
            β2 = beta1[2]
            σ2 = sqrt(bvar2)
            if is_calcu_pv
                df_test[xi, :] .= [β1, σ1, pval1, β2, σ2, pval2, NaN]
            else
                df_test[xi, :] .= [β1, σ1, tstats1, β2, σ2, tstats2, NaN]
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test_approx!(df_test::AbstractMatrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}; test_method::String="wald", init_eta::Union{T,Nothing}=nothing, max_iter::Int=100, thre::T=T.(1e-6), is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    n_tests = length(SQ2) + 1
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c2)
        xtx_init = zeros(FloatT, _c2, _c2)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c2)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    broadcast!(*, SQ, SQ, h1inv_init)
    [broadcast!(*, SQ2[ii], SQ2[ii], h1inv_init) for ii in 1:n_tests-1]
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        begin
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
        end
        X .= X_init
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            X[:, 1] .= SQ[:, xi]
            [X[:, ii+1] .= SQ2[ii][:, xi] for ii in 1:n_tests-1]
            xtx .= xtx_init
            xtx[1:n_tests, :] .= X[:, 1:n_tests]' * X
            xtx[:, 1:n_tests] .= xtx[1:n_tests, :]'
            xtx[diagind(xtx)] .+= FloatT(add_to_diag)
            xty[1:n_tests] .= X[:, 1:n_tests]' * y1_init
            try
                LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
            catch e
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            mul!(beta1, xtx, xty) 
            mul!(r2, X, beta1)
            r2 .= y1_init .- r2
            r2 .*= r2 
            tauinv = sum(r2) / (_ni - _c2)
            bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
            β_ = beta1[1:n_tests]
            if any(bvar_ .<= 0)
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            σ_ = sqrt.(bvar_)
            if test_method == "wald"
                tstats_ = β_ .^ 2 ./ bvar_
                @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
            end
            if is_calcu_pv
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= pval_
            else
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= tstats_
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test_approx!(df_test::AbstractMatrix{T1}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}; test_method::String="wald", init_eta::Union{T,Nothing}=nothing, max_iter::Int=100, thre::T=T.(1e-6), is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    n_tests = length(SQ2)
    begin
        _ni, _ns = size(SQ2[1])
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c2)
        xtx_init = zeros(FloatT, _c2, _c2)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c2)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    [broadcast!(*, SQ2[ii], SQ2[ii], h1inv_init) for ii in 1:n_tests]
    thread_local_storage = [Dict("X" => zeros(FloatT, _ni, _c2), "xtx" => zeros(FloatT2, _c2, _c2), "xty" => zeros(FloatT, _c2), "beta1" => zeros(FloatT, _c2), "r2" => zeros(FloatT, _ni)) for _ in 1:nthreads()] 
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects) 
        tid = blocki
        begin
            X = thread_local_storage[tid]["X"]
            xtx = thread_local_storage[tid]["xtx"]
            xty = thread_local_storage[tid]["xty"]
            beta1 = thread_local_storage[tid]["beta1"]
            r2 = thread_local_storage[tid]["r2"]
        end
        X .= X_init
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            [X[:, ii] .= SQ2[ii][:, xi] for ii in 1:n_tests]
            xtx .= xtx_init
            xtx[1:n_tests, :] .= X[:, 1:n_tests]' * X
            xtx[:, 1:n_tests] .= xtx[1:n_tests, :]'
            xtx[diagind(xtx)] .+= FloatT(add_to_diag)
            xty[1:n_tests] .= X[:, 1:n_tests]' * y1_init
            try
                LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
            catch e
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            mul!(beta1, xtx, xty) 
            mul!(r2, X, beta1)
            r2 .= y1_init .- r2
            r2 .*= r2 
            tauinv = sum(r2) / (_ni - _c2)
            bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
            β_ = beta1[1:n_tests]
            if any(bvar_ .<= 0)
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            σ_ = sqrt.(bvar_)
            if test_method == "wald"
                tstats_ = β_ .^ 2 ./ bvar_
                @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
            end
            if is_calcu_pv
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= pval_
            else
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= tstats_
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test_approx!(df_test::AbstractMatrix{T1}, yQ::Vector{T}, xQ::Matrix{T}, evals::Vector{T}, SQ2::Vector{Matrix{T}}, thread_local_storage::Vector{Dict{String, Array{T}}}; test_method::String="wald", init_eta::Union{T,Nothing}=nothing, max_iter::Int=100, thre::T=T.(1e-6), is_calcu_pv::Bool=true, add_to_diag::T=T.(1e-10)) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    if isnothing(init_eta)
        init_eta = FloatT(0.1)
    end
    n_tests = length(SQ2)
    begin
        _ni, _ns = size(SQ2[1])
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        h2inv_init = zeros(FloatT, _ni)
        h1inv_init = zeros(FloatT, _ni)
        X_init = zeros(FloatT, _ni, _c2)
        xtx_init = zeros(FloatT, _c2, _c2)
        y1_init = zeros(FloatT, _ni)
        xty_init = zeros(FloatT, _c2)
        h2inv_init .= init_eta .* evals .+ 1 
        h2inv_init .\= 1 
        h1inv_init .= sqrt.(h2inv_init)
        broadcast!(*, X_init, xW0, h1inv_init)
        mat_mul!(xtx_init, X_init', X_init) 
        broadcast!(*, y1_init, yQ, h1inv_init)
        mat_mul!(xty_init, X_init', y1_init)
    end
    [broadcast!(*, SQ2[ii], SQ2[ii], h1inv_init) for ii in 1:n_tests]
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        tid = blocki
        begin
            X = thread_local_storage[tid]["X"]
            xtx = thread_local_storage[tid]["xtx"]
            xty = thread_local_storage[tid]["xty"]
            beta1 = thread_local_storage[tid]["beta1"]
            r2 = thread_local_storage[tid]["r2"]
        end
        X .= X_init
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            [X[:, ii] .= SQ2[ii][:, xi] for ii in 1:n_tests]
            xtx .= xtx_init
            xtx[1:n_tests, :] .= X[:, 1:n_tests]' * X
            xtx[:, 1:n_tests] .= xtx[1:n_tests, :]'
            xtx[diagind(xtx)] .+= FloatT(add_to_diag)
            xty[1:n_tests] .= X[:, 1:n_tests]' * y1_init
            try
                LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
            catch e
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            mul!(beta1, xtx, xty) 
            mul!(r2, X, beta1)
            r2 .= y1_init .- r2
            r2 .*= r2 
            tauinv = sum(r2) / (_ni - _c2)
            bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
            β_ = beta1[1:n_tests]
            if any(bvar_ .<= 0)
                @debug string(xi, ": ", e)
                df_test[xi, :] .= NAN
                continue
            end
            σ_ = sqrt.(bvar_)
            if test_method == "wald"
                tstats_ = β_ .^ 2 ./ bvar_
                @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
            end
            if is_calcu_pv
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= pval_
            else
                df_test[xi, 1:3:n_tests*3] .= β_
                df_test[xi, 2:3:n_tests*3] .= σ_
                df_test[xi, 3:3:n_tests*3] .= tstats_
            end
        end
    end
    return df_test
end
function idul_two_assoc_test_0eta!(df_test::AbstractMatrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, SQ2::Matrix{T}; test_method::String="wald", is_calcu_pv::Bool=true) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + 2
        xW0 = [ones(FloatT, _ni, 2) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        xtx_init = zeros(FloatT, _c2, _c2)
        xty_init = zeros(FloatT, _c2)
        mat_mul!(xtx_init, xW0', xW0) 
        mat_mul!(xty_init, xW0', yQ)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        begin
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
        end
        X .= xW0
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            begin
                X[:, 1] .= SQ[:, xi]
                X[:, 2] .= SQ2[:, xi]
                xtx .= xtx_init
                xtx[1:2, :] .= X[:, 1:2]' * X
                xtx[:, 1:2] .= xtx[1:2, :]'
                xtx[diagind(xtx)] .+= FloatT(1e-10)
                xty[1:2] .= X[:, 1:2]' * yQ
            end
            begin 
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtx, xty) 
            end
            mul!(r2, X, beta1)
            r2 .= yQ .- r2
            r2 .*= r2 
            begin
                tauinv = sum(r2) / (_ni - _c2)
                bvar1 = tauinv * xtx[1, 1] 
                bvar2 = tauinv * xtx[2, 2] 
                begin
                    β1 = beta1[1]
                    σ1 = sqrt(bvar1)
                    β2 = beta1[2]
                    σ2 = sqrt(bvar2)
                end
                if test_method == "wald"
                    tstats1 = beta1[1] / bvar1 * beta1[1]
                    @runif is_calcu_pv pval1 = ccdf(d, Float64(tstats1))
                    tstats2 = beta1[2] / bvar2 * beta1[2]
                    @runif is_calcu_pv pval2 = ccdf(d, Float64(tstats2))
                end
                if is_calcu_pv
                    df_test[xi, :] .= [β1, σ1, pval1, β2, σ2, pval2, NaN]
                else
                    df_test[xi, :] .= [β1, σ1, tstats1, β2, σ2, tstats2, NaN]
                end
            end
        end
    end
    return df_test
end
function idul_multi_assoc_test_0eta!(df_test::AbstractMatrix{T1}, SQ::Matrix{T}, yQ::Vector{T}, xQ::Matrix{T}, SQ2::Vector{Matrix{T}}; test_method::String="wald", is_calcu_pv::Bool=true) where {T1<:AbstractFloat,T<:AbstractFloat}
    FloatT = T
    n_tests = length(SQ2) + 1
    begin
        _ni, _ns = size(SQ)
        _c = size(xQ, 2)
        _c2 = _c + n_tests
        xW0 = [ones(FloatT, _ni, n_tests) xQ]
        d = WaldTest(1, _ni - _c2)
    end
    begin 
        xtx_init = zeros(FloatT, _c2, _c2)
        xty_init = zeros(FloatT, _c2)
        mat_mul!(xtx_init, xW0', xW0) 
        mat_mul!(xty_init, xW0', yQ)
    end
    block_size = ceil(Int, _ns / nthreads())
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        begin
            X = zeros(FloatT, _ni, _c2)
            xtx = zeros(FloatT2, _c2, _c2)
            xty = zeros(FloatT, _c2)
            beta1 = zeros(FloatT, _c2)
            r2 = zeros(FloatT, _ni)
        end
        X .= xW0
        xty .= xty_init
        @views for xi in iter_collects[blocki]
            begin
                X[:, 1] .= SQ[:, xi]
                [X[:, ii+1] .= SQ2[ii][:, xi] for ii in 1:n_tests-1]
                xtx .= xtx_init
                xtx[1:n_tests, :] .= X[:, 1:n_tests]' * X
                xtx[:, 1:n_tests] .= xtx[1:n_tests, :]'
                xtx[diagind(xtx)] .+= FloatT(1e-10)
                xty[1:n_tests] .= X[:, 1:n_tests]' * yQ
            end
            begin 
                try
                    LinearAlgebra.inv!(cholesky!(Symmetric(xtx))) 
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                mul!(beta1, xtx, xty) 
            end
            mul!(r2, X, beta1)
            r2 .= yQ .- r2
            r2 .*= r2 
            begin
                tauinv = sum(r2) / (_ni - _c2)
                bvar_ = tauinv * diag(xtx[1:n_tests,1:n_tests])
                begin
                    β_ = beta1[1:n_tests]
                    σ_ = sqrt.(bvar_)
                end
                if test_method == "wald"
                    tstats_ = β_.^2 ./ σ_ 
                    @runif is_calcu_pv pval_ = ccdf(d, Float64.(tstats_))
                end
                if is_calcu_pv
                    df_test[xi, 1:3:n_tests*3] .= β_
                    df_test[xi, 2:3:n_tests*3] .= σ_
                    df_test[xi, 3:3:n_tests*3] .= pval_
                else
                    df_test[xi, 1:3:n_tests*3] .= β_
                    df_test[xi, 2:3:n_tests*3] .= σ_
                    df_test[xi, 3:3:n_tests*3] .= tstats_
                end
            end
        end
    end
    return df_test
end
function idul_two_assoc_test_perm!(df_test::AbstractMatrix{T}, SQ::Matrix{T}, y_perms::Matrix{T}, SQ2::Matrix{T}, xW0::Matrix{T}, xtx_init::Matrix{T}, xty_init::Matrix{T}, xtx_SQt_xQ::Matrix{T}, xtx_SQ2t_xQ::Matrix{T}, yQt_SQ::Matrix{T}, yQt_SQ2::Matrix{T}) where {T}
    FloatT = T
    _ni, _ns = size(SQ)
    _c2 = size(xtx_init, 1)
    n_perms = size(xty_init, 2)
    block_size = ceil(Int, _ns / (nthreads() * 10))
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        sub_task = @task begin
            X = similar(xW0)
            xtx = zeros(FloatT2, _c2, _c2)
            xtxi = similar(xtx)
            X .= xW0
            xty_perm = similar(xty_init)
            xty_perm .= xty_init
            beta1_perm = similar(xty_perm)
            r2_perm = zeros(FloatT, _ni, n_perms)
            sum_r2_perm = zeros(FloatT, 1, n_perms)
            tauinv = zeros(FloatT, n_perms)
            bvar = similar(tauinv)
            @views for xi in iter_collects[blocki]
                X[:, 1] .= SQ[:, xi]
                X[:, 2] .= SQ2[:, xi]
                xtx .= xtx_init
                xtx[1, :] .= xtx[:, 1] .= xtx_SQt_xQ[:, xi]
                xtx[2, :] .= xtx[:, 2] .= xtx_SQ2t_xQ[:, xi]
                xtx[diagind(xtx)] .+= FloatT(1e-10)
                try
                    xtxi .= inv(cholesky!(Symmetric(xtx)))
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                xty_perm[1, :] .= yQt_SQ[:, xi]
                xty_perm[2, :] .= yQt_SQ2[:, xi]
                mat_mul!(beta1_perm, xtxi, xty_perm)
                mat_mul!(r2_perm, X, beta1_perm)
                broadcast!(-, r2_perm, y_perms, r2_perm)
                broadcast!(abs2, r2_perm, r2_perm)
                sum!(sum_r2_perm, r2_perm)
                broadcast!(/, tauinv, vec(sum_r2_perm), _ni - _c2)
                broadcast!(*, bvar, tauinv, xtxi[1, 1])
                broadcast!(sqrt, bvar, bvar)
                broadcast!(/, df_test[:, xi*2-1], beta1_perm[1, :], bvar)
                broadcast!(*, bvar, tauinv, xtxi[2, 2])
                broadcast!(sqrt, bvar, bvar)
                broadcast!(/, df_test[:, xi*2], beta1_perm[2, :], bvar)
            end
        end
        schedule(sub_task)
        yield()
        wait(sub_task)
    end
    return df_test
end
function idul_two_assoc_test_perm!(df_test::AbstractMatrix{T}, SQ::Matrix{T}, y_perms::Matrix{T}, SQ2::Matrix{T}, xW0::Matrix{T}, xtx_init::Matrix{T}, xty_init::Matrix{T}, xtx_SQt_xQ::Matrix{T}, xtx_SQ2t_xQ::Matrix{T}) where {T}
    FloatT = T
    _ni, _ns = size(SQ)
    _c2 = size(xtx_init, 1)
    n_perms = size(xty_init, 2)
    block_size = ceil(Int, _ns / (nthreads() * 10))
    iter_collects = collect(Iterators.partition(1:_ns, block_size))
    @timeit to "run tasks-1" @notimeit to @inbounds @threads for blocki in 1:length(iter_collects)
        sub_task = @task begin
            X = similar(xW0)
            xtx = zeros(FloatT2, _c2, _c2)
            xtxi = similar(xtx)
            X .= xW0
            xty_perm = similar(xty_init)
            xty_perm .= xty_init
            beta1_perm = similar(xty_perm)
            r2_perm = zeros(FloatT, _ni, n_perms)
            sum_r2_perm = zeros(FloatT, 1, n_perms)
            tauinv = zeros(FloatT, n_perms)
            bvar = similar(tauinv)
            @views for xi in iter_collects[blocki]
                X[:, 1] .= SQ[:, xi]
                X[:, 2] .= SQ2[:, xi]
                xtx .= xtx_init
                xtx[1, :] .= xtx[:, 1] .= xtx_SQt_xQ[:, xi]
                xtx[2, :] .= xtx[:, 2] .= xtx_SQ2t_xQ[:, xi]
                xtx[diagind(xtx)] .+= FloatT(1e-10)
                try
                    xtxi .= inv(cholesky!(Symmetric(xtx)))
                catch e
                    @debug string(xi, ": ", e)
                    df_test[xi, :] .= NAN
                    continue
                end
                xty_perm[1:2, :] .= X[:, 1:2]' * y_perms
                mat_mul!(beta1_perm, xtxi, xty_perm)
                mat_mul!(r2_perm, X, beta1_perm)
                broadcast!(-, r2_perm, y_perms, r2_perm)
                broadcast!(abs2, r2_perm, r2_perm)
                sum!(sum_r2_perm, r2_perm)
                broadcast!(/, tauinv, vec(sum_r2_perm), _ni - _c2)
                broadcast!(*, bvar, tauinv, xtxi[1, 1])
                broadcast!(sqrt, bvar, bvar)
                broadcast!(/, df_test[:, xi*2-1], beta1_perm[1, :], bvar)
                broadcast!(*, bvar, tauinv, xtxi[2, 2])
                broadcast!(sqrt, bvar, bvar)
                broadcast!(/, df_test[:, xi*2], beta1_perm[2, :], bvar)
            end
        end
        schedule(sub_task)
        yield()
        wait(sub_task)
    end
    return df_test
end
function idul_herit(xQ, yQ, evals)
    algo = "mle"
    _ni, _c = size(xQ)
    v1 = zeros(_ni)
    d1 = zeros(_ni)
    d2 = zeros(_ni)
    X = zeros(_ni, size(xQ, 2))
    e1 = zeros(_ni)
    r2 = zeros(_ni)
    xtx = zeros(size(xQ, 2), size(xQ, 2))
    beta1 = zeros(size(xQ, 2))
    X2 = zeros(_ni, 2)
    xtx2 = zeros(2, 2)
    X3 = zeros(_ni, size(xQ, 2))
    steps = 0
    max_iter = 100
    eta0 = 1.0
    eta1 = 0.1
    diff = 1.0
    for steps in 1:max_iter
        steps += 1
        eta0 = eta1
        dd = eta0 .* evals .+ 1
        d2 = 1 ./ dd
        d1 = sqrt.(d2)
        X = xQ .* d1
        xtx = X' * X
        y1 = yQ .* d1
        xty = X' * y1
        beta1[:] = xtx \ xty
        e1 = y1 - X * beta1
        r2 = e1 .* e1
        if algo == "mle" 
            X2 = [d2 evals .* d2]
            xtx2 = X2' * X2
            xty2 = X2' * r2
            beta2 = xtx2 \ xty2
            t1 = beta2[2] / max(1e-20, beta2[1])
            t2 = beta2[2] / (sum(r2) / _ni)
            eta1 = (t1 + t2) / 2.0
            if eta1 >= 1e5
                eta1 = 1e5
            elseif eta1 <= 1e-5
                eta1 = 1e-5
            end
        elseif algo == "reml" 
            rh1r = sum(r2)
            rh2r = sum(r2 .* d2)
            trh1 = sum(d2)
            trh2 = sum(d2 .* d2)
            mh = trh1 / _ni
            vh = trh2 / _ni - mh * mh
            X3[:, :] = xQ .* d2
            xtx2[:, :] = X3' * X3
            mat = xtx \ xtx2
            fp1 = trh1 - trace(mat) - (_ni - _c) * rh2r / rh1r
            eta1 = eta0 + fp1 * eta0 / _ni / vh
            if eta1 >= 1e5
                eta1 = 1e5
            elseif eta1 <= 1e-5
                eta1 = 1e-5
            end
        end
        diff = eta1 / (1 + eta1) - eta0 / (1 + eta0)
        if (abs(diff) <= 1e-6) | (steps == max_iter)
            @debug string(eta1)
            @debug string(eta1 / (1 + eta1))
            break
        end
    end
end
