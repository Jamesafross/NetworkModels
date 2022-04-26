function fitR(modelFC,realFC)
    N = size(modelFC,1)
    modelFCLT = zeros(Int((N^2-N)/2))
    realFCLT = zeros(Int((N^2-N)/2))
    c = 1
    for i = 2:N
            for j = 1:i-1
                    modelFCLT[c] = modelFC[i,j]
                    realFCLT[c] = realFC[i,j]
                    c+=1
            end
    end
    return cor(modelFCLT,realFCLT)
end

function make_uhist(tgrid,u)
    sizeT = size(tgrid,1)
    t = tgrid[1]:(tgrid[end]-tgrid[1])/(sizeT -1):tgrid[end]
    interp = []
    for i = 1:size(u,1)
        if i == 1
            interp = [CubicSplineInterpolation(t,u[i,:])]
        else
            interp = cat(interp,[CubicSplineInterpolation(t,u[i,:])],dims=1) 
        end
    end
    return interp
end

function adapt_global_coupling(p,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64})
    @inbounds for ii = 1:N
        @inbounds for jj = 1:N
            if W[jj,ii] != 0.0
                if lags[jj,ii] == 0.0
                     W[jj,ii] += 0.1*u[ii]*(u[jj] - h(p,t-1.0;idxs=jj))
                else
                     W[jj,ii] += 0.1*h(p,t-lags[jj,ii];idxs=ii)*(u[jj] - h(p,t-1.0;idxs=jj))
                end
            end
        end
        W[W .< 0.0] .= 0.0
        if sum(W[:,ii]) != 0.0
        @views W[:,ii] = W[:,ii]./sum(W[:,ii])
        end

    end

     @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end

    return W
  
end

function h1(hparams,t;idxs = nothing)
    #history function used on first window
    u0 = hparams
        if t < 0
        return u0[idxs]
    end
end

function h2(hparams,t;idxs = nothing)
    #history function used on windows > 1
    u_hist= hparams
        if t < 0
            return u_hist[idxs](t)
        end
end

function normalise(W,N)

   for ii = 1:N
    W[W .< 0.0] .= 0.0
        if sum(W[:,ii]) != 0.0
        @views W[:,ii] = W[:,ii]./sum(W[:,ii])
        end

    end

     @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end
    return W
end


function makeInitConds(NGp)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    
    params = κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    rE0, rI0, vE0, vI0, gEE0, gEI0, gIE0, gII0 = init_conds_SS(params)
    perturb = 0.1*rand(8*N)

    u0 = zeros(8*N)
    u0[1:N] .= rE0
    u0[N+1:2N] .= rI0
    u0[2N+1:3N] .= vE0
    u0[3N+1:4N] .= vI0
    u0[4N+1:5N] .= gEE0
    u0[5N+1:6N] .= gIE0
    u0[6N+1:7N] .= gEI0
    u0[7N+1:8N] .= gII0

    u0 = u0 + perturb
end
