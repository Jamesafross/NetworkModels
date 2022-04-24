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

function adapt_global_coupling(hparams,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64},minSC::Float64,W_sum::Vector{Float64})
    @inbounds for ii = 1:N
        @inbounds for jj = 1:N
            if W[jj,ii]  > 0.0
                if lags[jj,ii] == 0.0
                     W[jj,ii] += 0.1*u[ii]*(u[jj] - h(hparams,t-1.0;idxs=jj))
                else
                     W[jj,ii] += 0.1*h(hparams,t-lags[jj,ii];idxs=ii)*(u[jj] - h(hparams,t-1.0;idxs=jj))
                end
            end
        end
        W[0. .< W .< minSC] .= minSC
        if sum(W[:,ii]) != 0.0
        @views W[:,ii] = W_sum[ii].*(W[:,ii]./sum(W[:,ii]))
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


f(x::Float64,β::Float64,θ::Float64) = 1/(1+exp(-β*(x-θ)))

function stim(t,i,stimNodes,Tstim,nRun,stimOpt)
    if i ∈ stimNodes && (Tstim[1] <t < Tstim[2]) && (stimOpt == "on" || stimOpt == "ON")
        return -1.
    else
        return 0.
    end
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




