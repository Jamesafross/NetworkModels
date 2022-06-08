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

function adapt_global_coupling(hparams,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64},minSC::Float64,W_sum::Vector{Float64})
    @inbounds for ii = 1:N

        @inbounds for jj = 1:N
         
            if W[ii,jj]  > 0.0
                if lags[ii,jj] == 0.0
                    
                     W[ii,jj] += 0.0000001*u[ii]*(u[jj] - h(hparams,t-1.0;idxs=jj))
                else
                     W[ii,jj] += 0.0000001*h(hparams,t-lags[ii,jj];idxs=ii)*(u[jj] - h(hparams,t-1.0;idxs=jj))
                end
                if W[jj,ii] < minSC
                    W[jj,ii] = minSC
                elseif W[jj,ii] > 0.11
                    W[jj,ii] = 0.11
                end
            end
        
        end
       
        
        if sum(W[:,ii]) != 0.0
        @views W[:,ii] = W_sum[ii].*(W[:,ii]./sum(W[:,ii]))
        end
        W[W .< 0.0] .= 0.0

    end

     @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end

    return W

end

function rij(cij,bsdp)
    return (cij + 1.0)^bsdp
end

function hill(x::Float64,hsdp::Float64,bsdp::Float64)
    return x/(x+hsdp^bsdp) - 1/2
end


function modHeaviside(x)
    if x < 0.0
        return 1.0
    else
        return -1.0
    end
end

function ΔWsdp(cij::Float64,hsdp::Float64,bsdp::Float64)
    return 0.0001*hill(rij(cij,bsdp),hsdp,bsdp)
end

function ΔWgtp(w,dist,cij,agdp,ηgdp)
    return agdp*modHeaviside(w - exp(cij*dist))*ηgdp
end

function getHistMat(HISTMAT,h,u,hparams,lags,t,N)
    for i = 1:N
        for j = 1:N
            if lags[i,j] > 0
                HISTMAT[i,j] = h(hparams,t-lags[i,j];idxs=j)#
            else
                HISTMAT[i,j] = u[j]
            end
        end
    end
end


function adapt_local_func(h,hparams,t,κS,NGp,rE,rI,i,N,c;type = "lim")
        @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
        κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
        @unpack κSEEv,κSIEv,κSEIv,κSIIv,κSUM = κS
        
        κSEEv[i] = (κSEEv[i] + c*rE*(rE - h(hparams,t-1.0;idxs = i)))
        κSIEv[i] = (κSIEv[i] + c*rE*(rI - h(hparams,t-1.0;idxs = i+N)))
        κSEIv[i] = (κSEIv[i] + c*rI*(rE - h(hparams,t-1.0;idxs = i)))
        κSIIv[i] = (κSIIv[i] + c*rI*(rI - h(hparams,t-1.0;idxs = i+N)))
        
        limEE = 0.2
        limEI = 0.2
        limIE = 0.2
        limII = 0.2

        if type == "lim"
            if κSEEv[i]  > κSEE + limEE
                κSEEv[i] = κSEE + limEE
            elseif κSEEv[i]  < κSEE - limEE
                κSEEv[i] = κSEE - limEE
            end

            if κSEIv[i]  > κSEI + limEI
                κSEIv[i] = κSEI + limEI
            elseif κSEIv[i]  < κSEI - limEI
                κSEIv[i] = κSEI - limEI
            end

            if κSIEv[i]  > κSIE + limIE
                κSIEv[i] = κSIE + limIE
            elseif κSIEv[i]  < κSIE - limIE
                κSIEv[i] = κSIE - limIE
            end

            if κSIIv[i]  > κSII + limII
                κSIIv[i] = κSII + limII
            elseif κSIIv[i]  < κSII - limII
                κSIIv[i] = κSII - limII
            end
        elseif type == "normalised"
            κSEEv[i], κSIEv[i],κSEIv[i], κSIIv[i] = κSUM*[κSEEv[i], κSIEv[i], κSEIv[i], κSIIv[i]]/(κSEEv[i] + κSIEv[i] + κSEIv[i] + κSIIv[i])
        end

    return κSEEv[i],κSIEv[i],κSEIv[i],κSIIv[i]
end



function adapt_global_coupling_cor(N::Int64,W::Matrix{Float64},dist::Matrix{Float64},minSC::Float64,W_sum::Vector{Float64},HIST::Array{Float64},hsdp,bsdp)
  
    @inbounds for ii = 1:N
        
        @inbounds for jj = 1:N 
            if W[jj,ii]  > 0.0
                cij = cor(HIST[ii,:],HIST[jj,:])
                W[jj,ii] += ΔWsdp(cij,hsdp,bsdp)
                if W[jj,ii] < minSC
                    W[jj,ii] = minSC
                elseif W[jj,ii] > 0.12
                    W[jj,ii] = 0.12
                end
            end
        
        end
       
        if sum(W[:,ii]) != 0.0
        @views W[:,ii] = W_sum[ii].*(W[:,ii]./sum(W[:,ii]))
        end
        W[W .< 0.0] .= 0.0

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


function makeInitConds(NGp,N)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    
    params = κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    rE0, rI0, vE0, vI0, gEE0, gEI0, gIE0, gII0 = init_conds_SS(params)
    perturb = 0.0*rand(8*N)

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

function stim(t,i,stimNodes,Tstim,nWindow,stimOpt,stimWindow,stimStr)

    if i ∈ stimNodes && (Tstim[1] <t < Tstim[2]) && (stimOpt == "on" || stimOpt == "ON") && nWindow == stimWindow
        return stimStr
    else 
        return 0.
    end
end

function findBestFit(R,FC_Array)
    fitRvec = zeros(size(FC_Array,3))
    for i = 1:size(FC_Array,3);
        fitRvec[i] = fitR(R,FC_Array[:,:,i]);
    end

    bestElementsTest = findfirst(x->x==sort(fitRvec, rev=true)[1],fitRvec)
    bestElements = findfirst(x->x==sort(fitRvec, rev=true)[1],fitRvec)

    sort(fitRvec, rev=true)
    for i = 1:size(FC_Array,3)
        currentFit = fitR(R,mean(FC_Array[:,:,bestElements],dims=3)[:,:,1])
        best = findfirst(x->x==sort(fitRvec, rev=true)[i],fitRvec)

        if fitR(R,mean(FC_Array[:,:,cat(bestElements,best,dims=1)],dims=3)[:,:,1]) > currentFit
            bestElements = cat(bestElements,findfirst(x->x==sort(fitRvec, rev=true)[i],fitRvec),dims=1)
        end
    end

    fit = fitR(R,mean(FC_Array[:,:,bestElements],dims=3)[:,:,1])
    return fit,bestElements
    
end
        
        

