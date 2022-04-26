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

function getHCPdata()
    HOMEDIR=homedir()
    DIR = "$HOMEDIR/PostDocFiles/connectivity_data/HCP"

    file = matopen("$DIR/HCPstr.mat");
    W= read(file, "W")
    close(file);
    
    file = matopen("$DIR/HCPdist.mat");
    dist = read(file, "dist")
    close(file);
    N = size(dist,1)

    return W,dist,N
end