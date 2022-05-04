#includes
using Distributed,SharedArrays
if nprocs() < 4
    addprocs(4)
    println("Number of Workers = ", nworkers())
end 

@everywhere begin 
    using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters
    HOMEDIR = homedir()
    WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    include("./functions/NextGenFunctions.jl")
    include("../../Balloon_Model/BalloonModel.jl")
    include("$InDATADIR/getData.jl")
    
    #load data and make struct & dist matrices
    SC,minSC,W_sum,dist,PaulFCmean,N = getData2(;normalise=0,delayDigits=2)
    SC = SC[1:size(SC,1) .!= 100,1:size(SC,1) .!= 100 ]
    W_sum = W_sum[1:size(W_sum,1) .!= 100]
    N = size(SC,1)
    PaulFCmean = PaulFCmean[1:size(PaulFCmean,1) .!= 100,1:size(PaulFCmean,1) .!= 100 ]
    PaulFCmean = PaulFCmean .- diagm(ones(N))
    W = zeros(N,N)
    W.=SC

    bP = ballonModelParameters()
    stimNodes = [21,39]
    Tstim = [60,90]
 
end

nWindows = 1
tWindows = 300

stimOpts = "off"
adapt = "off"
opts=modelOpts(stimOpts,adapt)
@everywhere nTrials0 = 1
@everywhere nTrials1 = 10
@everywhere nTrials2 = 30
cvec = [6000]

κVec = LinRange(0.09,0.12,nTrials1)
η_0EVec = LinRange(-14.5,-13.9,nTrials2)
@everywhere nTrials = nTrials0*nTrials1*nTrials2
@everywhere countn = reshape(LinRange(1,nTrials,nTrials),nTrials0,nTrials1,nTrials2)
@everywhere countn = Int.(round.(countn))

Rsave = SharedArray(zeros(N,N,nTrials0,nTrials1,nTrials2))
for k = 1:nTrials0
    global c=cvec[k]
        
  
    for i = 1:nTrials1
        @sync @distributed for j = 1:nTrials2
            lags = dist./c
           
            println("Trial: ",countn[k,i,j]," out of ", nTrials1*nTrials2*nTrials0 )
            NGp = NextGen2PopParams(κ=κVec[i],η_0E=η_0EVec[j])
            println(NGp.κ,NGp.η_0E)
            Rsave[:,:,k,i,j],Wsave = NGModelRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
            
        end
    end
end


fit = zeros(nWindows)

if nWindows > 1
    meanR = mean(Rsave[:,:,:,:,:],dims=3)[:,:,:,:,:]
    meanfit = fitR(meanR,PaulFCmean)
else
    meanfit = fit
    meanR = Rsave[:,:,:,:,:]
end

pSweep = Array{pSweepData}(undef,nTrials0,nTrials1,nTrials2)
for k = 1:nTrials0
    for i = 1:nTrials1
        for j = 1:nTrials2
            
            fitp = fitR(meanR[:,:,k,i,j],PaulFCmean)
            pSweep[k,i,j] = pSweepData(cvec[k],κVec[i],η_0EVec[j],fitp)
        end
    end
end


fitMat = zeros(nTrials0,nTrials1,nTrials2)
for k = 1:nTrials0
    for i = 1:nTrials1
        for j = 1:nTrials2
            fitMat[k,i,j] = pSweep[k,i,j].fit
        end
    end
end







save("./data/pSweep.jld","pSweep",pSweep)

