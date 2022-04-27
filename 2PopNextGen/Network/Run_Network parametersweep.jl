#includes
using Distributed,SharedArrays
if nprocs() < 2
    addprocs(1)
    println("Number of Workers = ", nworkers())
end 

@everywhere begin 
    using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters
    HOMEDIR = homedir()
    WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    include("./functions/stability.jl")
    include("./functions/DEfunctions.jl")
    include("./functions/functions.jl")
    include("./functions/parameters.jl")
    include("./modelFunc.jl")
    include("../../Balloon_Model/balloonModelFunctions.jl")
    include("../../Balloon_Model/balloonModelRHS.jl")
    include("../../Balloon_Model/parameter_sets.jl")
    normaliseSC = 0
    #load data and make struct & dist matrices
    SC = load("$InDATADIR/PaulSC.jld","C")
    dist = load("$InDATADIR/PaulDist.jld","dist")
    PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")
    N = size(SC,1) # number of nodes
    c =7000. # conductance velocity
    lags = round.(dist./c,digits=2) # axonal delays
    lags[lags .== 0.001] .= 0.0
    clags = reshape(lags[lags.>0.0],length(lags[lags.>0.0]))
    stimNodes = [21,39]
    Tstim = [60,90]
    if normaliseSC == 1
        SC = normalise(SC,N)
        W_sum = ones(N)
    else
        W_sum = zeros(N)
        for i = 1:N
            W_sum[i] = sum(SC[:,i])
        end 
    end
    minSC = minimum(SC[SC.>0.0])
    W = zeros(N,N)
    W.=SC
    NGp = NextGen2PopParams()
    W = W+(1/NGp.κ)diagm(ones(N))
    bP = ballonModelParameters()
end

nWindows = 1
tWindows = 300

stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)
nTrials1 = 20
nTrials2 = 20
κVec = LinRange(0.09,0.11,nTrials1)
η_0EVec = LinRange(-14.0,-14.4,nTrials2)

Rsave = SharedArray(zeros(N,N,nWindows,nTrials1,nTrials2))
for i = 1:nTrials1
    @sync @distributed for j = 1:nTrials2
        println("Trial: ",(i-1)*nTrials2 + j," out of ", nTrials1*nTrials2 )
        NGp = NextGen2PopParams(κ=κVec[i],η_0E=η_0EVec[j])
        Rsave[:,:,:,i,j],Wsave = NGModelRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
    end
end



fit = zeros(nWindows)

if nWindows > 1
    meanR = mean(Rsave[:,:,:,:,:],dims=3)[:,:,1,:,:]
    meanfit = fitR(meanR,PaulFCmean)
else
    meanfit = fit
    meanR = Rsave[:,:,1,:,:]
end

pSweep = Array{pSweepData}(undef,nTrials1,nTrials2)
for i = 1:nTrials1
    for j = 1:nTrials2
        
        fitp = fitR(meanR[:,:,i,j],PaulFCmean)
        pSweep[i,j] = pSweepData(κVec[i],η_0EVec[j],fitp)
    end
end





save("./data/pSweep.jld","pSweep",pSweep)