using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 8
    addprocs(8)
    println("Number of Workers = ", nworkers())
end
#includes


@everywhere begin 
    normaliseSC = 0
    using StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD
    include("functions/functions.jl")
    include("../Balloon_Model/balloonModelFunctions.jl")
    include("../Balloon_Model/balloonModelRHS.jl")
    include("../Balloon_Model/parameter_sets.jl")
    include("functions/parameters.jl")
    include("functions/DEfunctions.jl")
    include("functions/modelFunc.jl")
    HOMEDIR=homedir()
    WORKDIR="$HOMEDIR/NetworkModels/WilsonCowan_Distributed"
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    #load data and make struct & dist matrices

    SC = load("$InDATADIR/PaulSC.jld","C")
    dist = load("$InDATADIR/PaulDist.jld","dist")
    N = size(SC,1) # number of nodes
    c =7000. # conductance velocity
    lags = round.(dist./c,digits=2) # axonal delays
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
end
PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")

# get parameters and make structures


WCp= WCparams()
bP = ballonModelParameters()

# Stimulation Setup


nWindows = 8
tWindows = 300.0
nTrials = 40

R_Array = SharedArray(zeros(N,N,nWindows,nTrials))
fitArray = zeros(nWindows,nTrials)
W_save = SharedArray(zeros(N,N,nWindows,nTrials))
W = zeros(N,N)
W .= SC
stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)
etavec = LinRange(0.08,0.15,nTrials)

@sync @distributed for i = 1:nTrials
    WCp = WCparams(η= etavec[i])
    R_Array[:,:,:,i],W_save[:,:,:,i] = WCRun(WCp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
end

for i = 1:nWindows
	for j = 1:nTrials
		fitArray[i,j] = fitR(R_Array[:,:,i,j],PaulFCmean)
	end
end

global dataWC = []
for i = 1:nTrials
	global dataWC = cat(dataWC,dataStruct(R_Array[:,:,:,i],fitArray[:,i],W_save[:,:,:,i]),dims=1)
end 
save("$WORKDIR/data/dataWC.jld","dataWC",dataWC)

    
