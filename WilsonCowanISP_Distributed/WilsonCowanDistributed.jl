using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 2
    addprocs(1)
    println("Number of Workers = ", nworkers())
end
#includes
@everywhere begin 
    normaliseSC = 0
    using StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD
    HOMEDIR=homedir()
    WORKDIR="$HOMEDIR/NetworkModels/WilsonCowanISP_Distributed"
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    include("functions/functions.jl")
    include("../Balloon_Model/balloonModelFunctions.jl")
    include("../Balloon_Model/balloonModelRHS.jl")
    include("../Balloon_Model/parameter_sets.jl")
    include("functions/parameters.jl")
    include("functions/DEfunctions.jl")
    include("functions/modelFunc.jl")
    include("$InDATADIR/getData.jl")


    stimNodes = [21,39]
    Tstim = [60,90]
    c=7000
    SC,minSC,W_sum,lags,PaulFCmean,N = getData(c;delayDigits=3)
end


# get parameters and make structures
WCp= WCparams()
bP = ballonModelParameters()

# Stimulation Setup


nWindows = 4
tWindows = 300.0
nTrials = 1

R_Array = SharedArray(zeros(N,N,nWindows,nTrials))
fitArray = zeros(nWindows,nTrials)
W_save = SharedArray(zeros(N,N,nWindows,nTrials))
W = zeros(N,N)
W .= SC
stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)
#etavec = LinRange(0.08,0.15,nTrials)
WCp= WCparams(η=0.14)




R_Array,W_save= DistributedLoop(R_Array,W_save,WCp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)

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

heatmap(R_Array[:,:,1,1])

    
