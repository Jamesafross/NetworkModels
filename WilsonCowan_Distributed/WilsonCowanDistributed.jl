using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 2
    addprocs(2)
    println("Number of Workers = ", nworkers())
end
#includes


@everywhere begin 
    normaliseSC = 0
    using StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD,Interpolations
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
    include("$InDATADIR/getData.jl")
    #load data and make struct & dist matrices
    SC,minSC,W_sum,lags,PaulFCmean,N = getData(c;normalise=0,delayDigits=2)



  
end


# get parameters and make structures


WCp= WCparams()
bP = ballonModelParameters()

# Stimulation Setup


nWindows = 8
tWindows = 300.0
nTrials = 1
fitArray = zeros(nWindows,nTrials)
R_Array = SharedArray(zeros(N,N,nWindows,nTrials))
W_save = SharedArray(zeros(N,N,nWindows,nTrials))
W = zeros(N,N)
W .= SC
stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)


@sync @distributed for i = 1:nTrials
    println("working on Trial: ",i)

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

    
