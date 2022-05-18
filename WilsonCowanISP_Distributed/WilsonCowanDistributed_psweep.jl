using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 4
    addprocs(4)
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
    include("$InDATADIR/getData.jl")
    stimNodes = [21,39]
    Tstim = [60,90]
    #load data and make struct & dist matrices
    c=7000.
    SC_Array,FC_Array,dist = getData_nonaveraged(;SCtype="log")
    FC_Array = FC_Array
    PaulFCmean = mean(FC_Array,dims=3)
    SC = 0.01*SC_Array[:,:,1]
    lags = dist./c
    lags = round.(lags,digits=2) 
    lags[lags.<0.003] .= 0.000
    #lags[SC .< 0.018] .= 0  
    minSC,W_sum=getMinSC_and_Wsum(SC)
    N = size(SC,1)
    W = zeros(N,N)
    W.=SC
end


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

stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)
etavec = LinRange(0.08,0.15,nTrials)

@sync @distributed for i = 1:nTrials
    WCp = WCparams(η= etavec[i])
    R_Array[:,:,:,i],W_save[:,:,:,i] = WCModelRun(WCp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
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

    
