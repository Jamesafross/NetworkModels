using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 2
    addprocs(2)
    println("Number of Workers = ", nworkers())
end
#includes


@everywhere begin 
    normaliseSC = 1
    using StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD
    include("functions.jl")
    include("../Balloon_Model/balloonModelFunctions.jl")
    include("../Balloon_Model/balloonModelRHS.jl")
    include("../Balloon_Model/parameter_sets.jl")
    include("parameters.jl")
    include("DEfunctions.jl")
    include("modelFunc.jl")
    HOMEDIR=homedir()
    WORKDIR="$HOMEDIR/NetworkModels/WilsonCowan_Distributed"
    #load data and make struct & dist matrices

    SC = load("$WORKDIR/data/PaulSC.jld","C")
    dist = load("$WORKDIR/data/PaulDist.jld","dist")
    N = size(SC,1) # number of nodes
    c =7000. # conductance velocity
    lags = round.(dist./c,digits=2) # axonal delays
    stimNodes = [21,39]
    Tstim = [60,90]
    if normaliseSC == 1
        SC = normalise(SC,N)
    end
    minSC = minimum(SC[SC.>0.0])
end
PaulFCmean = load("$WORKDIR/data/PaulFCmean_140.jld","paulFCmean_140")

# get parameters and make structures


WCp= WCparams()
bP = ballonModelParameters()

# Stimulation Setup


nWindows = 1
tWindows = 300.0
nTrials = 1

R_Array = SharedArray(zeros(N,N,nWindows,nTrials))
W = zeros(N,N)
W .= SC
stimOpts = "on"
adapt = "off"


opts= modelOpts(stimOpts,adapt)

@sync @distributed for i = 1:nTrials

    R_Array[:,:,:,i],W[:,:] = WCModelRun(WCp,bP,nWindows,tWindows,W,lags,N,minSC,opts)
end

save("$WORKDIR/data/R_array.jld","R_Array","R_Array")    
