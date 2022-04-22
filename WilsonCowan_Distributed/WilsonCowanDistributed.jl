using Distributed,LinearAlgebra,SharedArrays,Plots

if nprocs() < 5
    addprocs(5)
    println("Number of Workers = ", nworkers())
end
#includes

@everywhere begin 
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

    C = load("$WORKDIR/data/PaulSC.jld","C")
    dist = load("$WORKDIR/data/PaulDist.jld","dist")
    N = size(C,1) # number of nodes
    c =7000. # conductance velocity
    lags = round.(dist./c,digits=2) # axonal delays
    stimNodes = [21,39]
    Tstim = [60,90]
end
PaulFCmean = load("$WORKDIR/data/PaulFCmean_140.jld","paulFCmean_140")


# get parameters and make structures

WCp= WCparams()
bP = ballonModelParameters()

# Stimulation Setup


nWindows = 5
tWindows = 500.0
nTrials = 10

R_Array = SharedArray(zeros(N,N,nWindows,nTrials))


@sync @distributed for i = 1:nTrials

    R_Array[:,:,:,i] = WCModelRun(WCp,bP,nWindows,tWindows,C,lags,N)
end

save("$WORKDIR/data/R_array.jld","R_Array","R_Array")    
