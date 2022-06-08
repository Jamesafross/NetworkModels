using LinearAlgebra,Plots,StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD,LinearAlgebra,Interpolations

    include("functions/functions.jl")
    include("../Balloon_Model/balloonModelFunctions.jl")
    include("../Balloon_Model/balloonModelRHS.jl")
    include("../Balloon_Model/parameter_sets.jl")
    include("functions/Structures.jl")
    include("functions/DEfunctions.jl")
    include("functions/modelFunc.jl")
    HOMEDIR=homedir()
    WORKDIR="$HOMEDIR/NetworkModels/WilsonCowan_Distributed"
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    include("$InDATADIR/getData.jl")
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

    stimOpt = "off"
    stimWindow = 2
    stimNodes = [37]
    Tstim = [30,60]
    adapt = "off"
    tWindows = 300.
    nWindows = 1
    ISP = "off"

    opts=solverOpts(stimOpt,stimWindow,stimNodes,Tstim,adapt,tWindows,nWindows,ISP)
    nP = networkParameters(W, dist,lags, N,minSC,W_sum)