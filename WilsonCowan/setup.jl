if parallel = "on"
    using Distributed,LinearAlgebra,SharedArrays,Plots

    if nprocs() < 4
        addprocs(4)
        println("Number of Workers = ", nworkers())
    end
    #includes

    @everywhere begin 
        normaliseSC = 0
        using StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD,LinearAlgebra,Interpolations
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
        c=7000.
        SC,minSC,W_sum,lags,PaulFCmean,N = getData(c;normalise=0,delayDigits=2)
        SC_Array,minSC_vec,W_sum_mat,lags,FC_Array,N = getData_nonaveraged(c;normalise=0,delayDigits=2)
       
        SC = SC[1:size(SC,1) .!= 100,1:size(SC,1) .!= 100 ]
        SC=SC_Array[:,:,2]
        W_sum = W_sum[1:size(W_sum,1) .!= 100]
        W_sum = W_sum_mat[:,2]
        N = size(SC,1)
        PaulFCmean = PaulFCmean[1:size(PaulFCmean,1) .!= 100,1:size(PaulFCmean,1) .!= 100 ]
        PaulFCmean = PaulFCmean .- diagm(ones(N))

    end
else 
    using LinearAlgebra,Plots,StochasticDelayDiffEq,Parameters,Statistics,StatsBase,DifferentialEquations,JLD,LinearAlgebra,Interpolations


