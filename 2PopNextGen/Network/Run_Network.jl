#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations,MKL

HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels_Data/StructDistMatrices"
OutDATADIR="$HOMEDIR/NetworkModels_Data/2PopNextGen_Data"
include("./functions/NextGenFunctions.jl")
include("../../Balloon_Model/BalloonModel.jl")
include("$InDATADIR/getData.jl")

Run_vec = [1]
plot_fit = "false"
save_data = "true"

# constants 
N = 139
u0 = zeros(8*N)
NGp = NextGen2PopParams2(η_0E = -14.19,κ=0.505)
LR = 0.00001 # learning rate adaptation
u0[:] = makeInitConds(NGp,N)

for jj = 1:length(Run_vec)

    for setstim = ["on","off"]

    stimNodes = [39]
    Tstim = [60,90]

    #load data and make struct & dist matrices

    c=7000.
    SC_Array,FC_Array,dist = getData_nonaveraged(;SCtype="log")
    FC_Array = FC_Array
    PaulFCmean = mean(FC_Array,dims=3)
    lags = dist./c
    lags = round.(lags,digits=2) 
    lags[lags.<0.003] .= 0.000
    #lags[SC .< 0.018] .= 0  
    SC = 0.01SC_Array[:,:,1]
    minSC,W_sum=getMinSC_and_Wsum(SC)
    N = size(SC,1)
    W = zeros(N,N)
    W.=SC
    NGp = NextGen2PopParams2(η_0E = -14.19,κ=0.505)
    #NGp = NextGen2PopParams3(η_0E = -5.,κ=0.5)
    nP = networkParameters(W, dist,lags, N,minSC,W_sum)

    Run = string(Int(round(Run_vec[jj])))
    nWindows = 22
    tWindows = 100.0
    stimOpt = setstim
    stimStr = 5.
    stimWindow = 5
    adapt = "on"
    synapses = "1stOrder"



    κSEEv = ones(N)*NGp.κSEE
    κSIEv = ones(N)*NGp.κSIE
    κSEIv = ones(N)*NGp.κSEI
    κSIIv = ones(N)*NGp.κSII
    κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]
    κS = weights(κSEEv, κSIEv, κSEIv, κSIIv, κSUM )
    wS = weightSave(κSEEv, κSIEv, κSEIv, κSIIv)
    bP = ballonModelParameters()
    opts=solverOpts(stimOpt,stimWindow,stimNodes,stimStr,Tstim,adapt,synapses,tWindows,nWindows)

    println("Running model ... ")
    @time out,weightSaved = NGModelRun(NGp,LR,bP,nP,κS,wS,opts,u0)

    BOLD_OUT=[]
    for ii = 1:nWindows
        if ii == 1
            BOLD_OUT = out[:,:,ii]
        else
            BOLD_OUT = cat(BOLD_OUT,out[:,:,ii],dims=2)
        end
    end
    

    savename = save1*save2
    dir0 = "LR_"
    dir1 = string(LR)
    dir2 = "_StimStr_"
    dir3 = string(stimStr)

    savedir = dir0*dir1*dir2*dir3

    

    println(fit)

    if save_data =="true"
        save("$OutDATADIR/$savedir/BOLD_$savename.jld","BOLD_$savename",BOLD_OUT)
        save("$OutDATADIR/$savedir/weights_$savename.jld","weights_$savename",weightSaved)
    end


    end

end

