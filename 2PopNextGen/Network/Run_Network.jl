#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations,MKL

HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels_Data/StructDistMatrices"
OutDATADIR="$HOMEDIR/NetworkModels_Data/2PopNextGen_Data"
include("./functions/NextGenFunctions.jl")
include("../../Balloon_Model/BalloonModel.jl")
include("$InDATADIR/getData.jl")

Run_vec = LinRange(1,20,20)
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
    tWindows = 300.0
    stimOpt = setstim
    stimSri = 5.
    stimWindow = 5
    adapt = "on"
    synapses = "1stOrder"



    κSEEv = ones(N)*NGp.κSEE
    κSIEv = ones(N)*NGp.κSIE
    κSEIv = ones(N)*NGp.κSEI
    κSIIv = ones(N)*NGp.κSII
    κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]
    κS = weights(κSEEv, κSIEv, κSEIv, κSIIv, κSUM )
    bP = ballonModelParameters()
    opts=solverOpts(stimOpt,stimWindow,stimNodes,stimStr,Tstim,adapt,synapses,tWindows,nWindows)

    println("Running model ... ")
    @time Rsave,Wsave,out = NGModelRun(NGp,LR,bP,nP,κS,opts,u0)

    BOLD_OUT=[]
    for ii = 1:nWindows
        if ii == 1
            BOLD_OUT = out[:,:,ii]
        else
            BOLD_OUT = cat(BOLD_OUT,out[:,:,ii],dims=2)
        end
    end
    FC_Array_stim= getFCstim_nonaverages()

    PaulFCmean_stim = mean(FC_Array_stim,dims=3)[:,:]

    fit = zeros(nWindows)
    fit2 = zeros(nWindows)

    fit_stim = zeros(nWindows)
    fit2_stim=zeros(nWindows)

    fitt,bestIdxs = findBestFit(Rsave[:,:,1],FC_Array)

    fitt,bestIdxs2 = findBestFit(Rsave[:,:,1].^2,FC_Array.^2)

    fitt,bestIdxs_stim = findBestFit(Rsave[:,:,1],FC_Array_stim)

    fitt,bestIdxs2_stim = findBestFit(Rsave[:,:,1].^2,FC_Array_stim.^2)


    for i = 1:nWindows

        fit[i] = fitR(Rsave[:,:,i],mean(FC_Array[:,:,bestIdxs],dims=3))
        fit2[i] = fitR(Rsave[:,:,i].^2,mean(FC_Array[:,:,bestIdxs2].^2,dims=3))
        fit_stim[i] = fitR(Rsave[:,:,i],mean(FC_Array_stim[:,:,bestIdxs_stim],dims=3))
        fit2_stim[i] = fitR(Rsave[:,:,i].^2,mean(FC_Array_stim[:,:,bestIdxs2_stim].^2,dims=3))
    end

    fitAll = [[fit],[fit2],[fit_stim],[fit2_stim]]
    dataSave = dataStruct(Rsave,fitAll,Wsave)

    if stimOpt == "on"
        save1 = "stim"
    else
        save1="NOstim"
    end
    if adapt == "on"
        save2 = "Adaptivity"
    else
        save2="NOadaptivity"
    end

    savename = save1*save2
    dir0 = "LR_"
    dir1 = string(LR)
    dir2 = "_StimStr_"
    dir3 = string(stimStr)

    savedir = dir0*dir1*dir2*dir3

    

    println(fit)

    if save_data =="true"
        save("$OutDATADIR/$savedir/dataSave_$savename.jld","data_save_$savename",dataSave)
        save("$OutDATADIR/$savedir/BOLD_$savename.jld","BOLD_$savename",BOLD_OUT)
    end

    end

end


if plot_fit == "true"
scatter(collect(1:1:nWindows),fit,label="FC fit ")

scatter!(collect(1:1:nWindows),fit2,label="FC fit2 ")

scatter!(collect(1:1:nWindows),fit_stim,label="FC stim fit ")
scatter!(collect(1:1:nWindows),fit2_stim,label="FC stim fit2 ")

plot!(collect(1:1:nWindows),fit,label="FC fit ")
plot!(collect(1:1:nWindows),fit2,label="FC fit2 ")
plot!(collect(1:1:nWindows),fit_stim,label="FC stim fit ")

scatterplot1 =  plot!(collect(1:1:nWindows),fit2_stim,label="FC stim fit2 ")
end