#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations

HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
include("./functions/NextGenFunctions.jl")
include("../../Balloon_Model/BalloonModel.jl")
include("$InDATADIR/getData.jl")

Run_vec = LinRange(6,10,5)
plot_fit = "false"


for jj = 1:length(Run_vec)
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
minSC,W_sum=getMinSC_and_Wsum(SC)
N = size(SC,1)
W = zeros(N,N)
W.=SC


Run = string(Int(round(Run_vec[jj])))
nWindows = 16
tWindows = 300.0
stimOpt = "on"
stimWindow = 2
adapt = "on"
synapses = "1stOrder"

NGp = get(ParSets,"Pset_2",1)
NGp = NextGen2PopParams2(η_0E = -14.19,κ=0.505)
κSEEv = ones(N)*NGp.κSEE
κSIEv = ones(N)*NGp.κSIE
κSEIv = ones(N)*NGp.κSEI
κSIIv = ones(N)*NGp.κSII
κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]
κS = weights(κSEEv, κSIEv, κSEIv, κSIIv, κSUM )
bP = ballonModelParameters()
opts=solverOpts(stimOpt,stimWindow,stimNodes,Tstim,adapt,synapses,tWindows,nWindows)
nP = networkParameters(W, dist,lags, N,minSC,W_sum)

println("Running model ... ")
@time Rsave,Wsave,out = NGModelRun(NGp,bP,nP,κS,opts)


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


#println("fit = ", fit[:])
#scatter(collect(1:1:nWindows),SCFCfit,label="SC fit")
#scatter!(collect(1:1:nWindows),fitκ2,label="FC (r²) fit (mean)")
#scatter!(collect(1:1:nWindows),fit,label="FC (r) fit (mean)")
#scatterplot2 = scatter!(real(fitArray),imag(fitArray),label="FC fit (windows)")
#heatmap(Rsave[:,:,1])

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

println(fit)


save("$HOMEDIR/NetworkModels/2PopNextGen/data/Run_$Run/dataSave_$savename.jld","data_save_$savename",dataSave)


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
en1