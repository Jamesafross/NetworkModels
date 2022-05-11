#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters,Interpolations

HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
include("./functions/NextGenFunctions.jl")
include("../../Balloon_Model/BalloonModel.jl")
include("$InDATADIR/getData.jl")
include("$InDATADIR/getDataNimg.jl")


stimNodes = [21,39]
Tstim = [60,90]



#load data and make struct & dist matrices
c=7000.
SC_Array,FC_Array_AL,FC_Array_DP,FC_Array_FL,FC_Array_VL,Dist_Array= getDataNimg(;normalise=0,delayDigits=2,SCtype="log")

dist = Dist_Array[:,:,1]
FC_Array = FC_Array_AL

PaulFCmean = mean(FC_Array,dims=3)

SC = 0.01*(SC_Array[:,:,1] .- diagm(diag(SC_Array[:,:,1])))


lags = dist./c

lags = round.(lags,digits=2)
#lags[lags.>0.005] .= 0.005
#lags[0.0.<lags.<0.0040] .=0
minSC,W_sum=getMinSC_and_Wsum(SC)
N = size(SC,1)

clags = reshape(lags[lags.>0.0],length(lags[lags.>0.0])) # lags cant be zero for solver
W = zeros(N,N)
W.=SC
NGp = get(ParSets,"Pset_2",1)
#NGp = NextGen2PopParams2(η_0E = -13.8)

bP = ballonModelParameters()
nWindows = 20
tWindows = 300.0

stimOpts = "off"
adapt = "off"
opts=modelOpts(stimOpts,adapt)

println("Running model ... ")
@time Rsave,Wsave = NGModelRun(NGp,bP,nWindows,tWindows,W,lags,dist,N,minSC,W_sum,opts)





fit = zeros(nWindows)
fit2 = zeros(nWindows)



fitt,bestIdxs = findBestFit(Rsave[:,:,1],FC_Array)

fitt,bestIdxs2 = findBestFit(Rsave[:,:,1].^2,FC_Array.^2)


for i = 1:nWindows

    fit[i] = fitR(Rsave[:,:,i],mean(FC_Array[:,:,bestIdxs],dims=3))
    fit2[i] = fitR(Rsave[:,:,i].^2,mean(FC_Array[:,:,bestIdxs2].^2,dims=3))
   
end

scatter(collect(1:1:nWindows),fit,label="FC fit ")

scatter!(collect(1:1:nWindows),fit2,label="FC fit2 ")


plot!(collect(1:1:nWindows),fit,label="FC fit ")


scatterplot1 = plot!(collect(1:1:nWindows),fit2,label="FC fit2 ")
#println("fit = ", fit[:])
#scatter(collect(1:1:nWindows),SCFCfit,label="SC fit")
#scatter!(collect(1:1:nWindows),fit2,label="FC (r²) fit (mean)")
#scatter!(collect(1:1:nWindows),fit,label="FC (r) fit (mean)")
#scatterplot2 = scatter!(real(fitArray),imag(fitArray),label="FC fit (windows)")
#heatmap(Rsave[:,:,1])