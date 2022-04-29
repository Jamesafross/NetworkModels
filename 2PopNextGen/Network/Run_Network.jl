#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters,Interpolations

HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
include("./functions/NextGenFunctions.jl")
include("../../Balloon_Model/BalloonModel.jl")
include("$InDATADIR/getData.jl")


stimNodes = [21,39]
Tstim = [60,90]


#load data and make struct & dist matrices
c=7000.
SC,minSC,W_sum,lags,PaulFCmean,N = getData(c;normalise=0,delayDigits=3)
lags[lags .<= 0.005] .= 0
SC = SC[1:size(SC,1) .!= 100,1:size(SC,1) .!= 100 ]
W_sum = W_sum[1:size(W_sum,1) .!= 100]
N = size(SC,1)
PaulFCmean = PaulFCmean[1:size(PaulFCmean,1) .!= 100,1:size(PaulFCmean,1) .!= 100 ]
PaulFCmean = PaulFCmean .- diagm(ones(N))

clags = reshape(lags[lags.>0.0],length(lags[lags.>0.0])) # lags cant be zero for solver
W = zeros(N,N)
W.=SC
NGp = NextGen2PopParams()

bP = ballonModelParameters()
nWindows = 1
tWindows = 300.0

stimOpts = "off"
adapt = "off"
opts=modelOpts(stimOpts,adapt)

println("Running model ... ")
@time Rsave,Wsave = NGModelRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)


fit = zeros(nWindows)
fit2 = zeros(nWindows)
for i = 1:nWindows
    fit[i] = fitR(Rsave[:,:,i],PaulFCmean)
    fit2[i] = fitR((Rsave[:,:,i].^2)./maximum(Rsave[:,:,i].^2),(PaulFCmean.^2)./maximum(PaulFCmean.^2))
end

if nWindows > 1
    meanR = mean(Rsave[:,:,:],dims=3)
    meanfit = fitR(meanR,PaulFCmean)
    meanfit2 = fitR(meanR.^2,PaulFCmean.^2)
else
    meanfit = fit
    meanR = Rsave
    meanfit2=fit2
end


p1 = heatmap(meanR[:,:,1],c=:jet)
p2 = heatmap(PaulFCmean,c=:jet)
p3 = heatmap((meanR[:,:,1].^2)./maximum(meanR.^2),c=:jet)
p4 = heatmap((PaulFCmean.^2)./maximum(PaulFCmean.^2),c=:jet)
p6 = heatmap((SC)/maximum(SC),c=:jet)
if nWindows > 1
    p5 = scatter(collect(1:1:nWindows),fit2)
    p = plot(p1,p2,p3,p4,p5,p6,layout=6)
else
    p = plot(p1,p2,p3,p4,p6,layout=5)
end


p[:plot_title]  = "fit = $meanfit, fit2 = $meanfit2 "
plot(p)