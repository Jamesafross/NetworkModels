#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters
HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
include("./functions/stability.jl")
include("./functions/DEfunctions.jl")
include("./functions/functions.jl")
include("./functions/parameters.jl")
include("./modelFunc.jl")
include("../../Balloon_Model/balloonModelFunctions.jl")
include("../../Balloon_Model/balloonModelRHS.jl")
include("../../Balloon_Model/parameter_sets.jl")


normaliseSC = 0
#load data and make struct & dist matrices
SC = load("$InDATADIR/PaulSC.jld","C")
dist = load("$InDATADIR/PaulDist.jld","dist")
PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")

N = size(SC,1) # number of nodes
c =7000. # conductance velocity
lags = round.(dist./c,digits=2) # axonal delays
stimNodes = [21,39]
Tstim = [60,90]
if normaliseSC == 1
    SC = normalise(SC,N)
    W_sum = ones(N)
else
    W_sum = zeros(N)
    for i = 1:N
        W_sum[i] = sum(SC[:,i])
    end 
end
minSC = minimum(SC[SC.>0.0])




c =7000
lags = round.(dist./c,digits=2)
lags[lags .== 0.001] .= 0.0
nonzeros_indx = findall(SC.> 0.0)
clags = reshape(lags[lags.>0.0],length(lags[lags.>0.0])) # lags cant be zero for solver
W = zeros(N,N)
W.=SC

NGp = NextGen2PopParams()
W = W+(1/NGp.κ)diagm(ones(N))
bP = ballonModelParameters()
nWindows = 1
tWindows = 300

stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)
nTrials = 10
κVec = LinRange(0.08,0.14,nTrials)
η_0EVec = LinRange(-14.0,-14.5,nTrials)

Rsave = zeros(N,N,nWindows,nTrials,nTrials)
for i = 1:nTrials
    for j = 1:nTrials
        NGp = NextGen2PopParams(κ=κVec[i],η_0E=η_0EVec[j])
        Rsave[:,:,:,i,j],Wsave = NGModelRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
    end
end



fit = zeros(nWindows)

if nWindows > 1
    meanR = mean(Rsave[:,:,:,:,:],dims=3)[:,:,1,:,:]
    meanfit = fitR(meanR,PaulFCmean)
else
    meanfit = fit
    meanR = Rsave[:,:,1,:,:]
end

pSweep = []
for i = 1:nTrials
    for j = 1:nTrials
        
        fitp = fitR(meanR[:,:,i,j],PaulFCmean)
        pSweep = cat(pSweep,pSweepData(κVec[i],η_0EVec[j],fitp),dims=1)
    end
end


p1 = heatmap(meanR[:,:,1],c=:jet)
p2 = heatmap(PaulFCmean,c=:jet)

if nWindows > 1
    p3 = scatter(collect(1:1:nWindows),fit)
    p = plot(p1,p2,p3,layout=3)
else
    p = plot(p1,p2,layout=2)
end


p[:plot_title]  = "fit = $meanfit"
plot(p)