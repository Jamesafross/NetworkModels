#includes 
using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,Random,NLsolve,Statistics,Parameters,Interpolations,BenchmarkTools,Profile

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
SC_Array,FC_Array,dist = getData_nonaveraged(;SCtype="log")

FC_Array = FC_Array

PaulFCmean = mean(FC_Array,dims=3)

SC = 0.01*SC_Array[:,:,1]


lags = dist./c

lags = round.(lags,digits=3) 
#lags[lags.<0.003] .= 0.000
#lags[SC .< 0.018] .= 0  
minSC,W_sum=getMinSC_and_Wsum(SC)
N = size(SC,1)

W = zeros(N,N)
W.=SC
NGp1 = get(ParSets,"Pset_2",1)
#NGp = NextGen2PopParams2(η_0E = -14.35)

bP = ballonModelParameters()
nWindows = 30
tWindows = 300.0

stimOpts = "off"
adapt = "off"
opts=modelOpts(stimOpts,adapt)

println("Running model ... ")
@unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp1

NGp = ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ 

nP = W,lags,N
   

       
u0 = zeros(8N)
u0[:] = makeInitConds(NGp1)

hparams = u0
      
        
tspan = (0.0,tWindows)
adpStops = collect(0.01:0.01:tWindows)
#println(adpTime)
clags = unique(reshape(lags[lags.>0.0],length(lags[lags.>0.0])))
println(clags)

h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(N)

p = (NGp,nP,hparams)


prob = DDEProblem(NextGen_test_,u0, h, tspan, p)
@time sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20);
