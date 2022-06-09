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

c=7000.
SC_Array,FC_Array,dist = getData_nonaveraged(;SCtype="log")
FC_Array = FC_Array
PaulFCmean = mean(FC_Array,dims=3)
lags = dist./c
lags = round.(lags,digits=2) 
#lags[lags.<0.003] .= 0.000
#lags[SC .< 0.018] .= 0  
SC = 0.01SC_Array[:,:,1]
minSC,W_sum=getMinSC_and_Wsum(SC)
N = size(SC,1)
W = zeros(N,N)
W.=SC
const NGp = NextGen2PopParams2(η_0E = -14.19,κ=0.505)

# constants 


stimNodes = [39]
Tstim = [60,90]

#load data and make struct & dist matrices
nWindows = 30
tWindows = 100.0
stimOpt = "off"
stimStr = -5.
stimWindow = 20
adapt = "off"
synapses = "1stOrder"



κSEEv = ones(N)*NGp.κSEE
κSIEv = ones(N)*NGp.κSIE
κSEIv = ones(N)*NGp.κSEI
κSIIv = ones(N)*NGp.κSII
κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]

const nP = networkParameters(W, dist,lags, N, minSC,W_sum)
const bP = ballonModelParameters()
const LR = 0.01 # learning rate adaptation
const u0 = makeInitConds(NGp,N)
const κS = weights(κSEEv, κSIEv, κSEIv, κSIIv, κSUM )
const wS = weightSave(κSEEv, κSIEv, κSEIv, κSIIv)
const opts=solverOpts(stimOpt,stimWindow,stimNodes,stimStr,Tstim,adapt,synapses,tWindows,nWindows)
const vP = variousPars(0.0, 100.0,0)
const aP = adaptParams(10.01,u0[1:N])