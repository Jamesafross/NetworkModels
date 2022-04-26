using LinearAlgebra,MAT,JLD,DifferentialEquations,Plots,StochasticDelayDiffEq,Random,NLsolve,Statistics,Parameters
HOMEDIR = homedir()
WORKDIR="$HOMEDIR/NetworkModels/2PopNextGen"
InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
include("./functions/stability.jl")
include("./functions/rhsFunctions.jl")
include("./functions/functions.jl")
include("./functions/parameters.jl")
include("../../Balloon_Model/balloonModelFunctions.jl")
include("../../Balloon_Model/balloonModelRHS.jl")
include("../../Balloon_Model/parameter_sets.jl")


normaliseSC = 0
#load data and make struct & dist matrices

#load data and make struct & dist matrices

SC = load("$InDATADIR/PaulSC.jld","C")
dist = load("$InDATADIR/PaulDist.jld","dist")

N = size(SC,1) # number of nodes
c =7000. # conductance velocity
lags = round.(dist./c,digits=1) # axonal delays
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
PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")



c =20000
lags = round.(dist./c,digits=3)
lags[lags .== 0.001] .= 0.0
nonzeros_indx = findall(SC.> 0.0)
clags = reshape(lags[lags.>0.0],length(lags[lags.>0.0])) # lags cant be zero for solver
W = zeros(N,N)
W.=SC

NGp = NextGen2PopParams()
W = W+(1/NGp.κ)diagm(ones(N))
Np = networkParameters(W,lags,N)
NoiseP = noiseParameters()

#set up initial conditions
@unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
params = κSEE,κSIE,κSEI,κSII,
αEE,αIE,αEI,αII,
κVEE,κVIE,κVEI,κVII,
VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

rE0, rI0, vE0, vI0, gEE0, gEI0, gIE0, gII0 = init_conds_SS(params)
perturb = 0.01*rand(8*N)
u0 = zeros(10*N)
u0 = zeros(8*N)
u0[1:N] .= rE0
u0[N+1:2N] .= rI0
u0[2N+1:3N] .= vE0
u0[3N+1:4N] .= vI0
u0[4N+1:5N] .= gEE0
u0[5N+1:6N] .= gIE0
u0[6N+1:7N] .= gEI0
u0[7N+1:8N] .= gII0

u0 = u0 + perturb

#time variables
Se = 1.0
T = 100*Se
tspan = (0.0,T)

#history
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? u0[idxs] :  u0
#FC mat init
R = zeros(N,N)
#solve!
p = NGp,Np,NoiseP
println("Running next-gen network")
prob = DDEProblem(NextGen,u0, h, tspan, p)
@time sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,saveat=0.05)


out_trans=(sol[1:N,100:end]')    
for n=1:N
   for m=n:N
        R[n,m]=(cor(out_trans[:,m],out_trans[:,n]))
        R[m,n]=R[n,m];
    end
end 




fit = fitR(R,PaulFCmean)

p1 = heatmap(R,c=:jet)
p2 = heatmap(PaulFCmean,c=:jet)

p = plot(p1,p2,layout=2)
p[:plot_title]  = "fit = $fit"
plot(p)