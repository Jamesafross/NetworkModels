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
SC_Array,FC_Array,dist = getData_nonaveraged()

PaulMeanFC = mean(FC_Array,dims=3)[:,:]
SC = SC_Array[:,:,1]
lags = dist./c
minSC,W_sum=getMinSC_and_Wsum(SC)
N = size(SC,1)

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


FC_Array_stim= getFCstim_nonaverages()

PaulFCmean_stim = mean(FC_Array_stim,dims=3)[:,:]

fit = zeros(nWindows)
fit2_stim = zeros(nWindows)
fit2 = zeros(nWindows)

fitArray = complex(zeros(nWindows,size(FC_Array,3)))
fitArray_stim = complex(zeros(nWindows,size(FC_Array_stim,3)))
FC_R = zeros(size(FC_Array,3),size(FC_Array,3))
SCFCfit = zeros(nWindows)
FC_fit = zeros(size(FC_Array,3))
for i = 1:nWindows
    for j = 1:size(FC_Array,3)
        
        fitArray[i,j] = i + im*fitR(Rsave[:,:,i].^2,FC_Array[:,:,j].^2)
        if j <= size(FC_Array_stim,3)
        fitArray_stim[i,j] = i + im*fitR(Rsave[:,:,i].^2,FC_Array_stim[:,:,j].^2)
        end
    end
    fit[i] = fitR(Rsave[:,:,i],PaulFCmean)
    fit2[i] = fitR((Rsave[:,:,i].^2),(PaulFCmean_stim.^2))
    fit2_stim[i] = fitR((Rsave[:,:,i].^2),(PaulFCmean.^2))
    SCFCfit[i] = fitR(SC*1,Rsave[:,:,i])
end

for i = 1:size(FC_Array,3)
    FC_fit[i] = fitR(PaulFCmean.^2,FC_Array[:,:,i].^2)
    for j = 1:size(FC_Array,3)
    FC_R[i,j] = fitR(FC_Array[:,:,i].^2,FC_Array[:,:,j].^2)
    end
end

if nWindows > 1
    meanR = mean(Rsave[:,:,:],dims=3)
    meanfit = fitR(meanR,PaulFCmean)
    meanfit2 = fitR(meanR.^2,PaulFCmean.^2)
else
    scatter(collect(1:1:nWindows),fit2)
    meanfit = round.(fit,digits=4)
    meanR = Rsave
    meanfit2=round.(fit2,digits=4)
end


p1 = heatmap(meanR[:,:,1],c=:jet)
p2 = heatmap(PaulFCmean,c=:jet)
p3 = heatmap((meanR[:,:,1].^2),c=:jet)
p4 = heatmap((PaulFCmean.^2),c=:jet)

if nWindows > 1
    p5 = scatter(collect(1:1:nWindows),fit2)
    p = plot(p1,p2,p3,p4,layout=4)
else
    p = plot(p1,p2,p3,p4,layout=4)
end

fitArray = reshape(fitArray,nWindows*28)
fitArray_stim = reshape(fitArray_stim,nWindows*25)

p[:plot_title]  = "fit (r) = $meanfit,  (r² ) = $meanfit2 "
plot(p)

#scatter(collect(1:1:nWindows),SCFCfit,label="SC fit")
scatter(collect(1:1:nWindows),fit2,label="FC fit (mean)")
scatterplot1 = scatter!(collect(1:1:nWindows),fit2_stim,label="FC stim fit (mean)")

p[:plot_title]  = "fit (r) = $meanfit,  (r² ) = $meanfit2 "
plot(p)

scatter(collect(1:1:nWindows),fit2,label="FC fit (mean)")
scatterplot1 = scatter!(collect(1:1:nWindows),fit2_stim,label="FC stim fit (mean)")

scatter(collect(1:1:nWindows),SCFCfit,label="SC fit")
scatter!(collect(1:1:nWindows),fit2,label="FC fit (mean)")

scatterplot2 = scatter!(real(fitArray),imag(fitArray),label="FC fit (windows)")