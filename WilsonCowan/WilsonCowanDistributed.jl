parallel = "on"
if parallel == "on"
    include("setupDistributed.jl")
else
    include("setup.jl")
end



# get parameters and make structures
if opts.ISP == "on"
    WCparamsISP()
else
    WCp= WCparams()
end
bP = ballonModelParameters()

# Stimulation Setup
nWindows = 20
tWindows = 300.0
nTrials = 1
fitArray = zeros(nWindows,nTrials)
Rsave= SharedArray(zeros(N,N,nWindows,nTrials))
W_save = SharedArray(zeros(N,N,nWindows,nTrials))
W = zeros(N,N)
W .= SC
stimOpts = "off"
adapt = "on"
opts=modelOpts(stimOpts,adapt)


if parallel == "on"
@sync @distributed for i = 1:nTrials
    println("working on Trial: ",i)
    Rsave[:,:,:,i],W_save[:,:,:,i] = WCRun(WCp,bP,opts)
end
else
    for i = 1:nTtrials
    println("working on Trial: ",i)
    Rsave[:,:,:,i],W_save[:,:,:,i] = WCRun(WCp,bP,opts)
    end
end


if nTrials == 1
    Rsave[:,:,:,1] = Rsave[:,:,:,1][:,:,:]
else
    Rsave = mean(Rsave,dims=4)[:,:,:]
end

FC_Array = getFC_nonaverages()
FC_Array_stim= getFCstim_nonaverages()
PaulFCmean = mean(FC_Array,dims=3)[:,:]
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
scatterplot = scatter!(collect(1:1:nWindows),fit2_stim,label="FC stim fit (mean)")

#scatterplot = scatter!(real(fitArray_stim),imag(fitArray_stim),label="FC stim fit  (windows)")
scatter(collect(1:1:nWindows),SCFCfit,label="SC fit")
scatter!(collect(1:1:nWindows),fit2,label="FC fit (mean)")

scatterplot2 = scatter!(real(fitArray),imag(fitArray),label="FC fit (windows)")


    
