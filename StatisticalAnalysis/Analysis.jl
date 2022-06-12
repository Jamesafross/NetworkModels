using JLD
using DelimitedFiles,Statistics,Plots

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels_Data/2PopNextGen_Data/"
OutDATADIR = "$HOMEDIR/NetworkModels_Data/Rcode_Data/2PopNextGen_Data"



function getFCwindows(SIG;type="r")
    FC = zeros(139,139)
    if type == "r"
        for ii = 1:139
            for jj = ii+1:139
                FC[ii,jj] = cor(SIG[ii,:],SIG[jj,:])
                FC[jj,ii] =  FC[ii,jj]
            end
        end
    
    elseif type == "phase"
            # Compute FC
            U_trans=imag(hilbert(SIG'))
            U_trans = U_trans'
            R=zeros(139,139);
            for n=1:139-1
                for m=n+1:139
                    FC[n,m]=abs((1/size(U_trans)[1]*sum(exp.(im*(U_trans[m,:]-U_trans[n,:])))));
                    FC[m,n]=FC[n,m];
                end
            end
    end
    return FC
end

counterT = 1
buffer = 100
LR = 0.01
stimStr = -5.
dir0 = "LR_"
dir1 = string(LR)
dir2 = "_StimStr_"
dir3 = string(stimStr)
savedir = dir0*dir1*dir2*dir3
run = 1

BOLD = load("$INDATADIR/$savedir/BOLD_NOstim_Run_$run.jld","BOLD_NOstim_Run_$run")[:,buffer:end]

step_i = 10
step_j = 300
for i = 1:step_i:size(BOLD,2)
    j = i + step_j
    if j < size(BOLD,2)
        counterT +=1
    else
        break
    end
end


counter = 1
TYPE = "r"


global FCstim = zeros(139,139,counterT) 
global FCnostim = zeros(139,139,counterT) 


runs = [1,2,3,4,5]
for run in runs
    for ii = 1

        BOLDstim = load("$INDATADIR/$savedir/BOLD_stim_Run_$run.jld","BOLD_stim_Run_$run")[:,buffer:end]
        BOLDnostim = load("$INDATADIR/$savedir/BOLD_NOstim_Run_$run.jld","BOLD_NOstim_Run_$run")[:,buffer:end]

        
        counter = 1
        for i = 1:step_i:size(BOLDstim,2)
            j = i + step_j
            if j < size(BOLDstim,2)
                global FCstim[:,:,counter] += getFCwindows(BOLDstim[:,i:j];type=TYPE)/length(runs)
                global FCnostim[:,:,counter] += getFCwindows(BOLDnostim[:,i:j];type=TYPE)/length(runs)
                counter += 1
                #println(counter)
                
            else
                break
            end
        end
    end
end

weights_stim = load("$INDATADIR/$savedir/weights_stim_Run_$run.jld","weights_stim_Run_$run")

weights_nostim = load("$INDATADIR/$savedir/weights_NOstim_Run_$run.jld","weights_NOstim_Run_$run")



anim1 = @animate for i = 1:1:counterT
    t = buffer*1.6 + round((i-1)*1.6*2,digits=2) + step_j*1.6 
    p1 = heatmap(FCstim[:,:,i].^2,c=:jet,title=t)
    p2 = heatmap(FCnostim[:,:,i].^2,c=:jet,title=t)
    p3 = heatmap(abs.(FCstim[:,:,i].^2 .-FCnostim[:,:,i].^2),c=:jet,title=t)
    plot(p1,p2,p3,title=i)
end

anim2 = @animate for i = 1:1:counterT
    t = buffer*1.6 + round((i-1)*1.6*2,digits=2) + step_j*1.6 
    p1 = heatmap(FCstim[:,:,i],c=:jet,title=t)
    p2 = heatmap(FCnostim[:,:,i],c=:jet,title=t)
    p3 = heatmap(abs.(FCstim[:,:,i] .-FCnostim[:,:,i]),c=:jet,title=t)
    plot(p1,p2,p3,title=i)
end

FCstim_rv = []
FCnostim_rv = []
rv=0
buff = 20
if rv == 1
    for i = buff+1:size(FCstim,3)
        j = i+buff
        jj = i-buff
        
        if j < size(FCstim,3)
            if i == buff+1
                FCstim_rv = mean(FCstim[:,:,jj:j],dims=3)
                FCnostim_rv = mean(FCnostim[:,:,jj:j],dims=3)
            else
                FCstim_rv = cat(FCstim_rv,mean(FCstim[:,:,jj:j],dims=3),dims=3)
                FCnostim_rv = cat(FCnostim_rv,mean(FCnostim[:,:,jj:j],dims=3),dims=3)
            end
        end
    end

    anim2 = @animate for i = 1:1:size(FCstim_rv,3)
        t = buffer*1.6 + round((i-1)*1.6*2,digits=2) + step_j*1.6 
        p1 = heatmap(FCstim_rv[:,:,i].^2,c=:jet,title=t)
        p2 = heatmap(FCnostim_rv[:,:,i].^2,c=:jet,title=t)
        p3 = heatmap(abs.(FCstim_rv[:,:,i].^2 .-FCnostim_rv[:,:,i].^2),c=:jet,title=t)
        plot(p1,p2,p3)
    end
end





gif(anim1,"anim1.gif",fps=20)




