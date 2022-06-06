using JLD
using DelimitedFiles,Statistics,Plots,Hilbert

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels_Data/2PopNextGen_Data/"
OutDATADIR = "$HOMEDIR/NetworkModels_Data/Rcode_Data/2PopNextGen_Data"


nWindows = 22



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
LR = 0.00001
stimStr = 10.
dir0 = "LR_"
dir1 = string(LR)
dir2 = "_StimStr_"
dir3 = string(stimStr)
savedir = dir0*dir1*dir2*dir3

BOLD = load("$INDATADIR/$savedir/BOLD_NOstimAdaptivity.jld","BOLD_NOstimAdaptivity")[:,buffer:end]

step_i = 5
step_j = 60
for i = 1:step_i:size(BOLD,2)
    j = i + step_j
    if j < size(BOLD,2)
        counterT +=1
    else
        break
    end
end

runs=1
counter = 1
TYPE = "r"



global FCstim = zeros(139,139,counterT) 
global FCnostim = zeros(139,139,counterT) 

for ii = 1

    BOLDstim = load("$INDATADIR/$savedir/BOLD_stimAdaptivity.jld","BOLD_stimAdaptivity")[:,buffer:end]
    BOLDnostim = load("$INDATADIR/$savedir/BOLD_NOstimAdaptivity.jld","BOLD_NOstimAdaptivity")[:,buffer:end]

    
    counter = 1
    for i = 1:step_i:size(BOLDstim,2)
        j = i + step_j
        if j < size(BOLDstim,2)
            global FCstim[:,:,counter] += getFCwindows(BOLDstim[:,i:j];type=TYPE)/runs
            global FCnostim[:,:,counter] += getFCwindows(BOLDnostim[:,i:j];type=TYPE)/runs
            counter += 1
            #println(counter)
            
        else
            break
        end
    end
end


@gif for i = 1:1:counterT
    t = 160 + round((i-1)*1.6,digits=2) + step_j*1.6
    p1 = heatmap(FCstim[:,:,i].^2,c=:jet,title=t)
    p2 = heatmap(FCnostim[:,:,i].^2,c=:jet,title=t)
    p3 = heatmap(abs.(FCstim[:,:,i].^2 .-FCnostim[:,:,i].^2),c=:jet,title=t)
    plot(p1,p2,p3,title=i)
end
    








