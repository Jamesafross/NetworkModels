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


BOLD_stim = load("$INDATADIR/$savedir/BOLD_stimAdaptivity.jld","BOLD_stimAdaptivity")[:,buffer:end]
BOLD_nostim = load("$INDATADIR/$savedir/BOLD_NOstimAdaptivity.jld","BOLD_NOstimAdaptivity")[:,buffer:end]

diffBold = BOLD_stim .- BOLD_nostim

stimIdx = findfirst(x->x>0,diffBold)[2]

