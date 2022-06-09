using JLD
using DelimitedFiles,Statistics,Plots,XLSX

HOMEDIR = homedir()
WORKDIR = "$HOMEDIR/NetworkModels/StatisticalAnalysis"

INDATADIR = "$HOMEDIR/NetworkModels_Data/2PopNextGen_Data"
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

BOLD_stim = BOLD_stim[:,stimIdx:end]
BOLD_nostim = BOLD_nostim[:,stimIdx:end]
size1 = size(BOLD_stim,2)

time = collect(0:1.6:size1*1.6)./60


ROI = string.(readdlm("$WORKDIR/ROI.txt"))
C = Matrix{String}(undef,size1,1)
fill!(C,"Control")
D = Matrix{String}(undef,size1,1)
fill!(D,"FUS")
Run = Matrix{String}(undef,size1,1)
fill!(Run,"Run1")

sub = Matrix{String}(undef,size1,1)
fill!(sub,"Model1")

touch("$WORKDIR/modelBOLD.xlsx")

XLSX.openxlsx("$WORKDIR/modelBOLD.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "Suj"
    sheet["B1"] = "Group"
    sheet["C1"] = "Run"
    sheet["D1"] = "Time"
    sheet["E1:EM1"] = ROI
    sheet["D2",dim=1] = zeros(1285)
    sheet["C2:C$(size1+1)"] = Run
    sheet["A2:A$(size1+1)"] = sub
    sheet["B2:B$(size1+1)"] = C
    sheet["E2:EM$(size1+1)"] = Matrix(BOLD_nostim[:,:]')

    sheet["D1287",dim=1] = time
    sheet["C$(size1+2):C$(2*size1+1)"] = Run
    sheet["A$(size1+2):A$(2*size1+1)"] = sub
    sheet["B$(size1+2):B$(2*size1+1)"] = D
    sheet["E$(size1+2):EM$(2*size1+1)"] = Matrix(BOLD_stim[:,:]')

end
