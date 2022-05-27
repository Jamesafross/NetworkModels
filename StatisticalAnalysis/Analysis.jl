using JLD
using DelimitedFiles,Statistics,Plots

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels_Data/2PopNextGen_Data/2PopNextGen_Data"
OUTDATADIR = "$HOMEDIR/NetworkModels_Data/Rcode_Data/2PopNextGen_Data"

BOLD = save("$OutDATADIR/Run_$Run/BOLD_$savename.jld","BOLD_$savename",BOLD_OUT)
nWindows = 16
for i = 1:m
    j = i + 4
 
end


function getFCwindows(nWindows,Run)
    FC = zeros(139,139,nWindows)
    BOLD = save("$INDATADIR/Run_$Run/BOLD_StimAdaptivity.jld","BOLD_StimAdaptivity")
    for i = 1:nWindows
        for ii = 1:139
            for jj = 1:139
                FC[ii,jj,nWindows] = cor(BOLD[ii,:],BOLD[jj,:])
            end
        end
    end
end

runs = 2

for i = 




