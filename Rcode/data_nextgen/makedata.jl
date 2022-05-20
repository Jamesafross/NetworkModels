using JLD
using DelimitedFiles

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels/2PopNextGen/data"
OUTDATADIR = "$HOMEDIR/NetworkModels/Rcode/data_nextgen"

runs = 2

FC_WINDOWS = []

ROIs = Matrix{String}(undef, 1,139)

for i = 1:139
    ROIs[i] = "X$i"
end


# control runs
for i = 1:runs
    FC_WINDOWS = load("$INDATADIR/Run_$i/dataSave_NOstimAdaptivity.jld","data_save_NOstimAdaptivity").modelR
    
    for j = 1:size(FC_WINDOWS,3)
        
        open("$OUTDATADIR/NO_STIM/control_window$j.txt", "w") do io
            writedlm(io,cat(ROIs,FC_WINDOWS[:,:,j],dims=1))
        end
    end
    
end
