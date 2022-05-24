using JLD
using DelimitedFiles,Statistics,Plots

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels/2PopNextGen/data"
OUTDATADIR = "$HOMEDIR/NetworkModels/Rcode/data_nextgen"

runs = 5

FC_WINDOWS_nostim = zeros(139,139,16,runs)
FC_WINDOWS_stim = zeros(139,139,16,runs)
SD_nostim = zeros(16)
SD_stim = zeros(16)

SC_Array,FC_Array,dist = getData_nonaveraged(;SCtype="log")
FC_Array_stim = getFCstim_nonaverages()


ROIs = Matrix{String}(undef, 1,139)

for i = 1:139
    ROIs[i] = "X$i"
end


# control runs
for i = 1:runs
    FC_WINDOWS_nostim[:,:,:,i] = load("$INDATADIR/Run_$i/dataSave_NOstimAdaptivity.jld","data_save_NOstimAdaptivity").modelR
    FC_WINDOWS_stim[:,:,:,i] = load("$INDATADIR/Run_$i/dataSave_stimAdaptivity.jld","data_save_stimAdaptivity").modelR
end

i1 = rand(1:139)
i2 = 39

for i = 1:16
    SD_nostim[i] = std(FC_WINDOWS_nostim[i1,i2,i,:])
    SD_stim[i] = std(FC_WINDOWS_stim[i1,i2,i,:])
end

plot(mean(FC_WINDOWS_nostim[i1,i2,:,:],dims=2),ribbon = SD_nostim,color="red",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])
plot!(mean(FC_WINDOWS_stim[i1,i2,:,:],dims=2),ribbon = SD_stim,color="blue",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])
#plot(FC_Array[i1,i2,:],color="red",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])
#plot!(FC_Array_stim[i1,i2,:],color="blue",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])






