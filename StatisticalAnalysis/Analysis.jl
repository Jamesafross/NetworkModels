using JLD
using DelimitedFiles,Statistics,Plots

HOMEDIR = homedir()
INDATADIR = "$HOMEDIR/NetworkModels_Data/2PopNextGen_Data"
OUTDATADIR = "$HOMEDIR/NetworkModels_Data/Rcode_Data/2PopNextGen_Data"

runs = 2

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

i1 = 1
i2 = 10

for i = 1:16
    SD_nostim[i] = std(FC_WINDOWS_nostim[i1,i2,i,:])
    SD_stim[i] = std(FC_WINDOWS_stim[i1,i2,i,:])
end

plot(mean(FC_WINDOWS_nostim[i1,i2,:,:],dims=2),ribbon = SD_nostim,color="red",fillalpha=0.2,alpha=1.0,legend=true,ylims=[-1,1],label="no stim")
plot!(mean(FC_WINDOWS_stim[i1,i2,:,:],dims=2),ribbon = SD_stim,color="blue",fillalpha=0.2,alpha=1.0,legend=true,ylims=[-1,1],label="stim")
#plot(FC_Array[i1,i2,:],color="red",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])
#plot!(FC_Array_stim[i1,i2,:],color="blue",fillalpha=0.2,alpha=1.0,legend=false,ylims=[-1,1])

plot(abs.(mean(FC_WINDOWS_nostim[i1,i2,:,:].-FC_WINDOWS_stim[i1,i2,:,:],dims=2)),ribbon = SD_nostim,color="red",fillalpha=0.2,alpha=1.0,legend=true,ylims=[-1,1],label="no stim")


stimWindow = 12
AbsDiff = zeros(139,139,16-stimWindow)
for i = stimWindow+1:16
    j = i-stimWindow
    println(j)
    AbsDiff[:,:,j] = abs.(mean(FC_WINDOWS_nostim[:,:,i,:].-FC_WINDOWS_stim[:,:,i,:],dims=3))[:,:]
end

mAbsDiff = mean(AbsDiff,dims=3)
mAbsDiff[39,39] = NaN



