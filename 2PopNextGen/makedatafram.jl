using StatsPlots, RDatasets
using DataFrames, IndexedTables

HOMEDIR = homedir()
DATADIR="$HOMEDIR/NetworkModels/2PopNextGen/data"
dataNoStimFC = load("$DATADIR/dataSave_NOstimAdaptivity.jld","data_save_NOstimAdaptivity").modelR
dataStimFC = load("$DATADIR/dataSave_stimAdaptivity.jld","data_save_stimAdaptivity").modelR



 MM = DataFrame(ROI = "[39,20]", a = dataNoStimFC[39,20,:],b = dataStimFC[39,20,:])

diff = dataStimFC.- dataNoStimFC





meanDiff = mean(diff[:,:,1:50],dims=3)[:,:]
meanDiff[abs.(meanDiff).<0.65] .= 0
meanDiff[abs.(meanDiff).>= 0.65] .= 1
stdDiff = std(diff,dims=3)[:,:]

sumVec = zeros(139)
for i = 1:139
    sumVec[i] = sum(meanDiff[i,:])
end





