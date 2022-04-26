using JLD

HOMEDIR = homedir()
DATADIR = "$HOMEDIR/NetworkModels/WilsonCowan_Distributed/data"
dataWC = load("$DATADIR/dataWC.jld","dataWC")

