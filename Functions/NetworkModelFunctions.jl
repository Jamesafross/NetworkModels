module NetworkModelFunctions
using Parameters,Distributed
export make_uhist,adapt_global_coupling,stim,h1,h2,normalise,fitR
export makeInitConds,NextGenRun
export f,WCRun
export WCISPRun,WCISPDistributedLoop
include("GlobalFunctions.jl")
include("NextGen.jl")
include("WilsonCowan.jl")
include("WilsonCowanISP.jl")
end