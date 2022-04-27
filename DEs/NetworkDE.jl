module NetworkDE
using Parameters
export NextGen, WC, WC_ISP,dW, dW_ISP

include("NextGen.jl")
include("WilsonCowan.jl")
include("WilsonCowanISP.jl")

end