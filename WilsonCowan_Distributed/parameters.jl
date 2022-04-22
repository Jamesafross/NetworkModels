@with_kw struct WCparams{R}
    cEE::R = 3.5
    cEI::R  = 3.75
    cIE::R  = -2.5
    cII::R  = 0.
    θE::R  = 1.
    θI::R  = 1.
    β::R = 4.
    η::R  = 0.14
    σ::R  = 0.001
    τE::R  = 0.01
    τI::R  =0.02
    τx::R = 0.1
    Pext::R  =0.31

end

mutable struct networkParameters
    W::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
end

