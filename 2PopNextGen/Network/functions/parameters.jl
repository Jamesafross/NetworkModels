using Parameters
@with_kw struct NextGen2PopParams{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = -14.78
    η_0I::R =-10.
    τE::R = 0.05
    τI::R = 0.06
    αEE::R = 10.2
    αIE::R = 10.5
    αEI::R = 10.6
    αII::R = 10.8
    κSEE::R = 5.
    κSIE::R = 3.3
    κSEI::R = 3.6
    κSII::R = 2.7
    κVEE::R = 0.
    κVIE::R = 0.
    κVEI::R = 0.
    κVII::R = 0.
    VsynEE::R = 4.0
    VsynIE::R = 1.5
    VsynEI::R = -0.5
    VsynII::R = -1.0
    κ::R = 0.10
end



@with_kw struct noiseParameters{R}
    τx::R = 0.1
    σE::R = 0.1
    σI::R = 0.1
end


mutable struct networkParameters
    W::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
end