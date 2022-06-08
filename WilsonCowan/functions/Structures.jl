@with_kw struct WCparams{R}
    cEE::R = 3.5
    cEI::R  = 3.75
    cIE::R  = -2.5
    cII::R  = 0.
    θE::R  = 1.
    θI::R  = 1.
    β::R = 4.
    η::R  = 0.12
    σ::R  = 0.005
    τE::R  = 0.01
    τI::R  =0.02
    τx::R = .01
    Pext::R  =0.315
end

@with_kw struct WCparamsISP{R}
    cEE::R = 3.5
    cEI::R  = 3.75
    cIE::R  = -2.5
    cII::R  = 0.
    θE::R  = 1.
    θI::R  = 1.
    β::R = 4.0
    η::R  = 0.08
    σ::R  = 0.005
    τE::R  = 0.01
    τI::R  =0.02
    τx::R = .01
    Pext::R  =0.31
    τISP::R = 20.
    ρ::R = 0.15
end

mutable struct networkParameters
    W::Matrix{Float64}
    dist::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
    minSC::Float64
    W_sum::Vector{Float64}
end


mutable struct dataStruct
	modelR
	fit
	SC
end

mutable struct variousPars
    tPrev::Real
    timeAdapt::Real
end

	
mutable struct solverOpts
    stimOpt::String
    stimWindow::Real
    stimNodes::Vector{Real}
    Tstim::Vector{Real}
    adapt::String
    tWindows::Real
    nWindows::Real
    ISP::String
end