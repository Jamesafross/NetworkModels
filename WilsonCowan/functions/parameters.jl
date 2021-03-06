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

mutable struct networkParameters
    W::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
end

mutable struct modelOpts
    stimOpt::String
    adapt::String
end

mutable struct dataStruct
	modelR
	fit
	SC
end

mutable struct variousPars
    tPrev::Float64
    timeAdapt::Float64
end

	
