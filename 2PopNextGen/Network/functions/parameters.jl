using Parameters

@with_kw mutable struct NextGen2PopParams{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = 1.
    η_0I::R =-2.
    τE::R = 0.05
    τI::R = 0.06
    αEE::R = 10.2
    αIE::R = 10.5
    αEI::R = 10.6
    αII::R = 10.1
    κSEE::R = 5.2
    κSIE::R = 3.6
    κSEI::R = 3.6
    κSII::R = 5.7
    κVEE::R = 0.
    κVIE::R = 0.
    κVEI::R = 0.
    κVII::R = 0.
    VsynEE::R = 2.0
    VsynIE::R = 1.5
    VsynEI::R = -1.5
    VsynII::R = -5.0
    κ::R = 0.1
end

@with_kw mutable struct NextGen2PopParams2{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = -14.19
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
    κ::R = 0.505
end
NGp1 = NextGen2PopParams()
NGp2 = NextGen2PopParams2()


ParSets = Dict("Pset_1"=>NGp1,"Pset_2"=>NGp2)

@with_kw struct noiseParameters{R}
    τx::R = 0.1
    σE::R = 0.1
    σI::R = 0.1
end


mutable struct networkParameters
    W::Matrix{Float64}
    dist::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
    minSC::Float64
    W_sum::Vector{Float64}
end


mutable struct solverOpts
    stimOpt::String
    stimWindow::Real
    stimNodes::Vector{Real}
    stimStr::Real
    Tstim::Vector{Real}
    adapt::String
    synapses::String
    tWindows::Real
    nWindows::Real
end

mutable struct dataStruct
	modelR
	fit
	SC
end

mutable struct pSweepData
    c
	κ
	ηE
	fit
end

mutable struct variousPars
    tPrev::Float64
    timeAdapt::Float64
    count::Int64
end

mutable struct adaptParams
    tP::Float64
    HIST::Array{Float64}
end


mutable struct weights
    κSEEv
    κSIEv
    κSEIv
    κSIIv
    κSUM
end


@with_kw mutable struct NextGen2PopParams3{R}
    ΔE::R = 0.5
    ΔI::R = 0.5
    η_0E::R = 0.1
    η_0I::R =0.1
    τE::R = 0.2
    τI::R = 0.3
    αEE::R = 10.
    αIE::R =10.0
    αEI::R = 5.0
    αII::R = 5.0
    κSEE::R = 0.4
    κSIE::R = 0.3
    κSEI::R = 0.1
    κSII::R = 0.1
    κVEE::R =0.01
    κVIE::R =0.0
    κVEI::R =0.0
    κVII::R =0.005
    VsynEE::R =15.0
    VsynIE::R =15.0
    VsynEI::R = -5.
    VsynII::R = -5.
    κ::R=0.5
 end