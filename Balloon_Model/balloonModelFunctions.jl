

function E_NL(fin::Float64, E_0::Float64)
    return 1.0 - (1.0 - E_0)^(1.0 / fin)
end

function  DeltS_NL(v,q,V_0::Float64,k_1::Float64,k_2::Float64,k_3::Float64)
    return  V_0 .* (k_1 .* (1.0 .- q) .+ k_2 .* (1.0 .- q ./ v) .+ k_3 .* (1.0 .- v))
end

function bold_convolve(saveat,u0,T,p1,p2)
    tspan = (0.0,T)
    prob = ODEProblem(hemo_func_NL,u0,tspan,p1)
    sol = solve(prob,saveat=saveat)
    ts = LinRange(0,T,size(sol,2))

    if minimum(sol[3, :]) < 0.0001
        vol_flag = 1
    else
        vol_flag = 0
    end

end

function make_In(tgrid,u)
    t = tgrid[1]:(tgrid[end]-tgrid[1])/(size(tgrid,1) -1):tgrid[end]
    interp = []
    for i = 1:size(u,1)
        if i == 1
            interp = [CubicSplineInterpolation(t,u[i,:])]
        else
            interp = cat(interp,[CubicSplineInterpolation(t,u[i,:])],dims=1) 
        end
    end

    return interp
end


function input_to_function(input,t,dt)
    return input[Int(round(t/dt))+1]
end
