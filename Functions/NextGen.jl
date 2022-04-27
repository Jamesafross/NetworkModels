function makeInitConds(NGp)
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    
    params = κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    rE0, rI0, vE0, vI0, gEE0, gEI0, gIE0, gII0 = init_conds_SS(params)
    perturb = 0.1*rand(8*N)

    u0 = zeros(8*N)
    u0[1:N] .= rE0
    u0[N+1:2N] .= rI0
    u0[2N+1:3N] .= vE0
    u0[3N+1:4N] .= vI0
    u0[4N+1:5N] .= gEE0
    u0[5N+1:6N] .= gIE0
    u0[6N+1:7N] .= gEI0
    u0[7N+1:8N] .= gII0

    u0 = u0 + perturb
end

function NextGenRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
    adpTime = 15.0
    R = zeros(N,N,nWindows)
    Wsave = zeros(N,N,nWindows)
    nP = networkParameters(W, lags, N)
    for j = 1:nWindows
        if j == 1
            u0 = zeros(8N)
            u0[:] = makeInitConds(NGp)
            
            hparams = u0
        else
            opts.stimOpt = "off"
            u0 = sol[:,end]
            iStart = findfirst(sol.t .> tWindows - 1.0)
            u_hist = make_uhist(sol.t[iStart:end] .- sol.t[end],sol[:,iStart:end])
            hparams = u_hist
            adpTime = 0.01
        end
        tspan = (0.0,tWindows)
        adpStops = collect(adpTime:0.01:tWindows)
        p = NGp,nP,adpTime,stimNodes,Tstim,hparams,j,minSC,W_sum,opts
        if j == 1
            prob = DDEProblem(NextGen,u0, h1, tspan, p)
        else
            prob = DDEProblem(NextGen,u0, h2, tspan, p)
        end
        global sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,tstops=adpStops,saveat=0.05)
        gEE = sol[4N+1:5N,:]
        BalloonIn= make_In(sol.t,gEE)
        tspanB = (sol.t[1],sol.t[end])
        balloonParams = bP,BalloonIn
        b0 =  cat(zeros(N),ones(3N),dims=1)
        out,v_save = runBalloon(b0,balloonParams,tspanB,collect(sol.t[1]+15:2:sol.t[end]))  
        out_trans=(out')    
        for n=1:N
            for m=n+1:N
                R[n,m,j]=(cor(out_trans[:,m],out_trans[:,n]))
                R[m,n,j]=R[n,m,j];
            end
        end
        Wsave[:,:,j] = nP.W
    end
    return R,Wsave
end

