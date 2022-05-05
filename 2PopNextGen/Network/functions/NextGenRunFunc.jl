function NGModelRun(NGp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
  
    
    R = zeros(N,N,nWindows)
    Wsave = zeros(N,N,nWindows)
    nP = networkParameters(W, lags, N)
   
    for j = 1:nWindows
        if nWindows > 1
            println("working on window  : ",j)
        end

     
        if j == 1
            u0 = zeros(8N)
            u0[:] = makeInitConds(NGp)
            vP = variousPars(0.0, 100.0)
            
            hparams = u0
        else
            opts.stimOpt = "off"
            u0 = sol[:,end]
            iStart = findfirst(sol.t .> tWindows - 1.0)
            u_hist = make_uhist(sol.t[iStart:end] .- sol.t[end],sol[:,iStart:end])
            hparams = u_hist
          
            vP = variousPars(0.0, 0.01)
            
        end
        
        tspan = (0.0,tWindows)
        adpStops = collect(0.01:0.01:tWindows)
        #println(adpTime)
       

        p = NGp,nP,vP,stimNodes,Tstim,hparams,j,minSC,W_sum,opts
        if j == 1
            prob = DDEProblem(NextGen,u0, h1, tspan, p)
        else
            prob = DDEProblem(NextGen,u0, h2, tspan, p)
        end
        
     
        global sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,tstops=adpStops,saveat=0.01)


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