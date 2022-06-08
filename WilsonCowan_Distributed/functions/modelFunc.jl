function WCRun(WCp,bP,nWindows,tWindows,W,lags,N,minSC,W_sum,opts)
    adpTime = 15.0
    Rvec = zeros(N,N,nWindows)
    Wvec = zeros(N,N,nWindows)
    nP = networkParameters(W, lags, N)
    
        
    for j = 1:nWindows
       
        println("working on window ",j)
        if j == 1
            u0 = 0.1rand(3N)
        
           
            hparams = u0
            vP = variousPars(0.0, 50.0)

        else
            opts.stimOpt = "off"
            u0 = sol[:,end]
            iStart = findfirst(sol.t .> tWindows - 1.0)
            u_hist = make_uhist(sol.t[iStart:end] .- sol.t[end],sol[:,iStart:end])
            hparams = u_hist
           
            vP = variousPars(0.0, 0.01)
        end
        
        tspan = (0.0,tWindows)
        #println(adpTime)
        println(nP.W[1,2])
        println(adpTime)

        p = WCp,nP,vP,stimNodes,Tstim,hparams,j,minSC,W_sum,opts
        if j == 1
            prob = SDDEProblem(WC, dW,u0, h1, tspan, p)
        else
            prob = SDDEProblem(WC, dW,u0, h2, tspan, p)
        end
        @time global sol = solve(prob,EM(),dt=0.001,maxiters = 1e20)
       
    
        E_interp = make_In(sol.t,sol[1:N,:])
        tspanB = (sol.t[1],sol.t[end])
        balloonParams = bP,E_interp
        b0 =  cat(zeros(N),ones(3N),dims=1)
        println("running balloon model...")
        out,v_save = runBalloon(b0,balloonParams,tspanB,collect(sol.t[1]+15:2:sol.t[end]))
        println("done")
        out_trans=(out')    

        for n=1:N
            for m=n+1:N
                Rvec[n,m,j]=(cor(out_trans[:,m],out_trans[:,n]))
                Rvec[m,n,j]=Rvec[n,m,j];
            end
        end

        Wvec[:,:,j] = nP.W

    end

    return Rvec,Wvec

end


