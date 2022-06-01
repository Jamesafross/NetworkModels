function NGModelRun(NGp,bP,nP,κS,opts,u0)
   @unpack W, dist,lags,N,minSC,W_sum = nP
   @unpack stimOpt,stimWindow,stimNodes,Tstim,adapt,synapses,tWindows,nWindows = opts
    
    R = zeros(N,N,nWindows)
    Wsave = zeros(N,N,nWindows)
    BOLD_saveat = collect(0:1.6:tWindows)
    size_out = length(BOLD_saveat)
    BOLD_out = zeros(N,size_out,nWindows)
    
   
   
    for j = 1:nWindows
       
 
         
        if nWindows > 1
            println("working on window  : ",j)
        end

     
        if j == 1
            
            global vP = variousPars(0.0, 100.0,0)
            global aP = adaptParams(100.01,u0[1:N])
            hparams = u0
          
        else
            u0 = sol[:,end]
            iStart = findfirst(sol.t .> tWindows - 1.1)
            rE0 = sol[1:N,:]
            u_hist = make_uhist(sol.t[iStart:end] .- sol.t[end],sol[1:2N,iStart:end])
            hparams = u_hist
            vP = variousPars(0.0, 0.01,0)
            aP.tP = 0.01
            println(size(aP.HIST))
         
        end
        
        tspan = (0.0,tWindows)
        adpStops = collect(0.01:0.01:tWindows)
        clags = cat(unique(reshape(lags[lags.>0.0],length(lags[lags.>0.0]))),1.0,dims=1)
        println(clags)
    
        global p = (NGp,nP,vP,aP,κS,hparams,j,opts)

        if opts.synapses == "1stOrder"
            probDDE = NextGen
        elseif opts.synapses == "2ndOrder"
            probDDE = NextGen2ndOrderSynapses
        end

        if j == 1
            prob = DDEProblem(probDDE,u0, h1, tspan, p)
            global sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,tstops=adpStops,saveat=0.01)
        else
            prob = DDEProblem(probDDE,u0, h2, tspan, p)
            global sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,tstops=adpStops,saveat=0.01)
        end
        vE = sol[2N+1:3N,:]
        vI = sol[3N+1:4N,:]
        gEE = sol[4N+1:5N,:]
        gIE = sol[5N+1:6N,:]
        gEI = sol[6N+1:7N,:]
        gII = sol[7N+1:8N,:]

       
        currentE =  (gEE.*(NGp.VsynEE .- vE) .+ gEI.*(NGp.VsynEI .- vE));
        currentI =  (gIE.*(NGp.VsynIE .- vI) .+ gII.*(NGp.VsynII .- vI));
        Current = currentE .+ currentI

    
        BalloonIn= make_In(sol.t,Current)
        tspanB = (sol.t[1],sol.t[end])
        balloonParams = bP,BalloonIn,N
        if j == 1
            b0 =  cat(zeros(N),ones(3N),dims=1)
        else
            b0 = endBM
        end
       
         global out,endBM = runBalloon(b0,balloonParams,tspanB,BOLD_saveat,N)
        
        out_trans=(out')    

        for n=1:N
            for m=n+1:N
                R[n,m,j]=(cor(out_trans[:,m],out_trans[:,n]))
                R[m,n,j]=R[n,m,j];
            end
        end

        Wsave[:,:,j] = nP.W
        BOLD_out[:,:,j] = out

    end

    return R,Wsave,BOLD_out


end