function WC(du,u,h,p,t)
    WCp,nP,vP,stimNodes,Tstim,hparams,nRun,minSC,W_sum,opts = p

    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ = WCp
    @unpack W,lags,N = nP
    @unpack tPrev,timeAdapt = vP
    @unpack stimOpt,adapt=opts

    if ((t==round(timeAdapt,digits=3)  && (t != tPrev)) || (t - tPrev > 0.01)) && (adapt == "on") && t >= 100.0

        nP.W .= adapt_global_coupling(hparams,N,W,lags,h,t,u,minSC,W_sum,vP)
        #vP.timeAdapt += 0.01
        #println("t = ", t, "adpTime = ",vP.timeAdapt)
        
    end
    @inbounds for i = 1:N
        d = 0.0
        #if t > 15.0
       
      

       
        @inbounds for j = 1:N
            if W[i,j] > 0.0
                if lags[i,j] > 0.0
                    d += W[i,j]*h(hparams,t-lags[i,j],idxs=j)
              
                else
                    d += W[i,j]*u[j]
                end
            end
        end
        #println(d)
        E = u[i]
        I = u[i+N]
        du[i] = (1/τE)*(-E + f(cEE*E + stim(t,i,stimNodes,Tstim,nRun,stimOpt)+  cIE*I + u[i+2N]+Pext + (η)*d,β,θE))
        du[i+N] =(1/τI)*( -I + f(cEI*E + cII*I+u[i+2N],β,θI) )
        du[i+2N] = (-1/τx)*u[i+2N]
    end
end

function dW(du,u,h,p,t)
    WCp,nP,vP,stimNodes,Tstim,hparams,nRum,minSC,W_sum,opts = p
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,η,σ = WCp
    @unpack W,lags,N = nP
    for i = 1:N
        du[i+2N] = σ
    end
end