function WCISP(du,u,h,p,t)
    WCp,nP,adpTime,stimNodes,Tstim,hparams,nRun,minSC,W_sum,opts = p

    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ,τISP,ρ = WCp
    @unpack W,lags,N = nP
    @unpack stimOpt,adapt=opts
  
    @inbounds for i = 1:N
        d = 0.0
        #if t > 15.0
        if adapt == "on" || adapt == "ON"
            if t ≈ adpTime
                nP.W .= adapt_global_coupling(hparams,N,W,lags,h,t,u,minSC,W_sum)
                adpTime += 0.01
            end
        end
      
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
        du[i] = (1/τE)*(-E + f(cEE*E + stim(t,i,stimNodes,Tstim,nRun,stimOpt) + u[i+2N]*I + u[i+3N]+Pext + (η)*d,β,θE))
        du[i+N] =(1/τI)*( -I + f(cEI*E + cII*I+u[i+3N],β,θI) )
        du[i+2N] = (1/τISP)*I*(E - ρ)
        du[i+3N] = (-1/τx)*u[i+3N]
    end
end

function dW(du,u,h,p,t)
    WCp,nP,adpTime,stimNodes,Tstim,hparams,nRum,minSC,W_sum,opts = p
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ,τISP,ρ = WCp
    @unpack W,lags,N = nP
    for i = 1:N
        du[i+3N] = σ
    end
end