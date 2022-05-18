#################################
# Functions for use in the DDE  #
# solver                        #
#################################
function NextGen(du,u,h,p,t)
    NGp,nP,vP,aP,stimNodes,Tstim,hparams,j,minSC,W_sum,opts = p


    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    @unpack κSEEv,κSIEv,κSEIv,κSIIv,κSUM = κS
    @unpack tPrev,timeAdapt,count = vP
    @unpack stimOpt,adapt = opts
    @unpack W,dist,lags,N = nP
    @unpack tP,HIST = aP 
    
    
 
   
   

    
    if  t >= tP && adapt == "on"
        #println(t)
      #  aP.HIST = hcat(aP.HIST,u[1:N])
       # if size(aP.HIST,2) > 100
        #    aP.HIST = aP.HIST[:, 1:end .!= 1]
        #end
        κS.κSEEv[i],κS.κSIEv[i],κS.κSEIv[i],κS.κSIIv[i] = adapt_local_func(h,hparams,t,κS,rE,rI,i,N,0.0000005)
        aP.tP += 0.01  
        aP.tP = round(aP.tP,digits=2)
        vP.count += 1
    end


 
      @inbounds for i = 1:N
          d = 0.0

          @inbounds for j = 1:N
              if W[j,i] != 0.
                  if lags[j,i] == 0.
                     d  += W[j,i]*u[j]   
                  else
                      d += W[j,i]*h(hparams,t-lags[j,i]; idxs=j)
                  end
              end
  
          end
  
          rE=u[i]
          rI=u[i+N]
          vE=u[i+2N]
          vI=u[i+3N]
          gEE=u[i+4N]
          gIE=u[i+5N]
          gEI=u[i+6N]
          gII=u[i+7N]
          


         # if ((t==round(timeAdapt,digits=3)  && (t != tPrev)) || (t - tPrev > 0.01)) && (adapt == "on") && t >= timeAdapt
            

            #κS.κSEEv[i],κS.κSIEv[i],κS.κSEIv[i],κS.κSIIv[i] = adapt_local_func(h,hparams,t,κS,rE,rI,i,N,0.0000005)
            #println(cor(aP.HIST[1,:],aP.HIST[100,:]))
            #if i == N
               # nP.W = adapt_global_coupling(hparams,N,W,lags,h,t,u,minSC,W_sum)
               # vP.timeAdapt += 0.01
               # vP.tPrev = maximum([vP.tPrev,t])
                
            #end
           
            #println(t)
        #end

            
         
          #rE
          du[i] =(1. /τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2. * rE * vE + (ΔE / (τE*pi)))
          #rI
          du[i+N] =(1. /τI)*(-gIE*rI - gII*rI -κVIE*rI - κVII*rI + 2. * rI * vI + (ΔI / (τI*pi)))
          #vE
          du[i+2N] =(1. /τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vI) + κVEI*(vI - vE) - (τE^2)*(pi^2) * (rE^2.) +  vE^2. + η_0E +stim(t,i,stimNodes,Tstim,2,stimOpt))
          #vI
          du[i+3N] =(1. /τI)*(gIE*(VsynIE - vE) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2)*(pi^2)* (rI^2.) + vI^2. + η_0I)
          #gEE
          du[i+4N] = αEE * (-gEE + κ*d + κSEEv[i]*rE)
          #gIE
          du[i+5N] = αIE * (-gIE + κSIEv[i] * rE)
          #gEI
          du[i+6N] = αEI * (-gEI + κSEIv[i] * rI)
          #gII
          du[i+7N] = αII * (-gII + κSIIv[i] * rI)
  
  
      end
  
  end


  function NextGen_test_(du,u,h,p,t)
    NGp,nP,hparams= p


    ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
 
    W,lags,N = nP
  
 
      @inbounds for i = 1:N
          d = 0.0

          @inbounds for j = 1:N
              if W[i,j] != 0.
                  if lags[i,j] == 0.
                     d  += W[i,j]*u[j]   
                  else
                      d += W[i,j]*h(hparams,t-lags[i,j]; idxs=j)
                  end
              end
  
          end
  
          rE=u[i]
          rI=u[i+N]
          vE=u[i+2N]
          vI=u[i+3N]
          gEE=u[i+4N]
          gIE=u[i+5N]
          gEI=u[i+6N]
          gII=u[i+7N]
         
          #rE
          du[i] =(1. /τE)*(-rE*(gEE + gEI + κVEE +  κVEI - 2 * vE) + (ΔE / (τE*pi)))
          #rI
          du[i+N] =(1. /τI)*(-rI*(gIE + gII + κVIE + κVII - 2 * vI) + (ΔI / (τI*pi)))
          #vE
          du[i+2N] =(1. /τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vI) + κVEI*(vI - vE) - (τE^2)*(pi^2) * (rE^2.) +  vE^2. + η_0E)
          #vI
          du[i+3N] =(1. /τI)*(gIE*(VsynIE - vE) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2)*(pi^2)* (rI^2.) + vI^2. + η_0I)
          #gEE
          du[i+4N] = αEE * (-gEE + κSEE * κ*d + κSEE*rE)
          #gIE
          du[i+5N] = αIE * (-gIE + κSIE * rE)
          #gEI
          du[i+6N] = αEI * (-gEI + κSEI * rI)
          #gII
          du[i+7N] = αII * (-gII + κSII * rI)
  
  
      end
  
  end