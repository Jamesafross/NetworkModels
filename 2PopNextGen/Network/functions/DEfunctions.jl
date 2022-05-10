#################################
# Functions for use in the DDE  #
# solver                        #
#################################
function NextGen(du,u,h,p,t)
    NGp,nP,vP,aP,stimNodes,Tstim,hparams,j,minSC,W_sum,opts = p
    
    @unpack ΔE,ΔI,η_0E,η_0I,τE,τI,αEE,αIE,αEI,αII,κSEE,κSIE,κSEI,
    κSII,κVEE,κVIE,κVEI,κVII,VsynEE,VsynIE,VsynEI,VsynII,κ = NGp
    @unpack tPrev,timeAdapt = vP
    @unpack stimOpt,adapt = opts
    @unpack W,dist,lags,N = nP
    @unpack tP,HIST = aP 
    
 
   
    if ((t==round(timeAdapt,digits=3)  && (t != tPrev)) || (t - tPrev > 0.01)) && (adapt == "on") && t >= timeAdapt
            
            #println(cor(aP.HIST[1,:],aP.HIST[100,:]))
             nP.W = adapt_global_coupling_cor(N,W,dist,minSC,W_sum,HIST,1.,2.,0.00001)
            vP.timeAdapt += 0.01
            vP.tPrev = maximum([vP.tPrev,t])
            #println(t)
    end
     
    if  t >= tP
            
        aP.HIST = hcat(aP.HIST,u[1:N])
        if size(aP.HIST,2) > 200
            aP.HIST = aP.HIST[:, 1:end .!= 1]
        end
        aP.tP += 0.01
        
       
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
         
          #rE
          du[i] =(1. /τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2. * rE * vE + (ΔE / (τE*pi)))
          #rI
          du[i+N] =(1. /τI)*(-gIE*rI - gII*rI -κVIE*rI - κVII*rI + 2. * rI * vI + (ΔI / (τI*pi)))
          #vE
          du[i+2N] =(1. /τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vI) + κVEI*(vI - vE) - (τE^2)*(pi^2) * (rE^2.) +  vE^2. + η_0E +stim(t,i,stimNodes,Tstim,1,stimOpt))
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