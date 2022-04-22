function bmRHS(du,u,p,t)
    bP,E = p 
    @unpack E_0,κ,γ,τ,α,V_0,TE,v_0,ϵ,r_0,k_1,k_2,k_3 = bP
    for i = 1:N
        s = u[i]
        fin = maximum([u[N+i],0.0])
        v= maximum([u[2N+i],0.0])
        q = u[3N+i]

        du[i] = 2.0*(E[i](t)) - κ*s - γ*(fin-1)
        du[N+i] = s
        du[2N+i] = (1/τ)*(fin-v^(1/α))
        du[3N+i] = (1/τ)*((fin*E_NL(fin,E_0))/(E_0)  - q*v^(1/α-1))
    end
end

function runBalloon(u0,balloonParams,tspan,saveat)
    bP,E = balloonParams

    prob = ODEProblem(bmRHS,u0,tspan,balloonParams)
    solBM = solve(prob,saveat = saveat)
   

    v_save = solBM[2N+1:3N,:]
    v_save[v_save .< 0] .= 0.0
    out = 100*DeltS_NL.(v_save,solBM[3N+1:4N,:],bP.V_0,bP.k_1,bP.k_2,bP.k_3)

    return out,v_save
end

