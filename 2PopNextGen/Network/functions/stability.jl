function f(F, x, p)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

    rE  = x[1]
    rI  = x[2]
    vE  = x[3]
    vI  = x[4]
    gEE = κSEE * rE
    gEI = κSEI * rI
    gIE = κSIE * rE
    gII = κSII * rI

    #rE
    F[1]=(-gEE*rE - gEI*rE - κVEE*rE - κVEI*rE + 2*rE*vE + (ΔE/(τE*pi)))
    #rI
    F[2]=(-gII*rI - gIE*rI - κVIE*rI - κVII*rI + 2*rI*vI + (ΔI/(τI*pi)))
    #vE
    F[3]=(κVEI*(vI-vE) + gEE*(VsynEE-vE)+gEI*(VsynEI-vE)-(τE^2)*(pi^2)*(rE^2)+(vE^2)+η_0E)

    #vI
    F[4]=(κVIE*(vE-vI) + gIE*(VsynIE-vI)+gII*(VsynII-vI)-(τI^2)*(pi^2)*(rI^2)+(vI^2)+η_0I)

end

function SteadyState(p)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI= p
    f!(F,x) = f(F,x,p)
    X = rand(4)
    SS = nlsolve(f!, [ X[1];X[2];X[3];X[4]])
    conds = false
    while conds == false
            X = rand(4)
            SS = nlsolve(f!, [ X[1];X[2];X[3];X[4]])
            if SS.f_converged == true  &&  (SS.zero[1] > 0 && SS.zero[2] > 0)
                conds = true
            end
    end
    
    return SS.zero
end



function init_conds_SS(p)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI= p
    x = SteadyState(p)
    rE  = x[1]
    rI  = x[2]
    vE  = x[3]
    vI  = x[4]
    gEE = κSEE * rE
    gEI = κSEI * rI
    gIE = κSIE * rE
    gII = κSII * rI

return [rE,rI,vE,vI,gEE,gIE,gEI,gII]

end


function find_steady_states(p)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI= p
    ssMat = []
    for i = 1:1000
        #print(i)
        global flagSS = false
            steadyStateTest =  SteadyState(p)'
            while steadyStateTest[1] < 0 || steadyStateTest[2] < 0
                steadyStateTest =  SteadyState(p)'
            end
        if i == 1
                ssMat = cat(ssMat,steadyStateTest,dims=1)
        else
            for j = 1:size(ssMat,1)
                if sum(abs.((steadyStateTest' .- ssMat[j,:]))) < 10^-1
                    global flagSS = true
                end
            end
            if flagSS  == false
                if steadyStateTest[1] > 0
                    if  steadyStateTest[3] > 0
                        ssMat =  cat(ssMat,steadyStateTest,dims=1)
                    end
                end
            end
        end
    end
    return ssMat
end
