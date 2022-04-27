function getData(c;normalise=0,delayDigits=2)
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices"
    #load data and make struct & dist matrices

    SC = load("$InDATADIR/PaulSC.jld","C")
    dist = load("$InDATADIR/PaulDist.jld","dist")
    PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")
    N = size(SC,1) # number of nodes
    lags = round.(dist./c,digits=delayDigits) # axonal delays
    
    if normalise == 1
        SC = normalise(SC,N)
        W_sum = ones(N)
    else
        W_sum = zeros(N)
        for i = 1:N
            W_sum[i] = sum(SC[:,i])
        end 
    end
    minSC = minimum(SC[SC.>0.0])

    return SC,minSC,W_sum,lags,PaulFCmean,N
end