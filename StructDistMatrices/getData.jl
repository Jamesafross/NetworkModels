function getData(c;normalise=0,delayDigits=2)
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices/Averaged"
    #load data and make struct & dist matrices

    SC = load("$InDATADIR/PaulStructmean_140.jld","PaulStructmean_140")
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


function getData_nonaveraged(c;normalise=0,delayDigits=2)
    InMeanDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices/Averaged"
    #load data and make struct & dist matrices
    FCDataDir = "$HOMEDIR/NetworkModels/StructDistMatrices/non-averaged/paulFunctional"
    SCDataDir = "$HOMEDIR/NetworkModels/StructDistMatrices/non-averaged/paulStructural"

    dataNamesFC = readdir("$FCDataDir")
    numDataFC = size(dataNamesFC,1)

    dataNamesSC = readdir("$SCDataDir")
    numDataSC = size(dataNamesSC,1)

    SC_Array = zeros(140,140,numDataSC)
    FC_Array = zeros(140,140,numDataFC)
    for i = 1:numDataSC
        SC_Array[:,:,i] = load("$SCDataDir/paulStruct_140_$i.jld", "paulStruct_140_$i")
    end

    for i = 1:numDataFC
        FC_Array[:,:,i] = load("$FCDataDir/paulFC_140_$i.jld", "paulFC_140_$i")
    end
  
  
    dist = load("$InMeanDATADIR/PaulDist.jld","dist")
    N = size(SC,1) # number of nodes
    lags = round.(dist./c,digits=delayDigits) # axonal delays
    

    minSC_vec = zeros(numDataSC) 
    W_sum_mat = zeros(N,numDataSC)
    for j = 1:numDataSC
        SC = SC_Array[:,:,j][:,:]
        minSC[j] = minimum(SC[SC.>0.0])
        for i = 1:N
            W_sum[i,j] = sum(SC[:,i])
        end 
    end

    

    return SC_Array,minSC_vec,W_sum_mat,lags,FC_Array,N
end