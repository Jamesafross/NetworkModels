function getData(c;normalise=0,delayDigits=2)
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices/Averaged"
    #load data and make struct & dist matrices

    SC = 0.1*load("$InDATADIR/PaulStructmean_140.jld","PaulStructmean_140")
    dist = load("$InDATADIR/PaulDist.jld","dist")
    PaulFCmean = load("$InDATADIR/PaulFCmean_140.jld","paulFCmean_140")
    N = size(SC,1) # number of nodes
    lags = round.(dist./c,digits=delayDigits) # axonal delays
    SC[SC .<0.01] .=0.0
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


function getData_nonaveraged(;normalise=0,delayDigits=2,SCtype="log")
    InMeanDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices/Paul/Averaged"
    #load data and make struct & dist matrices
    FCDataDir = "$HOMEDIR/NetworkModels/StructDistMatrices/Paul/non-averaged/paulFunctional"
    SCDataDir = "$HOMEDIR/NetworkModels/StructDistMatrices/Paul/non-averaged/paulStructural"

    dataNamesFC = readdir("$FCDataDir")
    numDataFC = size(dataNamesFC,1)

    dataNamesSC = readdir("$SCDataDir")
    numDataSC = size(dataNamesSC,1)

    SC_Array = zeros(139,139,numDataSC)
    FC_Array = zeros(139,139,numDataFC)
    for i = 1:numDataSC
        SC = load("$SCDataDir/paulStruct_140_$i.jld", "paulStruct_140_$i")
        SC_Array[:,:,i] = SC[1:size(SC,1) .!= 100,1:size(SC,1) .!= 100 ]
    end 

    if SCtype == "log"
        SC_Array = log.(SC_Array)
        SC_Array[SC_Array .< 0.0] .= 0.0
    end


    for i = 1:numDataFC
        FC = load("$FCDataDir/paulFC_140_$i.jld", "paulFC_140_$i")
        FC_Array[:,:,i] =  FC[1:size(FC,1) .!= 100,1:size(FC,1) .!= 100 ]
    end
  
  
    dist = load("$InMeanDATADIR/PaulDist.jld","dist")
   

    return SC_Array,FC_Array,dist
end


function getFCstim_nonaverages()

    #load data and make struct & dist matrices
    FCDataDir = "$HOMEDIR/NetworkModels/StructDistMatrices/Paul/non-averaged/paulFunctional-stim"
 

    dataNamesFC = readdir("$FCDataDir")
    numDataFC = size(dataNamesFC,1)


    FC_Array = zeros(139,139,numDataFC)
  
    for i = 1:numDataFC
        FC = load("$FCDataDir/paulFC_stim_140_$i.jld", "paulFC_stim_140_$i")
    
        FC_Array[:,:,i] =  FC[1:size(FC,1) .!= 100,1:size(FC,1) .!= 100 ]
    end

    return FC_Array

end

function getMinSC_and_Wsum(SC)

    minSC = minimum(SC[SC .> 0.0])
    N = size(SC,1)
    W_sum = zeros(N)
    for i = 1:N
        W_sum[i] = sum(SC[:,i])
    end 
    return minSC, W_sum
end


