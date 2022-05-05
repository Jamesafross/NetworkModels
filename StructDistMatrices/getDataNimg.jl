function getDataNimg(;normalise=0,delayDigits=2,logSC=1)
    HOMEDIR=homedir()
    InDATADIR="$HOMEDIR/NetworkModels/StructDistMatrices/Nimg_data"

    SCDataDir="$InDATADIR/Structural"
    FCDataDir="$InDATADIR/Functional"
    DistDataDir="$InDATADIR/Dists"

    dataNamesSC = readdir("$SCDataDir")
    dataNamesDist = readdir("$DistDataDir")
    dataNamesFC_AL = readdir("$FCDataDir/AL")
    dataNamesFC_DP = readdir("$FCDataDir/DP")
    dataNamesFC_FL = readdir("$FCDataDir/FL")
    dataNamesFC_VL = readdir("$FCDataDir/VL")
    
    numDataSC = size(dataNamesSC,1)
    numDataDist = size(dataNamesDist,1)
    numDataFC_AL = size(dataNamesFC_AL,1)
    numDataFC_DP = size(dataNamesFC_DP,1)
    numDataFC_FL = size(dataNamesFC_FL,1)
    numDataFC_VL = size(dataNamesFC_VL,1)

    
    

    SC_Array = zeros(184,184,numDataSC)
    Dist_Array = zeros(184,184,numDataDist)
    FC_Array_AL = zeros(184,184,numDataFC_AL)
    FC_Array_DP = zeros(184,184,numDataFC_DP)
    FC_Array_FL = zeros(184,184,numDataFC_FL)
    FC_Array_VL = zeros(184,184,numDataFC_VL)
    monkeys = ["AL","DP","FL","VL"]

    for i = 1:numDataDist
        monkey = monkeys[i]
        Dist_Array[:,:,i] = load("$DistDataDir/dist_$monkey.jld","dist_$monkey")
    end


    for i = 1:numDataFC_AL
        FC_Array_AL[:,:,i] = load("$FCDataDir/AL/AL_FC_$i.jld", "AL_FC_$i")
    end
    for i = 1:numDataFC_DP
        FC_Array_DP[:,:,i] = load("$FCDataDir/DP/DP_FC_$i.jld", "DP_FC_$i")
    end
    for i = 1:numDataFC_FL
        FC_Array_FL[:,:,i] = load("$FCDataDir/FL/FL_FC_$i.jld", "FL_FC_$i")
    end
    for i = 1:numDataFC_VL
        FC_Array_VL[:,:,i] = load("$FCDataDir/VL/VL_FC_$i.jld", "VL_FC_$i")
    end

    for i = 1:numDataSC
        monkey = monkeys[i]
        SC_Array = load("$SCDataDir/SL_$monkey.jld", "SL_$monkey")
    end

    #load data and make struct & dist matrices

   


    return SC_Array,FC_Array_AL,FC_Array_DP,FC_Array_FL,FC_Array_VL,Dist_Array
end
