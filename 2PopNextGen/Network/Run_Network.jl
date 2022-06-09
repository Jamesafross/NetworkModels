include("./setup.jl")

for jj = 1:length(Run_vec)

    κS.κSEEv = ones(N)*NGp.κSEE
    κS.κSIEv = ones(N)*NGp.κSIE
    κS.κSEIv = ones(N)*NGp.κSEI
    κS.κSIIv = ones(N)*NGp.κSII
    κS.κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]

    for setstim = ["on","off"]
        Run = string(Int(round(Run_vec[jj])))

    global opts.stimOpt = setstim

    println("Running model ... ")
    @time out,weightSaved = NGModelRun(u0,κS,wS)

    BOLD_OUT=[]
    for ii = 1:nWindows
        if ii == 1
            BOLD_OUT= out[:,:,ii]
        else
            BOLD_OUT = cat(BOLD_OUT,out[:,:,ii],dims=2)
        end
    end


    if stimOpt == "on"
        save1 = "stim"
    else
        save1="NOstim"
    end
    if adapt == "on"
        save2 = "Adaptivity"
    else
        save2="NOadaptivity"
    end
    

    savename = save1*save2
    dir0 = "LR_"
    dir1 = string(LR)
    dir2 = "_StimStr_"
    dir3 = string(stimStr)

    savedir = dir0*dir1*dir2*dir3

    

   

    if save_data =="true"
        save("$OutDATADIR/$savedir/BOLD_$savename.jld","BOLD_$savename",BOLD_OUT)
        save("$OutDATADIR/$savedir/weights_$savename.jld","weights_$savename",weightSaved)
    end


    end

end

