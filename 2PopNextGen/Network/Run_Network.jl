include("./setup.jl")
LinearAlgebra.BLAS.set_num_threads(4)

Run_vec = LinRange(1,5,5)

for jj = 1:length(Run_vec)
   
   

    for setstim = ["on","off"]
        Run = string(Int(round(Run_vec[jj])))
        nP.W .= SC
        κS.κSEEv = ones(N)*NGp.κSEE
        κS.κSIEv = ones(N)*NGp.κSIE
        κS.κSEIv = ones(N)*NGp.κSEI
        κS.κSIIv = ones(N)*NGp.κSII
        κS.κSUM = κSEEv[1]+κSIEv[1]+κSEIv[1]+κSIIv[1]
        IC.u0 = init0
        
        
        wS.κSEEv[:,1] = κSEEv
        wS.κSIEv[:,1] = κSIEv
        wS.κSEIv[:,1] = κSEIv
        wS.κSIIv[:,1] = κSIIv
        wS.count = 2
        
        global opts.stimOpt = setstim

        println("Running model ... ")
        @time out,weightSaved = NGModelRun(κS,wS,start_adapt)

        BOLD_OUT=[]
        for ii = 1:nWindows
                if ii == 1
                    BOLD_OUT= out[:,:,ii]
                else
                    BOLD_OUT = cat(BOLD_OUT,out[:,:,ii],dims=2)
                end
        end


        if opts.stimOpt == "on"
            save1 = "stim"
        else
            save1="NOstim"
        end
    
        
        run = "_moredelays_Run_$Run_vec[jj]"
        savename = save1*run
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

