#=
FUNCTIONS FOR DATA SAVING, READING AND DATA ANALYSIS (TULLY-LIKE PLOTS FOR FINAL POPULATIONS)
=#

function oldfile(filename::String)
    if isfile(filename)
        oldname = filename*".old"
        if isfile(oldname)
            oldcount = 1
            while isfile(oldname*"$oldcount")
                oldcount += 1
            end
            oldcount -= 1
            for (i,oldnum) in enumerate(oldcount:-1:1)
                run(`mv -f $oldname$oldnum $oldname$(oldnum+1)`)
                println("Moved file $oldname$oldnum to $oldname$(oldnum+1)")
            end
            run(`mv -f $oldname $(oldname)1`)
            println("Moved file $oldname to $(oldname)1")
        end
        run(`mv -f $filename $filename.old`)
        println("Moved file $filename to $oldname")
    end
end

function folder_merge(folder,dyn)
    FILES = readdir(folder)
    F1 = FILES[1]
    for fName in FILES[2:end]
        merge_save(fName,F1,[dyn])
    end
    println("Finished merging $folder !")
end


function merge_save(oldFile,newFile,DYN_LIST=false,remove=true) #merge save only useful for save files of multiple trajectories!
    oldFid = h5open(oldFile,"cw")
    oldNames = names(oldFid)
    newFid = h5open(newFile,"cw")
    newNames = names(newFid)
    if DYN_LIST == false
        dListOld = read(oldFid,"DYN_LIST")
        dListNew = read(newFid,"DYN_LIST")
    else
        dListOld = DYN_LIST
        dListNew = DYN_LIST
    end
    
    for dyn in dListOld
        println("Starting $dyn transfer from $oldFile to $newFile...")
        if !(dyn in dListNew)
            println("$dyn not in $newFile, migrating to new file...")
            h5write(newFile,"R_"*dyn,read(oldFid,"R_"*dyn))
            h5write(newFile,"P_"*dyn,read(oldFid,"P_"*dyn))
            if !(dyn in CL_LIST)
                h5write(newFile,"C_REAL_"*dyn,read(oldFid,"C_REAL_"*dyn))
                h5write(newFile,"C_IMAG_"*dyn,read(oldFid,"C_IMAG_"*dyn))
                if dyn in SH_LIST
                    h5write(newFile,"AST"*dyn,read(oldFid,"AST"*dyn))
                end
            end
        else
            Rold = read(oldFid,"R_"*dyn)
            Rnew = read(newFid,"R_"*dyn)
            Pold = read(oldFid,"P_"*dyn)
            Pnew = read(newFid,"P_"*dyn)
            flags = length(Rold[:,1,1])
            NDOFs = length(Rold[1,:,1])
            NtrajsOld = length(Rold[1,1,:]) # R[flags,NDOFs,Ntrajs]
            NtrajsNew = length(Rnew[1,1,:])
            @show NtrajsOld, NtrajsNew
            Rmerged = zeros(flags,NDOFs,NtrajsOld + NtrajsNew)
            Rmerged[:,:,1:NtrajsNew] = Rnew
            Rmerged[:,:,NtrajsNew+1:end] = Rold
            Pmerged = zeros(flags,NDOFs,NtrajsOld + NtrajsNew)
            Pmerged[:,:,1:NtrajsNew] = Pnew
            Pmerged[:,:,NtrajsNew+1:end] = Pold
            o_delete(newFid,"R_"*dyn)
            o_delete(newFid,"P_"*dyn)
            write(newFid,"R_"*dyn,Rmerged)
            write(newFid,"P_"*dyn,Pmerged)
            println("Transfered R and P...")
            if !(dyn in CL_LIST)
                Cold = read(oldFid,"C_REAL_"*dyn) + 1im * read(oldFid,"C_IMAG_"*dyn)
                Cnew = read(newFid,"C_REAL_"*dyn) + 1im * read(newFid,"C_IMAG_"*dyn)
                nsts = length(Cold[1,:,1])
                Cmerged = zeros(Complex,flags,nsts,NtrajsOld + NtrajsNew) #C[flags,nsts,Ntrajs]
                Cmerged[:,:,1:NtrajsNew] = Cnew
                Cmerged[:,:,NtrajsNew+1:end] = Cold
                o_delete(newFid,"C_REAL_"*dyn)
                o_delete(newFid,"C_IMAG_"*dyn)
                write(newFid,"C_REAL_"*dyn,real.(Cmerged))
                write(newFid,"C_IMAG_"*dyn,imag.(Cmerged))
                println("Transfered C")
                if dyn in SH_LIST
                    ASTold = read(oldFid,"AST_"*dyn)
                    ASTnew = read(newFid,"AST_"*dyn)
                    ASTmerged = zeros(flags,NtrajsOld + NtrajsNew) #AST[flags,Ntrajs]
                    ASTmerged[:,1:NtrajsNew] = ASTnew
                    ASTmerged[:,NtrajsNew+1:end] = ASTold
                    o_delete(newFid,"AST_"*dyn)
                    write(newFid,"AST_"*dyn,ASTmerged)
                    println("Transfered AST")
                end
            end
            println("Finished for $dyn")
        end
    end
    close(oldFid)
    close(newFid)

    if remove
        rm(oldFile)
        println("Erasing $oldFile")
    end
end


function migrate_save(oldFile,newFile,remove = true)
    oldFid = h5open(oldFile,"cw")
    oldNames = names(oldFid)
    newFid = h5open(newFile,"cw")
    newNames = names(newFid)

    for name in oldNames
        if name in newNames
            println("Warning, $name appears in both $oldFile and $newFile")
        else
            h5write(newFile,name,h5read(oldFile,name))
            println("Transfered $name from $oldFile to $newFile")
        end
    end

    println("Finished transfer")
    close(newFid)
    close(oldFid)
    if remove
        rm(oldFile)
        println("Erasing $oldFile")
    end
end

function CL_save(Tf,Rvec,pvec,file,prefix)
    while isfile(file*"OPEN") #make sure no other subprocess is accessing the file
        sleep(0.1)
    end
    touch(file*"OPEN")
    fid = h5open(file,"cw")
    namelist = names(fid)
    if !("T_$prefix" in namelist) #Time save does not exist for prefix, data creation routine
        h5write(file,"T_"*prefix,Tf)
        N = d_create(fid, "Ntrajs_"*prefix, Int64, ((1,),(1,)), "chunk", (1))
        N[1] = 1
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        R = d_create(fid, "R_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        R[:,:,1] = Rvec
        P = d_create(fid, "P_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        P[:,:,1] = pvec
        println("Created R_$prefix and P_$prefix CL save")
    else
        N = d_open(fid,"Ntrajs_"*prefix)
        N[1] += 1
        Told = h5read(file,"T_"*prefix)
        if Told != Tf
            error("Warning, time arrays in savefile and simulation are not equal, FIRST_RUN should be set true to start new save")
        end
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        R = d_open(fid,"R_"*prefix)
        P = d_open(fid,"P_"*prefix)
        set_dims!(R,(flags,NDOFs,N[1]))
        set_dims!(P,(flags,NDOFs,N[1]))
        R[:,:,N[1]] = Rvec
        P[:,:,N[1]] = pvec
        println("Saved new trajectory on top of $(N[1] - 1) already existing ones")
    end
    close(fid)
    rm(file*"OPEN")
end


function MF_save(Tf,Rvec,pvec,C,file,prefix)
    while isfile(file*"OPEN") #make sure no other subprocess is accessing the file
        sleep(0.1)
    end
    touch(file*"OPEN")
    fid = h5open(file,"cw")
    namelist = names(fid)
    if !("T_$prefix" in namelist) #Time save does not exist for prefix, data creation routine
        h5write(file,"T_"*prefix,Tf)
        N = d_create(fid, "Ntrajs_"*prefix, Int64, ((1,),(1,)), "chunk", (1))
        N[1] = 1
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        trSts = length(C[1,:])
        R = d_create(fid, "R_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        R[:,:,1] = Rvec
        P = d_create(fid, "P_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        P[:,:,1] = pvec
        Creal = d_create(fid, "C_REAL_"*prefix, Float64, ((flags,trSts,1),(flags,trSts,-1)), "chunk", (flags,trSts,1))
        Creal[:,:,1] = real.(C)
        Cimag = d_create(fid, "C_IMAG_"*prefix, Float64, ((flags,trSts,1),(flags,trSts,-1)), "chunk", (flags,trSts,1))
        Cimag[:,:,1] = imag.(C)
        println("Created R_$prefix, P_$prefix and C_$prefix MF save")
    else
        N = d_open(fid,"Ntrajs_"*prefix)
        N[1] += 1
        NtrajsOld = N[1] - 1
        Told = h5read(file,"T_"*prefix)
        if Told != Tf
            error("Warning, time arrays in savefile and simulation are not equal, FIRST_RUN should be set true to start new save")
        end
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        trSts = length(C[1,:])
        R = d_open(fid,"R_"*prefix)
        P = d_open(fid,"P_"*prefix)
        Creal = d_open(fid,"C_REAL_"*prefix)
        Cimag = d_open(fid,"C_IMAG_"*prefix)
        set_dims!(R,(flags,NDOFs,NtrajsOld + 1))
        set_dims!(P,(flags,NDOFs,NtrajsOld + 1))
        set_dims!(Creal,(flags,trSts,NtrajsOld + 1))
        set_dims!(Cimag,(flags,trSts,NtrajsOld + 1))
        R[:,:,NtrajsOld+1] = Rvec
        P[:,:,NtrajsOld+1] = pvec
        Creal[:,:,NtrajsOld+1] = real.(C)
        Cimag[:,:,NtrajsOld+1] = imag.(C)
        println("Saved new trajectory on top of $NtrajsOld already existing ones")
    end
    close(fid)
    rm(file*"OPEN")
end

function SH_save(Tf,Rvec,pvec,C,AST,file,prefix)
    while isfile(file*"OPEN") #make sure no other subprocess is accessing the file
        sleep(0.1)
    end
    touch(file*"OPEN")
    fid = h5open(file,"cw")
    namelist = names(fid)
    if !("T_$prefix" in namelist) #Time save does not exist for prefix, data creation routine
        h5write(file,"T_"*prefix,Tf)
        N = d_create(fid, "Ntrajs_"*prefix, Int64, ((1,),(1,)), "chunk", (1))
        N[1] = 1
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        trSts = length(C[1,:])
        R = d_create(fid, "R_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        R[:,:,1] = Rvec
        P = d_create(fid, "P_"*prefix, Float64, ((flags,NDOFs,1),(flags,NDOFs,-1)), "chunk", (flags,NDOFs,1))
        P[:,:,1] = pvec
        Creal = d_create(fid, "C_REAL_"*prefix, Float64, ((flags,trSts,1),(flags,trSts,-1)), "chunk", (flags,trSts,1))
        Creal[:,:,1] = real.(C)
        Cimag = d_create(fid, "C_IMAG_"*prefix, Float64, ((flags,trSts,1),(flags,trSts,-1)), "chunk", (flags,trSts,1))
        Cimag[:,:,1] = imag.(C)
        Ast = d_create(fid, "AST_"*prefix, Float64, ((flags,1),(flags,-1)), "chunk", (flags,1))
        Ast[:,1] = AST 
        println("Created R_$prefix, P_$prefix, C_$prefix and AST_$prefix SH save")
    else
        N = d_open(fid,"Ntrajs_"*prefix)
        N[1] += 1
        NtrajsOld = N[1] - 1
        Told = h5read(file,"T_"*prefix)
        if Told != Tf
            error("Warning, time arrays in savefile and simulation are not equal, FIRST_RUN should be set true to start new save")
        end
        flags = length(Rvec[:,1])
        NDOFs = length(Rvec[1,:])
        trSts = length(C[1,:])
        R = d_open(fid,"R_"*prefix)
        P = d_open(fid,"P_"*prefix)
        Creal = d_open(fid,"C_REAL_"*prefix)
        Cimag = d_open(fid,"C_IMAG_"*prefix)
        Ast = d_open(fid,"AST_"*prefix)
        NtrajsOld = size(R)[3]
        set_dims!(R,(flags,NDOFs,NtrajsOld + 1))
        set_dims!(P,(flags,NDOFs,NtrajsOld + 1))
        set_dims!(Creal,(flags,trSts,NtrajsOld + 1))
        set_dims!(Cimag,(flags,trSts,NtrajsOld + 1))
        set_dims!(Ast,(flags,NtrajsOld + 1))
        R[:,:,NtrajsOld+1] = Rvec
        P[:,:,NtrajsOld+1] = pvec
        Creal[:,:,NtrajsOld+1] = real.(C)
        Cimag[:,:,NtrajsOld+1] = imag.(C)
        Ast[:,NtrajsOld + 1] = AST
        println("Saved new trajectory on top of $NtrajsOld already existing ones")
    end
    close(fid)
    rm(file*"OPEN")
end

# = Full save format, check for K_SIMULATIONS!
function CL_FULL_SAVE(Tf,Rvec,pvec,file,prefix)
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
end

function MF_FULL_SAVE(Tf,Rvec,pvec,C,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
end

function SH_FULL_SAVE(Tf,Rvec,pvec,C,AST,file,prefix)
    Creal=Float64.(real.(C))
    Cim=Float64.(imag.(C))
    h5write(file,"T_"*prefix,Tf)
    h5write(file,"R_"*prefix,Rvec)
    h5write(file,"P_"*prefix,pvec)
    h5write(file,"C_REAL_"*prefix,Creal)
    h5write(file,"C_IMAG_"*prefix,Cim)
    h5write(file,"AST_"*prefix,AST)
end
# =#


function CL_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    return T,R,P
end

function MF_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C
end

function SH_read(file,prefix)
    T=h5read(file,"T_"*prefix)
    R=h5read(file,"R_"*prefix)
    P=h5read(file,"P_"*prefix)
    C_REAL=h5read(file,"C_REAL_"*prefix)
    C_IMAG=h5read(file,"C_IMAG_"*prefix)
    AST=h5read(file,"AST_"*prefix)
    C=C_REAL+1im.*C_IMAG
    return T,R,P,C,AST
end

function CL_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FP_"*prefix,A[2])
end

function MF_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
end

function SH_K_save(A,file,prefix)
    h5write(file,"FR_"*prefix,A[1])
    h5write(file,"FM_"*prefix,A[2])
    h5write(file,"Fpop_"*prefix,A[3])
    h5write(file,"Fast_"*prefix,A[4])
end

function general_K_save(A,file,prefix)
    if length(A)==2
        CL_K_save(A,file,prefix)
    elseif length(A)==3
        MF_K_save(A,file,prefix)
    elseif length(A)==4
        SH_K_save(A,file,prefix)
    else
        error("Check our input, it should have 2, 3 or 4 entries!")
    end
end

function general_K_read(file,prefix)
    A=Array[]
    if prefix in METHOD_LIST
        FR=h5read(file,"FR_"*prefix)
        FM=h5read(file,"FM_"*prefix)
        push!(A,FR)
        push!(A,FM)
        if prefix in MF_LIST || prefix in SH_LIST
            Fpop=h5read(file,"Fpop_"*prefix)
            push!(A,Fpop)
            if prefix in SH_LIST
                Fast=h5read(file,"Fast_"*prefix)
                push!(A,Fast)
            end
        end
    end

    return A
end

function full_memory_save(T,S_ARRAY,file)
    prefix=S_ARRAY[1].prefix
    steps=length(S_ARRAY)
    NDOFs=S_ARRAY[1].cl.NDOFs
    if NDOFs==1
        R=zeros(steps)
        p=zeros(steps)
    else
        R=zeros(steps,NDOFs)
        p=zeros(steps,NDOFs)
    end

    if !(prefix in CL_LIST) #not classical, need electronic coefficients tracking
        eff_sts=length(S_ARRAY[1].el.C)
        C=zeros(Complex,steps,eff_sts)
        U=zeros(steps,nsts,nsts)
    end

    for k in 1:steps
        S=S_ARRAY[k]
        if NDOFs==1
            R[k]=S.cl.R
            p[k]=S.cl.p
        else
            R[k,:].=S.cl.R
            p[k,:].=S.cl.p
        end
        if !(prefix in CL_LIST) #not classical, need electronic coefficients tracking
            C[k,:].=S.el.C
            U[k,:,:].=S.el.Ua
        end
    end
    Creal=zeros(size(C))
    Cimag=zeros(size(C))
    Creal.=real.(C)
    Cimag.=imag.(C)

    h5write(file,"R_"*prefix,R)
    h5write(file,"p_"*prefix,p)
    if !(prefix in CL_LIST)
        h5write(file,"Creal_"*prefix,Creal)
        h5write(file,"Cimag_"*prefix,Cimag)
        h5write(file,"U_"*prefix,U)
    end
    if prefix in SH_LIST
        AST=[S_ARRAY[k].ast for k in 1:steps]
        h5write(file,"AST_"*prefix,AST)
    end
    h5write(file,"T_"*prefix,T)
end

function full_memory_load(prefix,file,flags=100)
    R=h5read(file,"R_"*prefix)
    p=h5read(file,"p_"*prefix)
    T=h5read(file,"T_"*prefix)
    steps=length(T)
    NDOFs=length(R[1,:])

    kmap=Int.(round.(collect(range(1,stop=steps,length=flags))))

    if !(prefix in CL_LIST)
        C=h5read(file,"Creal_"*prefix)+1im.*h5read(file,"Cimag_"*prefix)
        U=h5read(file,"U_"*prefix)
        if prefix in SH_LIST
            AST=h5read(file,"AST_"*prefix)
        end
    end

    if NDOFs==1
        if prefix in CL_LIST
            S_ARRAY=[builder_CL_state(R[k],p[k],prefix,NDOFs) for k in kmap]
        elseif prefix in MF_LIST
            S_ARRAY=[builder_MF_state(R[k],p[k],C[k,:],prefix,U[k,:,:],NDOFs) for k in kmap]
        else
            S_ARRAY=[builder_SH_state(R[k],p[k],C[k,:],AST[k],prefix,U[k,:,:],NDOFs) for k in kmap]
        end
    else
        if prefix in CL_LIST
            S_ARRAY=[builder_CL_state(R[k,:],p[k,:],prefix,NDOFs) for k in kmap]
        elseif prefix in MF_LIST
            S_ARRAY=[builder_MF_state(R[k,:],p[k,:],C[k,:],prefix,U[k,:,:],NDOFs) for k in kmap]
        else
            S_ARRAY=[builder_SH_state(R[k,:],p[k,:],C[k,:],AST[k],prefix,U[k,:,:],NDOFs) for k in kmap]
        end
    end

    return T[kmap],S_ARRAY
end
