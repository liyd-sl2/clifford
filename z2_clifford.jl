#jun_09, Z2 model, large system, magnetization (distribution and time dependence)
#random t_s (quenched) random noise
function main_D_1(λ, μ, Λ, τ, T_CHAR, Δ, Δh, Δϕ, λX, RUN, dir)

    #global α = Α # project onto |1> + α |0> + |-1>
    #global β = 1 # initial state |1> + β |0> + |-1>
    #global p0 = π0 # prob of bond forming
    #global p1 = π1 # prob of bond breaking
    #global st = [0 0 1 0 α 0 1 0 0]/√(2+α^2)
    #global M0 = √(p0) * reshape(kron(conj(st), st), 9, 9)
    #global M1 = √(1 - p0) * reshape(kron(conj(st), st), 9, 9) + (diagm([1,1,1,1,1,1,1,1,1]) - reshape(kron(conj(st), st), 9, 9))
    #global P0 = p0 * reshape(kron(conj(st), st), 9, 9) # P0 = M0†M0
    #global P1 = diagm([1,1,1,1,1,1,1,1,1]) - P0  # P1 = M1†M1
    #global P0 = diagm([1,1,1,0,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0])
    #global P1 = diagm([0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1,0,1,1,1,0,1])
    #global J = 1.
    #global P0 = diagm([1-J,1+J,1+J,1-J])/√(2*(1+J^2))
    #global P1 = diagm([1+J,1-J,1-J,1+J])/√(2*(1+J^2))
    #global U
    global L = λ    # size of lattice
    global M = μ    # number of qutrits
    global LOOP = Λ # runs of averages
    global T = τ    # time slices in each run
    #global J = ζ    # strength of interaction
    #global p0 = .5
    global B = 2    # dim. of qudit space, B(its)

    global counter

    #global g = Array{Int64}(M, M)
    #global mark = Array{Bool}(M)
    #global visited = Array{Bool}(M)
    #global res = Array{Int64}(T, M)
    #global acc = Array{Float64}(T, 6)
    #global Qacc = Array{Float64}(T, 2)
    #global dimer = Array{Float64}(T, 3)
    #global C = Array{Float64}(T)
    global Sacc = Array{Float64}(T, div(M, 2))
    global Macc = Array{Float64}(T)
    #global Macc = Array{Float64}(T*div(M, 2)*3, 8)
    #global Rini = Array{Float64}(LOOP, 2)
    #global Rfin = Array{Float64}(LOOP, 2)
    #global Mini = Array{Float64}(LOOP)
    #global Mfin = Array{Float64}(LOOP)
    #global pos, wfn
    global wfn = Array{Complex{Float64}}(2, B^M)
    #global sch = Array{Complex{Float64}}(B^M, B^M)
    #global visited = Array{Int64}(2*M)

    #global MA = MB = div(M, 2)
    #global ρ = Array{Complex{Float64}}(2, B^M, B^M)
    #global σ = Array{Complex{Float64}}(2, B^M, B^M)

    function single_trit_amp(x)  #
        return 1
    end

    function normalize(v)
        return v/norm(v)
    end

    function HaarUMatrix(n)
        z = (randn(n,n) + im*randn(n,n) )/√2
        q,r = qr(z)
        d = Diagonal(r)
        ph = d/abs.(d)
        q = q*ph
        return q
    end

#function main()
    I = diagm([1, 1])
    tchar = T_CHAR # random real number, characteristic time scale
    #Θchar = 0 # random real angle, characteristic interaction
    #δh = [0 cos(Θchar)-sin(Θchar)*im; cos(Θchar)+sin(Θchar)*im 0]
    X = [0 1; 1 0]
    Y = [0 -im; +im 0]
    Z = [1 0; 0 -1]
    #δhs = [0 1; 1 0]
    #δhf = [0 -im; +im 0]
    #δhf = [0 1; 1 0]
    #δhf = [1 0; 0 -1]
    #Us = expm(im*(kron(X, I) + kron(I, X))*tchar)
    #Us1 = kron(X, I) / √2
    #Us2 = kron(I, X) / √2
    #Us = HaarUMatrix(4)
    #Uf = expm(im*(kron(δhf, I) + kron(I, δhf))*tchar*λX)
    #Uf = kron(I, I)

    #λMp = [1 +λX; +λX 1]/√(2(1+λX^2))
    #λMm = [1 -λX; -λX 1]/√(2(1+λX^2))

    #pp = λX^2/(1+λX^2)

    #λMp = √λX X
    #λMm = √(1-λX) I

        #(pos, wfn) = init(L, M)
        #if abs(Δϕ) < 0.0001
        #    factor = 1.
        #else
        #    factor = sin(Δϕ/2)/(Δϕ/2)
        #end

        Θ = (π/4 + Δh) + 2*(rand()-0.5) * Δ

        for time = 1:T
            for MA=1:div(M,2)
                Sacc[time, MA] = 0
            end
            #for i = 1:2*M*M
            #    Macc[time, i] = 0
            #end
            Macc[time] = 0
        end


#for loop = 1:LOOP

        now = 1
        nex = 2

        psi1 = [cos(Θ) sin(Θ)]
        psi2 = [0 0]
        tmp1 = 1
        tmp2 = 1
        for i = 1:M
            tmp1 = kron(tmp1, psi1)
            tmp2 = kron(tmp2, psi2)
        end
        for k = 1:B^M
            wfn[now, k] = tmp1[k]
            wfn[nex, k] = tmp2[k]
        end

        #println(norm(wfn[now, :]))
        #println(norm(wfn[nex, :]))

        depth = 0

    ts = Array{Float64}(M)
    for i = 1:M
        ts[i] = rand()*2π
    end
    for time = 1:T
            #print("time = ", time, "\n")
            for i = 1:M
                #if mod(i-time, 2) == 0 && (mod(M, 2) == 0 || i < M)
                if i < M && mod(i-time, 2) == 0 #open boundary condition
                    if i < M
                        j = i+1
                    else
                        j = 1
                    end
                    n = M

                    #Us = HaarUMatrix(4)
                    #Uf = HaarUMatrix(4)

                    #ts = rand()*2π
                    #ts = tchar
                    Us = expm(im*(kron(X, I) + kron(I, X))*ts[i])
                    Uf = kron(I, I)

                    for k = 1:B^M
                        wfn[nex,k] = 0
                    end

                    probf = 0
                    probs = 0
                    for k = 1:B^M
                        ikbit = ((k-1) & (1<<(n-i))) >> (n-i)
                        jkbit = ((k-1) & (1<<(n-j))) >> (n-j)

                        if (ikbit == jkbit)
                            probf += abs(wfn[now, k])^2
                        else
                            probs += abs(wfn[now, k])^2
                        end
                    end

                    #println("probf: ", probf)

                    if rand() < probf
                        for k = 1:B^M
                            ikbit = ((k-1) & (1<<(n-i))) >> (n-i)
                            jkbit = ((k-1) & (1<<(n-j))) >> (n-j)

                            if (ikbit == jkbit)
                                kbit = ikbit*B+jkbit+1
                                kk = k - ((k-1) & (1<<(n-i))) - ((k-1) & (1<<(n-j)))
                                for ik = 0:B-1
                                    for jk = 0:B-1
                                        wfn[nex, kk + ik<<(n-i) + jk<<(n-j)] += Uf[ik*B+jk+1, kbit] * wfn[now,k]
                                    end
                                end
                            end
                        end
                    else
                        for k = 1:B^M
                            ikbit = ((k-1) & (1<<(n-i))) >> (n-i)
                            jkbit = ((k-1) & (1<<(n-j))) >> (n-j)

                            if (ikbit != jkbit)
                                kbit = ikbit*B+jkbit+1
                                kk = k - ((k-1) & (1<<(n-i))) - ((k-1) & (1<<(n-j)))
                                for ik = 0:B-1
                                    for jk = 0:B-1
                                        wfn[nex, kk + ik<<(n-i) + jk<<(n-j)] += Us[ik*B+jk+1, kbit] * wfn[now,k]
                                    end
                                end
                            end
                        end
                    end

                    #println(norm(wfn[nex, :]))
                    wfn[nex, :] = normalize(wfn[nex, :])
                    #println("wfn: ", wfn[nex, :])
                    now = 3-now
                    nex = 3-nex
                end
            end

            #-----random-noise-------------------------------------------------
            if (mod(time, 2) == 0)
                sflip = Array{Int64}(M)
                for i = 1:M
                    if rand() < λX
                        sflip[i] = 1
                    else
                        sflip[i] = 0
                    end
                end
                for k = 1:B^M
                    wfn[nex,k] = 0
                end
                for k = 1:B^M
                    kbit = Array{Int64}(M)
                    ktmp = k-1
                    for i = 1:M
                        kbit[M+1-i] = mod(ktmp, B)
                        ktmp = div(ktmp, B)
                    end
                    ltmp = 0
                    for i = 1:M
                        tmp = mod(kbit[i]+sflip[i], 2)
                        ltmp = ltmp*2+tmp
                    end
                    l = ltmp+1
                    wfn[nex,l] = wfn[now, k]
                end
                now = 3-now
                nex = 3-nex
            end

            #-----measurement-of-observable-------------------------------------------------
            #if (mod(time, 2) == 0)
                depth += 1

                                #Sρtmp = 0
                                #Sσtmp = 0
                                #for k = 1:B^M
                                #    for l = 1:B^M
                                #        Sρtmp += abs(ρ[now, k, l])^2
                                #        Sσtmp += abs(σ[now, k, l])^2
                                #    end
                                #end
                                #Sacc[depth, 1] = Sρtmp
                                #Sacc[depth, 2] = Sσtmp
                                #diff = (ρ[now, :, :]-σ[now, :, :])
                                #diff = diff*diff
                                #Sacc[depth, 3] = abs(trace(sqrtm(diff)))
                                #Sacc[depth, 4] = abs(trace(diff))

                                #Macc[depth, 1] = Macc[depth, 2] = Macc[depth, 3] = Macc[depth, 4] = 0
                                #Macc[depth, 5] = Macc[depth, 6] = Macc[depth, 7] = Macc[depth, 8] = 0

                                #ff = abs(trace(ρ[now,:,:]))
                                for k = 1:B^M
                                    mtmp = 0
                                    for pos = 1:M
                                        btmp = ((k-1) & (1<<(M-pos))) >> (M-pos)
                                        mtmp += (2.0*btmp-1.0)
                                    end
                                    mtmp = mtmp / M
                                #    Macc[depth, 1] += mtmp   * abs(ρ[now, k, k]) / ff
                                    Macc[depth] += mtmp^2 * abs(wfn[now, k])^2
                                #    Macc[depth, 3] += mtmp^3 * abs(ρ[now, k, k]) / ff
                                #    Macc[depth, 4] += mtmp^4 * abs(ρ[now, k, k]) / ff
                                    #Macc[depth, 5] += mtmp   * abs(σ[now, k, k]) / ff
                                    #Macc[depth, 6] += mtmp^2 * abs(σ[now, k, k]) / ff
                                    #Macc[depth, 7] += mtmp^3 * abs(σ[now, k, k]) / ff
                                    #Macc[depth, 8] += mtmp^4 * abs(σ[now, k, k]) / ff
                                end
                #for MA=1:div(M,2)
                #    MB = M-MA
                #    S2 = 0
                #    sch = Array{Complex{Float64}}(B^MA, B^MB)
                #    for k1 = 1:B^MA
                #        for k2 = k1:B^MA
                #            tmp = 0
                #            for l = 1:B^MB
                #                tmp += wfn[now, (k1-1)<<MB+l] * conj(wfn[now, (k2-1)<<MB+l])
                #            end
                #            if k1 == k2
                #                S2 += real(tmp^2)
                #            else
                #                S2 += 2*real(tmp^2)
                #            end
                #        end
                #    end
                #    Sacc[time, MA] += -log(2,S2)/LOOP
                #end
                #for MA=1:div(M,2)
                #    MB = M-MA
                #    S2 = 0
                #    sch = Array{Complex{Float64}}(B^MA, B^MB)
                #    for k = 1:B^MA
                #        for l = 1:B^MB
                #            sch[k, l] = wfn[now, (k-1)<<MB+l]
                #        end
                #    end
                #    rho = sch * ctranspose(sch)
                #    S2 = real(trace(rho*rho))
                #    Sacc[depth, MA] += -log(2,S2)/LOOP
                #end

                #for k = 1:B^M
                #    mtmp = 0
                #    for pos = 1:M
                #        btmp = ((k-1) & (1<<(M-pos))) >> (M-pos)
                #        mtmp += (2.0*btmp-1.0)
                #    end
                #    mtmp = mtmp / M
                #    Macc[depth] += mtmp * abs(wfn[now, k])^2 / LOOP
                #end
                #for MA=1:div(M,2)
                #    MB = M-MA
                #    S2 = 0
                #    sch = Array{Complex{Float64}}(B^MA, B^MB)
                #    for k = 1:B^MA
                #        for l = 1:B^MB
                #            sch[k, l] = wfn[now, (k-1)<<MB+l]
                #        end
                #    end
                #    rho = sch * ctranspose(sch)
                #    S2 = real(trace(rho*rho))
                #    Sacc[depth, MA] += -log(2,S2)/LOOP
                #end

                #for k = 1:B^M
                #    for i = 1:M
                #        for j = 1:M
                #            ii = mod(i, M) + 1
                #            jj = mod(j, M) + 1
                #            ibit  = ((k-1) & (1<<(M-i))) >> (M-i)
                #            jbit  = ((k-1) & (1<<(M-j))) >> (M-j)
                #            iibit = ((k-1) & (1<<(M-ii))) >> (M-ii)
                #            jjbit = ((k-1) & (1<<(M-jj))) >> (M-jj)
                #            Macc[depth, (i-1)*M+j] += (2*ibit-1) * (2*jbit-1) * abs(wfn[now, k])^2
                #            Macc[depth, (i-1)*M+j+M*M] += (2*ibit-1) * (2*iibit-1) * (2*jbit-1) * (2*jjbit-1) * abs(wfn[now, k])^2
                #        end
                #    end
                #end
            #end
            #-------------------------------------------------------------------------------
    end
    #println("run = ", λX)
    #for i = 1:depth
    #    println(Δh, " ", i, " ", Macc[i])
    #end
#end
    #print(Macc)
    #writedlm(string(dir, "_N_",  M, "_Us_XX_Uf_I_Sacc", "_run_",  @sprintf("%3d", λX), ".csv"), Sacc, " ")
    writedlm(string(dir, "_N_",  M, "_Us_XX_Uf_I_Macc", "_h_", @sprintf("%.3f", λX), "_run_",  @sprintf("%3d", RUN), ".csv"), Macc, " ")
    #writedlm(string(dir, "_N_",  M, "_deltah_", @sprintf("%.3f",Δh), "_run_",  @sprintf("%3d", λX), "_Us_XX_Uf_I_Sacc",  ".csv"), Sacc, " ")
    #writedlm(string(dir, "_N_",  M, "_deltah_", @sprintf("%.3f",Δh), "_run_",  @sprintf("%3d", λX), "_Us_XX_Uf_I_Macc",  ".csv"), Macc, " ")
    #writedlm(string(dir, "_N_",  M, "_Jy_", @sprintf("%.3f", Θchar), "_t_", @sprintf("%.3f", tchar), "_delta_", @sprintf("%.3f",Δ), "_deltah_", @sprintf("%.3f",Δh), "_deltaphi_", @sprintf("%.3f",Δϕ), "_tf_", @sprintf("%.3f",λX), "_Macc_SPARSE", ".csv"), Macc, " ")
    #writedlm(string(dir, "_N_",  M, "_Jy_", @sprintf("%.3f", Θchar), "_t_", @sprintf("%.3f", tchar), "_delta_", @sprintf("%.3f",Δ), "_deltah_", @sprintf("%.3f",Δh), "_deltaphi_", @sprintf("%.3f",Δϕ), "_tf_", @sprintf("%.3f",λX), "_Rhof_SPARSE", ".csv"), ρ[now, :, :], " ")
    #writedlm("./_wfn.txt", map(abs, wfn), " ")
    #writedlm("./_dimers.txt", dimer, " ")

    return 0
end
