#tchar = 0.30901699437494745

    #println(parse(Int64, ARGS[1]))
    #println(parse(Float64, ARGS[2]))
    include("z2_prod_D_1_RANDOM.jl");
    Δ  = [0.0]
    Δt = 0.30901699437494745 * [1.0]
    Δϕ = [0.0]
    #Δh = [0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.321751,0.34,0.36,0.387699,0.41,0.44,0.463648,0.5,0.53,0.559885,0.63,0.7,0.785398]
    Δh = [0.0]

    n = parse(Int64, ARGS[1])
    h = parse(Float64, ARGS[2])
    run = parse(Int64, ARGS[3])

    for δ in Δ
        for δh in Δh
            for δϕ in Δϕ
                for δt in Δt
                    main_D_1(2, n, 1, 1000, δt, δ, δh, δϕ, h, run, "./D_1_prod_init/");
                    #print("done: N = ", parse(Int64, ARGS[1]), ", δ = ", δ, ", δh = ", δh, ", δϕ = ", δϕ, ", δt = ", δt, ", λ = ", parse(Float64, ARGS[2]), "\n")
                end
            end
        end
    end
