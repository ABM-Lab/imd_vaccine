# This file implements the Bayesian model found at 
# https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-021-06906-x#Sec2
# https://doi.org/10.1186/s12879-021-06906-x

using Turing, MCMCChains, MCMCDiagnosticTools
using Random 

function get_data() 
    cc, pp = get_raw_data()
    case_counts = Array{Int64}(undef, 8, 21, 2)
    case_counts[:, :, 1] = cc[:, :, 1]
    case_counts[:, :, 2] = cc[:, :, 2] + cc[:, :, 3]

    pop_counts = Array{Int64}(undef, 8, 21, 2)
    pop_counts[:, :, 1] = pp[:, :, 1]
    pop_counts[:, :, 2] = pp[:, :, 2] + pp[:, :, 3]

    return case_counts, pop_counts
end 

function get_consol_data() 
    # merges age groups 4 and 5 together
    C, K = get_data() 
    a, t, _ = size(C)

    case_counts = Array{Int64}(undef, 7, 21, 2)
    pop_counts = Array{Int64}(undef, 7, 21, 2)

    for d in 1:2 
        case_counts[:, :, d] .=  [C[1:3, :, d]; sum(C[4:5, :, d]; dims=1); C[6:end, :, d]]
        pop_counts[:, :, d] .= [K[1:3, :, d]; sum(K[4:5, :, d]; dims=1); K[6:end, :, d]]
    end
    case_counts, pop_counts
end

function get_raw_data() 
    case_counts = Array{Int64}(undef, 8, 21, 3)
    case_counts[:, :, 1] .= 
    [79 77 48 22 31 36 49 47 23 34 19 18 20 8 8 6 4 9 8 6 2;
    95 46 47 36 37 31 43 42 21 32 26 9 14 17 10 17 5 15 13 7 8;
    49 48 24 21 29 24 21 27 15 18 16 8 5 2 1 4 5 5 6 7 2;
    96 89 66 39 31 30 23 21 14 5 1 3 0 0 1 2 2 1 2 0 0;
    181 161 156 101 81 96 89 70 71 72 53 32 22 7 8 12 2 2 4 6 4;
    251 185 207 98 108 121 106 111 101 96 99 78 92 61 47 49 46 38 53 28 46;
    175 64 104 44 65 72 69 87 74 65 71 56 37 41 23 45 30 36 46 35 24;
    268 171 165 119 93 107 97 137 107 101 108 63 71 51 40 42 38 46 58 27 21;]

    case_counts[:, :, 2] .= 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 2 3 3 1 1 2 1 2 0 1 0 0 2 1 0;
    0 0 0 0 0 2 5 5 4 6 5 5 4 3 2 2 2 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 2 1 1 1 4 2 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]


    case_counts[:, :, 3] .= 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 1 2 0 1 3 3 0 1 2 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 1 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]
    

    pop_counts = Array{Int64}(undef, 8, 21, 3)
    pop_counts[:, :, 1] .= 
    [4012658 3951461 3975871 4014258 4004393 4041738 4147997 4132735 4003587 3951495 3963264 3926731 3931411 3954973 3984144 3963268 3882437 3826908 3762227 3735010 3564493;
    15285559 15477731 15616575 15771627 15913007 15897145 15977965 16138392 16240931 16238083 16162947 16054205 15923773 15923833 15940562 15973469 16012579 15951619 15809112 15566282 15262845;
    24450611 24116233 23825893 23606188 23539189 23665793 23790329 24023282 24315589 24518896 24471257 24519060 24649877 24642601 24611852 24597701 24524520 24399364 24288025 24311884 24395037;
    20785307 21115481 21382965 21524592 18788689 17196407 13749796 11465642 9398806 7802853 6840326 6368484 5686098 5404001 5093124 4813640 4299102 4112247 3827749 3960208 4207511;
    32380312 32936731 33272014 33515402 33379208 33126673 31117522 28895913 25831596 22856343 19654066 16252214 12992600 11212849 9384969 8005308 6703566 5995153 5269497 4556754 4102840;
    108971637 109062111 109064863 109299940 109389704 109555557 109507695 109376382 109318696 109071085 108781933 108636495 108384675 106793925 105299011 103633220 101657062 98994345 95982707 92659702 89701182;
    43792580 45443238 47106223 48869972 50720230 52500986 54268612 55796537 57410443 59134390 60632510 61147520 61826041 62495048 62997741 63188106 63203451 63075500 62901109 62799204 63858481;
    35290291 35522207 35863529 36203319 36649798 37164107 37825711 38777621 39623176 40478249 41350891 43132211 44632337 46161005 47655870 49208479 50757639 52354605 54036735 55659365 54438296;]

    pop_counts[:, :, 2] .= 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 2727780 4119230 7377263 9463594 11385896 12919363 13894944 14400245 15042857 15315123 15647472 15871905 16400713 16727076 17101841 16893567 17490297;
    0 0 0 0 399594 1106499 3462000 6016819 9235048 12347940 13576037 15822932 17616028 17541662 17446300 16615847 16422668 15557280 14859855 13975351 13396257;
    0 0 0 0 0 0 0 0 0 0 57316 205995 591423 2592408 4704014 7185327 9045913 11395575 13803343 16071273 18199514;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]
    

    pop_counts[:, :, 3] .= 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 2957 2911 2760 2721 2679 2647 2606 2565 2523 2505 2490 2453 2436 2434 2433 2395 2413;
    0 0 0 0 2050 2866 3557 4328 5083 5799 2194706 3407611 4778199 6343533 7967826 10008670 11347959 12684020 13955059 15488866 16953753;
    0 0 0 0 0 0 0 0 0 0 678 1394 2105 2863 3619 4362 862083 1762073 2730261 3804262 5059660;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]
    return case_counts, pop_counts
end

@model function model(data, K)
    max_ag, max_t, max_d = size(data) # maximum age group, time, vaccine  
    ρ ~ filldist(Normal(-10, 0.5), max_ag)
    b ~ filldist(Normal(0, 0.5), max_t)
    σ ~ filldist(Truncated(Normal(0, 1), 0, Inf), max_t)
    for i = 1:max_t
        b[i] ~ Normal(0, σ[i]) # overwrite b prior distributions
    end
    θ₁ ~ truncated(Normal(-1.171, 10), -Inf, 0)
    θ₂ ~ truncated(Normal(-1.171, 10), -Inf, 0)
    θ = [θ₁, θ₂]
    for d = 1:max_d
        for t = 1:max_t
            for ag = 1:max_ag
                if d == 1 # no vaccine 
                    λ = log(K[ag, t, d]) + ρ[ag] + b[t]
                    data[ag, t, d] ~ Poisson(exp(λ)) 
                else 
                    if K[ag, t, d] > 0 
                        λ = log(K[ag, t, d]) + ρ[ag] + b[t] + θ[d-1]
                        data[ag, t, d] ~ Poisson(exp(λ)) 
                    else 
                        data[ag, t, d] ~ Poisson(0) 
                    end 
                end
            end
        end
    end
end

function run_model(;progress=false, iter=5000, ag45=false) 
    Random.seed!(1418)
    println("running model - params: iter: $iter, ag45: $ag45")
    if ag45  # we don't report this in the paper anymore so this branch is never run
        case_counts, pop_counts = get_consol_data()
    else 
        case_counts, pop_counts = get_raw_data() #get_data()
    end
    model_fun = model(case_counts, pop_counts);
    sampler = NUTS(1000, 0.65)
    chain = sample(model_fun, sampler, iter; progress=progress)
    return chain, case_counts, pop_counts
end

function infer(modelresults) 
    # describe(group(modelresults, :data)) .|> DataFrame
    _, pop_counts = get_data()    
    case_counts = Array{Union{Missing, Float64}}(undef, 8, 21, 2)
    prior_pred = predict(model(case_counts, pop_counts), modelresults)
    return prior_pred

    # infer's the model output 
    # I think this samples from the posterior distribution of all of parameters 
    # and then generates the quantities of interest
    # _df = infer(modelresults)
    # meandf = describe(_df)[1] |> DataFrame
    # quantdf = describe(_df)[2] |> DataFrame
    # @select!(meandf, :parameters, :mean)
    # mean_and_quants = outerjoin(meandf, quantdf[:, [1, 2, 6]], on=:parameters)
    # rename!(mean_and_quants,[:parameters,:mean, :lower, :upper])
    #CSV.write("model_theta_on.csv", mean_and_quants)

    # DONT SEND SEYED THIS DATA -- THE MANUAL CALCULATION IS SLIGHTLY DIFFERENT (For Quantiles)
end
