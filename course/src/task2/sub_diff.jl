using LinearAlgebra

include("calcs.jl")

function init_point(C, d)
    midC = mid.(C)
    C̃ = [pos.(midC) neg.(midC)
         neg.(midC) pos.(midC)]
    return C̃ \ sti(d)
end

function sub_diff(C, d, x0, eps)
    𝒢(x) = sti(C * sti_inv(x)) - sti(d)

    x = x0
    𝒢_val = 𝒢(x)
    count = 0

    while(norm(𝒢_val) >= eps)
        println("x ", x)
        try
            x -= inv(D(C, x)) * 𝒢_val
        catch;
            println("Subgradient D is singular")
            break
        end
        𝒢_val = 𝒢(x)
        count += 1
    end

    return (sti_inv(x), count)
end
