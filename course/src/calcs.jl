using ModalIntervalArithmetic

pos(x) = x > 0 ? x : 0
neg(x) = x < 0 ? -x : 0

function sti(𝐱)
    n = size(𝐱, 1)
    x = zeros(2 * n)
    for i = 1:n
        x[i] = -𝐱[i].inf
        x[i+n] = 𝐱[i].sup
    end
    return x
end

function sti_inv(x)
    n = size(x, 1) ÷ 2
    𝐱 = zeros(ModalInterval, n)
    for i = 1:n
        𝐱[i] = ModalInterval(-x[i], x[i+n])
    end
    return 𝐱
end

function ∂x_pos(x)
    if x > 0
        return 1
    elseif x == 0
        return 0.5 # actually something ∈ [0, 1]
    else
        return 0
    end
end

function ∂x_neg(x)
    if x < 0
        return -1
    elseif x == 0
        return -0.5 # actually something ∈ [-1, 0]
    else
        return 0
    end
end

function ∂max_1(C, i, j, x)
    n = size(x, 1) ÷ 2
    prod_1 = pos(C[i, j].sup)pos(x[j])
    prod_2 = neg(C[i, j].inf)pos(x[j+n])
    if prod_1 > prod_2
        return (pos(C[i, j].sup), 0)
    elseif prod_1 == prod_2
        return (0.5pos(C[i, j].sup), 0.5neg(C[i, j].inf))
    else
        return (0, neg(C[i, j].inf))
    end
end

function ∂max_2(C, i, j, x)
    n = size(x, 1) ÷ 2
    prod_1 = pos(C[i, j].sup)pos(x[j+n])
    prod_2 = neg(C[i, j].inf)pos(x[j])
    if prod_1 > prod_2
        return (0, pos(C[i, j].sup))
    elseif prod_1 == prod_2
        return (0.5neg(C[i, j].inf), 0.5pos(C[i, j].sup))
    else
        return (neg(C[i, j].inf), 0)
    end
end

function ∂F(C, i, x)
    n = size(x, 1) ÷ 2
    res = zeros(2 * n)
    if 1 <= i <= n
        for j = 1:n
            temp = ∂max_1(C, i, j, x)
            res_1 = pos(C[i, j].inf)∂x_neg(x[j]) + neg(C[i, j].sup)∂x_neg(x[j+n]) - temp[1]
            res_2 = pos(C[i, j].inf)∂x_neg(x[j]) + neg(C[i, j].sup)∂x_neg(x[j+n]) - temp[2]
            res[j] -= res_1
            res[j+n] -= res_2
        end
    else
        i -= n
        for j = 1:n
            temp = ∂max_2(C, i, j, x)
            res_1 = temp[1] - pos(C[i, j].inf)∂x_neg(x[j+n]) - neg(C[i, j].sup)∂x_neg(x[j])
            res_2 = temp[2] - pos(C[i, j].inf)∂x_neg(x[j+n]) - neg(C[i, j].sup)∂x_neg(x[j])
            res[j] += res_1
            res[j+n] += res_2
        end
    end
    return res
end

function D(C, x)
    n = size(x, 1)
    D = zeros(n, n)
    for i = 1:n
        D[i, :] = ∂F(C, i, x)
    end
    return D
end
