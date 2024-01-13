using ModalIntervalArithmetic
using Plots

include("sub_diff.jl")

function task_1()
    A = [ModalInterval(3.0, 4.0) ModalInterval(5.0, 6.0)
         ModalInterval(-1.0, 1.0) ModalInterval(-3.0, 1.0)]
    b = [ModalInterval(-2.7, 3.9), ModalInterval(-1.0, 1.7)]

    x0 = [0.0, 1.0, 0.0, 1.0]#init_point(A, b)
    println("x0    = ", x0)
    (x, count) = sub_diff(A, b, x0, 1.e-9)
    println("x     = ", x)
    println("count = ", count)

    # n = 8
    # iters = []
    # for k = 1:n
    #     prec = 10.0^-k
    #     (x, count) = sub_diff(A, b, x0, prec)
    #     push!(iters, count)
    # end
    # plot(1:n, iters,
    #     xlabel = "-lg(ε)",
    #     ylabel = "Iterations number",
    #     xticks = [0:1:n;],
    #     yticks = [0:2:max(count);],
    #     dpi=300,
    #     legend = false)
    # savefig("plot.png")
end
