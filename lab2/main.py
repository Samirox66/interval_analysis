import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def line(ax, coefs, res, x, y, mod):
    if coefs[0] == 0 and coefs[1] == 0:
        return 
    if coefs[1] == 0:
        ax.plot([res / coefs[0]] * len(y), y, mod, linewidth=1.5)
    else:
        ax.plot(x, (res - coefs[0] * x) / coefs[1], mod, linewidth=1.5)


def start_linear_system_plot(A, b):
    colors = ['c', 'y', 'r', 'g']
    x = np.linspace(-1, 5, 100)
    y = [-0.5, 4]
    fig, ax = plt.subplots()
    for coefs, res, color in zip(A, b, colors):
        line(ax, coefs.mid, res.mid, x, y, color + '-')
        line(ax, coefs.a, res.a, x, y, color + '-')
        line(ax, coefs.b, res.b, x, y, color + '-')


def linear_system_plot(A, b, title):
    start_linear_system_plot(A, b)
    plt.title(title)
    plt.grid()
    plt.xlim(-1, 5)
    plt.ylim(-0.5, 4)
    plt.show()

def draw_res(ax, x):
    ax.plot([x[0], x[0]], [x[2], x[3]], color="r")
    ax.plot([x[1], x[1]], [x[2], x[3]], color="r")
    ax.plot([x[0], x[1]], [x[2], x[2]], color="r")
    ax.plot([x[0], x[1]], [x[3], x[3]], color="r")

def start_tol_plot(A, b, needVe=False):
    # x, y = np.mgrid[-1:5:100j, -1:4:45j]
    # z = np.zeros(x.shape)
    # for i in range(0, x.shape[0]):
    #     for j in range(0, x.shape[1]):
    #         z[i][j] = ip.linear.Tol(A, b, [x[i][j], y[i][j]])
    # max = ip.linear.Tol(A, b, maxQ=True)
    fig, ax = plt.subplots(figsize=(9, 5.7))
    # cs = ax.contour(x, y, z, levels = 20, cmap=cm.gray)
    # fig.colorbar(cs, ax=ax)
    # ax.clabel(cs)

    x, y = np.mgrid[-1:1:40j, -0.7:0.5:40j]
    z = np.zeros(x.shape)
    for i in range(0, x.shape[0]):
        for j in range(0, x.shape[1]):
            f = ip.linear.Tol(A, b, [x[i][j], y[i][j]])
            if(f >= 0):
                z[i][j] = f
    ax.contour(x, y, z, levels = 0, linewidths=1, colors='#178913')
    x1 = [-0.18, 0.35, 0.95, 0.02]
    x2 = [0.14, 0.5, 0.5, 0.33]
    # x3 = [0.15, 0.35, 0.95, 0.02]
    # x4 = [0.3, 0.25, 0.65, 0.22]

    draw_res(ax, x1)
    draw_res(ax, x2)
    

    # colors = ['c', 'y', 'r', 'g']
    # x = np.linspace(-1, 5, 100)
    # y = [-0.5, 4]
    # for coefs, res, color in zip(A, b, colors):
    #     line(ax, coefs.mid, res.mid, x, y, color + '-')
    #     line(ax, coefs.a, res.a, x, y, color + '-')
    #     line(ax, coefs.b, res.b, x, y, color + '-')

    # ax.plot(max[1][0], max[1][1], 'b*', label='Максимум ({:.4f}, {:.4f}), значение: {:.4f}'.format(max[1][0], max[1][1], max[2]), markersize=3)
    
    # if needVe:
    #     ive = ip.linear.ive(A, b)
    #     rve = ive * np.linalg.norm(b.mid) / np.linalg.norm([max[1][0], max[1][1]])
    #     print("ive: {}".format(ive))
    #     print("rve: {}".format(rve))
    #     iveRect = plt.Rectangle((max[1][0] - ive, max[1][1] - ive), 2 * ive, 2 * ive, edgecolor='black', facecolor='none', label='Брус ive', linewidth=1.5)
    #     plt.gca().add_patch(iveRect)
    #     rveRect = plt.Rectangle((max[1][0] - rve, max[1][1] - rve), 2 * rve, 2 * rve, edgecolor='#7E2F8E', facecolor='none', label='Брус rve', linewidth=1.5)
    #     plt.gca().add_patch(rveRect)

    return 2#max[2]


def tol_plot(A, b, title, needVe=False):
    tol = start_tol_plot(A, b, needVe)
    plt.title(title)
    # plt.legend()
    # plt.grid()
    plt.xlim(-1, 1.2)
    plt.ylim(-0.7, 0.5)

    return tol


def b_correction_uneven(b, K, weights):
    return b + K * ip.Interval(-1, 1) * weights

def b_correction_even(b, K):
    return b + K * ip.Interval(-1, 1)


def A_correction(A, K, weights, E):
    mul = K * weights * E
    newA = A.a - mul.a
    newB = A.b - mul.b
    return ip.Interval(newA, newB)


midA = np.array([[3,5], [3, -5], [1, 0], [0, 1]])
radA = np.array([[1, 1], [0, 1], [0.5, 0], [0, 0.5]])
A = ip.Interval(midA, radA, midRadQ=True)

midb = np.array([7, 0, 2, 1])
radb = np.array([2, 1, 1, 1])
b = ip.Interval(midb, radb, midRadQ=True)

print(np.linalg.cond(midA))

# plot_brus(A, b, 'tet')

            # wb_1 = [0.2, 0.3, 1, 0.2]
            # x = np.zeros(A.shape[1])
            # print(min(wb_1*(b.rad - abs(b.mid - (A @ x)))))



# linear_system_plot(A, b, 'начальная ИСЛАУ')
# maxTol = tol_plot(A, b, 'Tol для начальной ИСЛАУ')
# print("maxTol: {}\n".format(maxTol))

my_A = ip.Interval([
    [[3, 4], [5, 6]],
    [[-1, 1], [-3, 1]],
])

my_b = ip.Interval([[-2.7, 3.9], [-1, 1.7]])

tol_plot(my_A, my_b, "bachenich", False)

# K = 1
# bEvenCorrected = b_correction_even(b, K)
# print("Равномерное уширение b: {}".format(bEvenCorrected.data))
# maxTolBCorrected = tol_plot(A, bEvenCorrected, 'Tol для системы с равномерным уширением правой части', True)
# print("maxTol для равномерного уширения: {}\n".format(maxTolBCorrected))

# K = 3
# weightsB = [0.5, 0.1, 1, 0.5]
# bUnevenCorrected = b_correction_uneven(b, K, weightsB)
# print("Неравномерное уширение b: {}".format(bUnevenCorrected.data))
# maxTolBCorrected = tol_plot(A, bUnevenCorrected, 'Tol для системы с неравномерным уширением правой части', True)
# print("maxTol для неравномерного уширения: {}\n".format(maxTolBCorrected))

# K = 1
# midE = np.zeros((4, 2))
# radE = np.array([[0.5, 0.4], [0.3, 0.5], [0.5, 0], [0, 0.3]])
# E = ip.Interval(midE, radE, midRadQ=True)
# weightsA = np.ones((4, 2))
# ACorrected = A_correction(A, K, weightsA, E)
# print("Изначальная матрица A: \n{}\n".format(A.data))
# print("Скорректированная матрица A: \n{}".format(ACorrected.data))
# maxTolACorrected = tol_plot(ACorrected, b, 'Tol для системы со скорректированной левой частью', True)


plt.show()