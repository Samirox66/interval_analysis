import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def draw_res(ax, x):
    ax.plot([x[0], x[0]], [x[2], x[3]], color="r")
    ax.plot([x[1], x[1]], [x[2], x[3]], color="r")
    ax.plot([x[0], x[1]], [x[2], x[2]], color="r")
    ax.plot([x[0], x[1]], [x[3], x[3]], color="r")
    ax.plot([0, 0], [-0.7, 0.5], color="k", linewidth=1)
    ax.plot([-1.0, 1.5], [0, 0], color="k", linewidth=1)

def start_tol_plot(A, b, needVe=False):
    fig, ax = plt.subplots(figsize=(9, 5.7))

    x, y = np.mgrid[-1:1.2:400j, -0.7:0.5:400j]
    z = np.zeros(x.shape)
    for i in range(0, x.shape[0]):
        for j in range(0, x.shape[1]):
            f = ip.linear.Tol(A, b, [x[i][j], y[i][j]])
            if(f >= 0):
                z[i][j] = f
    ax.contour(x, y, z, levels = 0, linewidths=1, colors='#178913')
    x1 = [-0.255, 0.75, -0.25, 0.085]
    x2 = [-0.15, 0.8, -0.4, 0.133]

    draw_res(ax, x2)
    #draw_res(ax, x2)

    return


def plot(A, b, title, needVe=False):
    tol = start_tol_plot(A, b, needVe)
    plt.title(title)
    plt.xlim(-1, 1.5)
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

my_A = ip.Interval([
    [[3, 4], [5, 6]],
    [[-0.75, 1], [-3, 1]],
])

my_b = ip.Interval([[-3, 4], [-1, 2]])

plot(my_A, my_b, "Графичек", False)

plt.show()