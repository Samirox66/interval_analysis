import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt


def line(ax, coefs, res, x, y, color):
    if coefs[0] == 0 and coefs[1] == 0:
        return
    if coefs[1] == 0:
        ax.plot([res / coefs[0]] * len(y), y, color=color, linestyle='dashed')
    else:
        ax.plot(x, (res - coefs[0] * x) / coefs[1], color=color, linestyle='dashed')


def start_linear_system_plot(A, b):
    colors = ['pink', 'indigo', 'royalblue', 'deeppink']
    x = np.linspace(-1, 5, 100)
    y = [-0.5, 4]
    fig, ax = plt.subplots()
    for coefs, res, color in zip(A, b, colors):
        line(ax, coefs.mid, res.mid, x, y, color=color)


def linear_system_plot(A, b, title):
    start_linear_system_plot(A, b)
    plt.title(title)
    plt.grid()
    plt.show()


def start_tol_plot(A, b, needVe=False):
    x, y = np.mgrid[-1:5:100j, -0.5:3:45j]
    z = np.zeros(x.shape)
    for i in range(0, x.shape[0]):
        for j in range(0, x.shape[1]):
            z[i][j] = ip.linear.Tol(A, b, [x[i][j], y[i][j]])
    max = ip.linear.Tol(A, b, maxQ=True)

    fig, ax = plt.subplots()
    cs = ax.contour(x, y, z, levels=20, cmap='plasma')
    fig.colorbar(cs, ax=ax)
    ax.clabel(cs)
    ax.plot(max[1][0], max[1][1], 'r*', label=f'Максимум ({max[1][0]}, {max[1][1]}), значение: {max[2]}')
    if needVe:
        ive = ip.linear.ive(A, b)
        rve = ive * np.linalg.norm(b.mid) / np.linalg.norm([max[1][0], max[1][1]])
        print(f"ive: {ive}")
        print(f"rve: {rve}")
        iveRect = plt.Rectangle((max[1][0] - ive, max[1][1] - ive), 2 * ive, 2 * ive, edgecolor='deeppink', facecolor='none', label='Брус ive')
        plt.gca().add_patch(iveRect)
        rveRect = plt.Rectangle((max[1][0] - rve, max[1][1] - rve), 2 * rve, 2 * rve, edgecolor='navy', facecolor='none', label='Брус rve')
        plt.gca().add_patch(rveRect)
    return max[2]


def start_tol_plot_with_lines(A, b, needVe=False):
    x, y = np.mgrid[-1:5:200j, -0.5:4:90j]
    z = np.zeros(x.shape)
    for i in range(0, x.shape[0]):
        for j in range(0, x.shape[1]):
            z[i][j] = ip.linear.Tol(A, b, [x[i][j], y[i][j]])
    max = ip.linear.Tol(A, b, maxQ=True)

    fig, ax = plt.subplots()
    cs = ax.contour(x, y, z, levels=20, cmap='plasma')
    fig.colorbar(cs, ax=ax)
    ax.clabel(cs)
    ax.plot(max[1][0], max[1][1], 'r*', label='Максимум ({}, {}), значение: {}'.format(max[1][0], max[1][1], max[2]))
    if needVe:
        ive = ip.linear.ive(A, b)
        rve = ive * np.linalg.norm(b.mid) / np.linalg.norm([max[1][0], max[1][1]])
        print("ive: {}".format(ive))
        print("rve: {}".format(rve))
        iveRect = plt.Rectangle((max[1][0] - ive, max[1][1] - ive), 2 * ive, 2 * ive, edgecolor='aqua', facecolor='none', label='Брус ive')
        plt.gca().add_patch(iveRect)
        rveRect = plt.Rectangle((max[1][0] - rve, max[1][1] - rve), 2 * rve, 2 * rve, edgecolor='coral', facecolor='none', label='Брус rve')
        plt.gca().add_patch(rveRect)

    colors = ['pink', 'indigo', 'royalblue', 'deeppink']
    x = np.linspace(-1, 5, 100)
    y = [-0.5, 4]
    for coefs, res, color in zip(A, b, colors):
        line(ax, coefs.mid, res.mid, x, y, color=color)

    return max[2]


def end_plot(title):
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()


def tol_plot(A, b, title, needVe=False):
    tol = start_tol_plot(A, b, needVe)
    end_plot(title)
    return tol


def b_correction(b, K, weights):
    return b + K * ip.Interval(-1, 1) * weights


def A_correction(A, K, weights, E):
    mul = K * weights * E
    newA = A.a - mul.a
    newB = A.b - mul.b
    return ip.Interval(newA, newB)