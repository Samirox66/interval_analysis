import intvalpy as ip
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
# Task 1

def calculateDet(delta):
        midA = [[1.05, 1],
                [0.95, 1]]
        radA = [[delta, delta], [delta, delta]]
        A = ip.Interval(midA, radA, midRadQ=True)
        detA = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return detA

# Find delta that 0 is in det(A)
def task1():
    deltaArray = np.linspace(0, 1, num=10)
    detArrayHight = []
    detArrayLow = []
    for i in range(len(deltaArray)):
        delta = deltaArray[i]
        detA = calculateDet(delta)
        detArrayHight.append(detA.b)
        detArrayLow.append(detA.a)
    
    plt.figure()
    plt.title("Task 1")
    plt.grid()
    plt.plot(deltaArray, detArrayHight, label = "det up border")
    plt.plot(deltaArray, detArrayLow, label = "det down border")
    plt.plot(0.025, calculateDet(0.025).a, 'b*', label = "delta = 0.025")
    plt.xlabel('delta')
    plt.ylabel('det(A)')
    plt.legend()
    plt.show()


def det2(A):
    return A[0][0] * A[1][1] - A[1][0] * A[0][1]

def bauman(vertices):
    for i in range(len(vertices)):
        for j in range(i + 1):
            if det2(vertices[i]) * det2(vertices[j]) <= 0:
                return False
    return True

def task2():
    A = np.array([[1.05, 0.95], [ 1.0, 1.0]])
    
    deltas = np.linspace(0.04, 0.06, 10)
    # Регрессия
    print("\nЗадача регрессии")
    for delta in deltas:
        A1 = A.copy()
        A1[0][0] -= delta
        A1[1][0] -= delta
        A2 = A.copy()
        A2[0][0] += delta
        A2[1][0] -= delta
        A3 = A.copy()
        A3[0][0] += delta
        A3[1][0] += delta
        A4 = A.copy()
        A4[0][0] -= delta
        A4[1][0] += delta

        verticesA = [A1, A2, A3, A4]
        print("{:2.6f}\t{}".format(delta, "неособенная" if bauman(verticesA) else "особенная"))

    # Томография
    deltas = np.linspace(0.02, 0.029, 10)
    print("\nЗадача томографии")
    for delta in deltas:
        A1 = A.copy()
        A1[0][0] -= delta
        A1[0][1] -= delta
        A1[1][0] -= delta
        A1[1][1] -= delta

        A2 = A.copy()
        A2[0][0] += delta
        A2[0][1] -= delta
        A2[1][0] -= delta
        A2[1][1] -= delta
        A3 = A.copy()
        A3[0][0] -= delta
        A3[0][1] += delta
        A3[1][0] -= delta
        A3[1][1] -= delta
        A4 = A.copy()
        A4[0][0] -= delta
        A4[0][1] -= delta
        A4[1][0] += delta
        A4[1][1] -= delta
        A5 = A.copy()
        A5[0][0] -= delta
        A5[0][1] -= delta
        A5[1][0] -= delta
        A5[1][1] += delta

        A6 = A.copy()
        A6[0][0] += delta
        A6[0][1] += delta
        A6[1][0] -= delta
        A6[1][1] -= delta
        A7 = A.copy()
        A7[0][0] -= delta
        A7[0][1] += delta
        A7[1][0] += delta
        A7[1][1] -= delta
        A8 = A.copy()
        A8[0][0] -= delta
        A8[0][1] -= delta
        A8[1][0] += delta
        A8[1][1] += delta
        A9 = A.copy()
        A9[0][0] += delta
        A9[0][1] -= delta
        A9[1][0] -= delta
        A9[1][1] += delta
        A10 = A.copy()
        A10[0][0] += delta
        A10[0][1] -= delta
        A10[1][0] += delta
        A10[1][1] -= delta
        A11 = A.copy()
        A11[0][0] -= delta
        A11[0][1] += delta
        A11[1][0] -= delta
        A11[1][1] += delta


        A12 = A.copy()
        A12[0][0] += delta
        A12[0][1] += delta
        A12[1][0] += delta
        A12[1][1] -= delta
        A13 = A.copy()
        A13[0][0] -= delta
        A13[0][1] += delta
        A13[1][0] += delta
        A13[1][1] += delta
        A14 = A.copy()
        A14[0][0] += delta
        A14[0][1] -= delta
        A14[1][0] += delta
        A14[1][1] += delta
        A15 = A.copy()
        A15[0][0] += delta
        A15[0][1] += delta
        A15[1][0] -= delta
        A15[1][1] += delta

        A16 = A.copy()
        A16[0][0] += delta
        A16[0][1] += delta
        A16[1][0] += delta
        A16[1][1] += delta

        verticesA = [A1, A2, A3, A4, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16]
        print("{:2.5f}\t{}".format(delta, "неособенная" if bauman(verticesA) else "особенная"))


def task3():
    midA = [[1.05, 1],
            [0.95, 1]]
    U, S, V = np.linalg.svd(midA)
    print("Сингулярные числа " + str(S))


if __name__ == "__main__":
    print("Task 1:")
    task1()
    print("Task 2:")
    task2()
    print("Task 3:")
    task3()