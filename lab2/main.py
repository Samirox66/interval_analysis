from linear_system import *

#midA = np.array([[1, 1.5], [1, -2], [1, 0], [0, 1]])
#radA = np.array([[0.5, 1], [0, 1], [0.1, 0], [0, 0.1]])
midA = np.array([[1, 1.5], [1, -2], [1, 0], [0, 1]])
radA = np.array([[0.5, 1], [0, 1], [0.1, 0], [0, 0.1]])
A = ip.Interval(midA, radA, midRadQ=True)
print(f"A values: {A}")

# b1, b2 = np.random.uniform(1, 5), np.random.uniform(1, 5)
b1, b2 = 3.2, 1.8
print(f"b values: {b1}, {b2}")
midb = np.array([5, 0, 4, 1])
radb = np.array([2, 0.5, 1, 1])
b = ip.Interval(midb, radb, midRadQ=True)

linear_system_plot(A, b, 'mid исходной СЛАУ')
maxTol = tol_plot(A, b, 'Tol для исходной системы')
K = 1.5 * maxTol
print(f"K: {K}")

weightsB = np.ones(len(b))
bCorrected = b_correction(b, K, weightsB)
print(f"corrected b: {bCorrected}")
maxTolBCorrected = tol_plot(A, bCorrected, 'Tol для системы со скорректированной правой частью', True)

midE = np.zeros((4, 2))
radE = np.array([[0.3, 0.6], [0, 0.6], [0.06, 0], [0, 0.06]])
E = ip.Interval(midE, radE, midRadQ=True)
weightsA = np.ones((4, 2))
ACorrected = A_correction(A, K, weightsA, E)
print(f"corrected A: {ACorrected}")
maxTolACorrected = tol_plot(ACorrected, b, 'Tol для системы со скорректированной левой частью', True)

start_tol_plot_with_lines(ACorrected, b)
end_plot('Tol и mid для системы со скорректированной левой частью')

ACorrected1 = ip.Interval(midA, radA, midRadQ=True)
row = 0
ACorrected1[row, 0] = ip.Interval(midA[row, 0], midA[row, 0])
ACorrected1[row, 1] = ip.Interval(midA[row, 1], midA[row, 1])
start_tol_plot_with_lines(ACorrected1, b)
end_plot('Tol и mid для системы со скорректированной первой строкой')
print(f"corrected A1: {ACorrected1}")

ACorrected2 = ip.Interval(midA, radA, midRadQ=True)
row = 1
ACorrected2[row, 0] = ip.Interval(midA[row, 0], midA[row, 0])
ACorrected2[row, 1] = ip.Interval(midA[row, 1], midA[row, 1])
start_tol_plot_with_lines(ACorrected2, b)
end_plot('Tol и mid для системы со скорректированной второй строкой')
print(f"corrected A2: {ACorrected2}")

ACorrected3 = ip.Interval(midA, radA, midRadQ=True)
row = 2
ACorrected3[row, 0] = ip.Interval(midA[row, 0], midA[row, 0])
ACorrected3[row, 1] = ip.Interval(midA[row, 1], midA[row, 1])
start_tol_plot_with_lines(ACorrected3, b)
end_plot('Tol и mid для системы со скорректированной третьей строкой')
print(f"corrected A3: {ACorrected3}")

ACorrected4 = ip.Interval(midA, radA, midRadQ=True)
row = 3
ACorrected4[row, 0] = ip.Interval(midA[row, 0], midA[row, 0])
ACorrected4[row, 1] = ip.Interval(midA[row, 1], midA[row, 1])
start_tol_plot_with_lines(ACorrected4, b)
end_plot('Tol и mid для системы со скорректированной четвёртой строкой')
print(f"corrected A4: {ACorrected3}")