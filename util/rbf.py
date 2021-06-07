import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import product


def kernel(x):
    # return np.abs(x)**2 * np.log(np.abs(x + 0.000001))
    return np.exp(-(1000 * x) ** 2)


def compute_distance_matrix(x):
    distance_matrix = np.zeros((len(x), len(x)))

    for i in range(len(x)):
        for j in range(len(x)):
            dx = x[i, 0] - x[j, 0]
            dy = x[i, 1] - x[j, 1]
            distance_matrix[i, j] = math.sqrt(dx * dx + dy * dy)

    return distance_matrix


class Rbf:
    def __init__(self, xs, fs):
        self.xs = xs
        self.fs = fs
        self.weights = None

    def solve(self):
        # build kernel matrix
        phi = compute_distance_matrix(self.xs)
        phi = kernel(phi)

        # build and solve lgs
        self.weights = np.matmul(np.linalg.inv(phi), self.fs)

    def eval(self, x):
        phi = np.full(self.xs.shape, x) - self.xs
        phi = np.linalg.norm(phi, axis=-1)
        phi = kernel(phi)
        return np.matmul(phi, self.weights)


# define function we want to interpolate

# f_xs = np.array([[x, y] for x, y in product(np.linspace(0, 4, 5), np.linspace(0, 4, 5))])
f_xs = np.array([
    [0, 0],
    [0, 1],
    [1, 0],
    [1, 1]
])
print(f_xs)
f_ys = np.array([
    [0, 0],
    [1, 1],
    [2, 0],
    [0, 1]
])
print(f_ys)

print(compute_distance_matrix(f_xs))

rbf = Rbf(f_xs, f_ys)
rbf.solve()

test_xs = np.array([
    [-100, -100],
    [-100, 0],
    [0, -100],
    [0, 100],
    [100, 0],
    [100, 100]
])

for x in f_xs:
    print(x, rbf.eval(x))

for x in test_xs:
    print(x, rbf.eval(x))


# sample_xs = np.linspace(0, 40, 200)
# sample_ys = np.array([rbf.eval(x) for x in sample_xs])

# plt.plot(sample_xs, sample_ys)
# plt.show()
