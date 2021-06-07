import math
from sys import argv

def sigmoid(c, lam):
    x = c[0]
    for i in range(1, len(c)):
        x = x * lam + c[i]

    y = 1.0 / math.sqrt(x * x + 1)
    return 0.5 * x * y + 0.5


min_lambda = 380
max_lambda = 780

c = []
for i in range(1, len(argv)):
    c.append(float(argv[i]))

for i in range(min_lambda, max_lambda + 1, 5):
    c_l = (i - min_lambda) / (max_lambda - min_lambda)
    print(i, sigmoid(c, c_l))
