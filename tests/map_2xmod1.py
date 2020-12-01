import numpy as np
X = np.zeros(10000)
for i in range(10000):
    if i == 0:
        X[i] = 0.1
    else:
        X[i] = (2 * X[i - 1]) % 1

print(X)

import matplotlib.pyplot as plt
plt.plot(X)
plt.ylabel('some numbers')
plt.show()
