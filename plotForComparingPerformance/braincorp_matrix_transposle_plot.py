import matplotlib.pyplot as plt
import numpy as np

x_axis = np.arange(1, 8, 1)

# collected from running matrix transpose
tiling_bf_ratio = [1.000000, 1.250000, 1.785714, 1.300000, 1.866142, 1.651572, 1.644214]
tiling_thread_bf_ratio = [3.000000, 2.500000, 6.250000, 4.333333, 6.155844, 8.931973, 9.023697]

# Plot the results
plt.figure()
plt.scatter(x_axis, tiling_bf_ratio, c="darkorange", label="tiling/bruteforce")
plt.scatter(x_axis, tiling_thread_bf_ratio, c="green", label="tiling_thread/bruteforce")

plt.plot(x_axis, tiling_bf_ratio, color="cornflowerblue",  linewidth=2)
plt.plot(x_axis, tiling_thread_bf_ratio, color="red", linewidth=2)
plt.xlabel("Matrices with growing dimentions")
plt.ylabel("Fold increase in efficiency over brute force")
plt.title("tiling+thread and tiling only both have increased efficiency for transpose")
plt.legend()
plt.show()