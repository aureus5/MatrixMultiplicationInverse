import matplotlib.pyplot as plt
import numpy as np

x_axis = np.arange(1, 8, 1)

# collected from running matrix transpose
tiling_bf_ratio = [0.777778, 1.026936, 0.963702, 1.074975, 1.060190, 2.727796, 3.195113]
tiling_thread_bf_ratio = [0.777778, 3.112245, 1.382812, 2.637028, 1.457329, 9.538142, 11.271130]

# Plot the results
plt.figure()
plt.scatter(x_axis, tiling_bf_ratio, c="darkorange", label="tiling/bruteforce")
plt.scatter(x_axis, tiling_thread_bf_ratio, c="green", label="tiling_thread/bruteforce")

plt.plot(x_axis, tiling_bf_ratio, color="cornflowerblue",  linewidth=2)
plt.plot(x_axis, tiling_thread_bf_ratio, color="red", linewidth=2)
plt.xlabel("Matrices with growing dimentions")
plt.ylabel("Fold increase in efficiency over brute force")
plt.title("tiling+thread and tiling only both have increased efficiency for multiply")
plt.legend()
plt.show()