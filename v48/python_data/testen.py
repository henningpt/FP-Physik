import numpy as np


def rate(arr):
    val = 0.0
    for i in range(1, len(arr)):
        val += abs(arr[i] - arr[i-1])
    return(val / (len(arr) - 1))


arr1 = np.array([1.0, 2.0, 4.0, 7.0])
print("array 1: ", arr1)

print("Raten: ", rate(arr1))
print("Raten2: ", abs(arr1[0] - arr1[len(arr1) - 1]) / (len(arr1) - 1))
