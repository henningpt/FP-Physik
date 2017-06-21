import numpy as np

arr1 = np.array([1.0, 4.3, -2.1, 2.0])
arr2 = np.array([0, 1, 2, 3])
arr1_p = arr1[arr1 > 0]
arr2_p = arr2[arr1 > 0]


print("arr1: ", arr1)
print("arr2: ", arr2)
print("arr1_p: ", arr1_p)
print("arr2_p: ", arr2_p)
