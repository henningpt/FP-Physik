import numpy as np
arr = np.array([1, 0, 2, 3])
arr_2 = np.array([0, 1, 2, 3])

arr_new = arr_2[arr > 0]
print(arr_new)
