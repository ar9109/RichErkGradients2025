# mymod.py
"""Python module demonstrates passing MATLAB types to Python functions"""
def fit_func_py(B,x):
    import numpy as np
    B = np.asarray(B,'float64')
    x = np.asarray(x,'float64')
    return B[0]*x ** B[1]