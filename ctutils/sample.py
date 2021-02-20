import numpy as np

# input: uncertainty data file, number of samples, number of 
def sample_k_factor(filename, n):
    # load data
    d = np.genfromtxt(filename, dtype=np.dtype('i4,f8'), names='idx,f')
    # get dimension
    ndim = d.size
    #
    idx = d['idx'] - 1
    lnf = np.log(d['f'])
    # x = ln(k/k_o)/ln(f)
    x = np.random.normal(0,1,size=(ndim,n))
    # multiplier factors
    factors = np.exp(lnf * x.transpose())

    return idx, factors