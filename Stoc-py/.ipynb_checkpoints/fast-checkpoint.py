import numpy as np
import time
from finite import finite_difference
from scipy.stats import norm


def blsprice(r, v, K, S, T):
    # https://gist.github.com/carljohanrehn/4f20cc0433e10dc1c312
    d_1 = 1 / v / np.sqrt(T) * (np.log(S / K) + (r + v ** 2 / 2) * T)
    d_2 = d_1 - v * np.sqrt(T)
    return norm.cdf(d_1) * S - norm.cdf(d_2) * K * np.exp(-r * T)


r = 0.05
sigma = 0.2
K = 100
S_0 = 100
F = 2
dS = 5
T = 1/12

print("""Parameters
r:      {}
σ:      {}
K:      {}
S₀:     {}
Sₘₐₓ:   {}
dS:     {}
T:      {}""".format(r, sigma, K, S_0, F*S_0, dS, T))

s = time.time()
F, FF = finite_difference(r, sigma, K, S_0, F, dS, T)
BL = blsprice(r, sigma, K, S_0, T)
e = time.time()

print("""
Benchmark
time:   {}
shape:  {}

Results
finite: {}
theory: {}""".format(e - s, F.shape, FF, BL))
