import numpy as np
from scipy.special import expi

class wellcalculator:
    def __init__(self):
        self.dt = 86400
        self.N = 300

    def calculate_flow(self, h, fi, k, mu, c, ro,  p_0, r_w, p_zab):      
        kf = k / (fi * mu * c)
        B = - (r_w**2) / (4 * kf * self.dt)
        Ei0 = -expi(B)
        A = 4 * kf * np.pi * fi * c * h / Ei0
        kV = np.arange(1, self.N + 1)
        tV = kV * self.dt
        Q = np.zeros(self.N)
        
        Ei = np.zeros(self.N)
        for n in range(self.N):
            current_k = n + 1
            Ei[n] = -expi(B/current_k) - (-expi(B/(current_k + 1)))

        E = Ei / Ei0
        
        Q[0] = A * (p_0 - pressure_profile[0]) 

        for n in range (1, self.N):
            ii = np.arange(0, n)
            ii_rev = ii[::-1]
            sum_term = np.dot(E[ii], Q[ii_rev])
            Q[n] = A * (p_0 - pressure_profile[n]) + sum_term

        return tV, Q, kf