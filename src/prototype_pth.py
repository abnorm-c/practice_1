import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expi
import pandas as pd

plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['font.size'] = 12

h = float(input("Высота (м): "))
fi = float(input("Пористость (ед): "))
k = float(input("Проницаемость (м^2): "))
mu = float(input("Вязкость динамическая (Па * с): "))
ro = float(input("Плотность (кг/м^3): "))
c = float(input("сжимаемость полная (Па^-1): "))
r_w = float(input("радиус скважины (м): "))
p_0 = float(input("Начальное давление (Па): "))
p_zab = float(input("Давление забойное (Па): "))

p_inj = p_0

if any(param <=0 for param in [h, fi, k, mu, ro, c, r_w, p_0, p_zab]):
    print("n/a")

class wellcalculator:
    def __init__(self):
        self.dt = 86400
        self.N = 300

    def calculate_flow(self, h, fi, k, mu, c, p_0, r_w, p_zab):
        kf = k / (fi * mu * c)
        print(f"Коэффициент пьезопроводности: {kf:.2e} м^2/c")
        B = - (r_w**2) / (4 * kf * self.dt)
        print(f"Параметр B: {B:.6f}")
        Ei0 = -expi(B)
        print(f"Ei0: {Ei0:.6f}")
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

        return tV, Q

calculator = wellcalculator() 

pressure_profile = np.ones(calculator.N) * p_zab
pressure_profile[0:50] = p_inj
pressure_profile[50:] = p_zab

time_seconds, flow_rate = calculator.calculate_flow(h, fi, k, mu, c, p_0, r_w, pressure_profile)
time_days = time_seconds / 86400
flow_rate_m3_day = flow_rate * 86400
print ("Расчет завершен")

print("РЕЗУЛЬТАТЫ:")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
ax1.plot(time_days, np.ones(calculator.N) * p_0 / 1e5, '--b', label='Пластовое давление')
ax1.plot(time_days, pressure_profile / 1e5, '-r', label='Забойное давление')
ax1.set_ylabel('Давление, бар')
ax1.set_ylim([0,500])
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_title('Профиль давления в скважине')
    
# ax2.plot(time_days, flow_rate_m3_day, '-g')
ax2.plot(time_days[:50], flow_rate_m3_day[:50], '-g', label='До скачка')
ax2.plot(time_days[49:52], flow_rate_m3_day[49:52], '--g', linewidth=2, label='Скачок')
ax2.plot(time_days[51:], flow_rate_m3_day[51:], '-g', label='После скачка')
ax2.set_ylabel('Дебит, м³/сут')
ax2.set_xlabel('Время, дни')
ax2.grid(True, alpha=0.3)
ax2.set_title('Дебит скважины во времении')

ax1.set_xlim([0, time_days[-1]])
ax2.set_xlim([0, time_days[-1]])

plt.tight_layout()
plt.show()