import numpy as np
import matplotlib.pyplot as plt
import math

def F(x):
  return 1.0 + np.sin(0.8*x*math.pi)

#====================================== Generate Fx
L = 5.0
x = np.linspace(0.0,L,100)
F_x = F(x)

#====================================== Monte Carlo integration
Is = 0.0
Np = 10000
for i in range(0,Np):
  xr = np.random.uniform()*L
  Is += F(xr)
Is = L*Is/Np

print("Monte Carlo integration:",Is)

#====================================== Numerical integration (Riemann)
In = 0.0
dx = L/Np
for i in range(0,Np):
  xn = 0.5*dx + dx*i
  In += F(xn)*dx

print("Numerical integration:",In)

# ===================================== Plotting
plt.plot(x,F_x)
plt.show()