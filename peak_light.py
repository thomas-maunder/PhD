import numpy as np

with open("lightcurve.txt") as f:
    n = len(list(f))

lightcurve = np.loadtxt('lightcurve.txt', skiprows=2, max_rows=n-4)

print(np.round(lightcurve[np.argmax(lightcurve[:,1]),0]), np.max(lightcurve[:,1]))
