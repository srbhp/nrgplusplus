"""
#!/usr/bin/env python3
DOC STRING
"""

# -*- coding: utf-8 -*-
import h5py
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# plt.figure(figsize=[10, 4.6])
import numpy as np
gamma = 0.01
cm = plt.get_cmap('jet')
pp = PdfPages("out.pdf")
nr_num = 30  # half of the nr iterations
nr_energy = 100
noStates = 85000
# ----------------------------------------------------------------
ifl = "resultsTc.h5"
f = h5py.File(ifl, "r")
nrgMax = 39
nrgMin = 4
x = np.array(f["GreenFnEnergyPoints"])
Emin = -10
Emax = 1
xr = np.logspace(Emin, Emax, 100)
bpar = 0.7
prefac = 1.0 / (bpar * np.sqrt(np.pi))
for ii in range(2):
    z = np.array(
        f["GreenFnNegativeWeight{}".format(ii)])
    y = np.array(
        f["GreenFnPositiveWeight{}".format(ii)])
    yr = 0 * xr
    zr = 0 * xr
    for i in range(len(xr)):
        yr[i] = (1 / xr[i]) * np.sum(
            y * np.exp(-np.power(-bpar * 0.25 + (np.log(xr[i] / x) / bpar), 2))
        )
        zr[i] = (1 / xr[i]) * np.sum(
            z * np.exp(-np.power(-bpar * 0.25 + (np.log(xr[i] / x) / bpar), 2))
        )
    yr = yr * prefac * gamma * np.pi
    zr = zr * prefac * gamma * np.pi
    #
    yr = yr/np.sum(yr)
    zr = zr/np.sum(zr)
    plt.plot(xr, yr,  label="Positive : {}".format(ii))
    plt.plot(xr, zr,  label="Negative: {}".format(ii))
    # plt.plot(rx[:-1], ry, label="Positive : {}".format(i))
    # plt.plot(rx[:-1], ry, label="Negative : {}".format(i))
    print("Sum : ", np.sum(y))
    print("Sum : ", np.sum(z))
# plt.plot(xr, yr*0 + 1, '--')
#
#
# plt.xlim(-15, 0)
plt.xscale('log')
plt.legend()
pp.savefig()
plt.clf()
# ----------------------------------------------------------------
xr = np.full((nr_num*2, noStates), np.inf)
yr = np.full((nr_num*2, noStates), np.inf)

for ic in np.arange(nr_num):
    try:
        x = np.array(f["Eigenvalues{}".format(ic*2)])
        y = np.array(f["Eigenvalues{}".format(ic*2+1)])
        print("Eigenvalues{}:{}".format(2*ic + 1, x.size))
        xr[ic][0:x.size] = x - x[0]
        yr[ic][0:y.size] = y - y[0]
    except Exception as e:
        print("Error: {}".format(e))
        pass
xr = xr.T
yr = yr.T
y = np.arange(nr_num)
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
for enum in np.arange(nr_energy):
    ax1.plot(2*y, xr[enum][0:nr_num], "--", color='red', lw=0.25)
    ax2.plot(2*y+1, yr[enum][0:nr_num], "--", color='red', lw=0.25)


pp.savefig()
pp.close()
