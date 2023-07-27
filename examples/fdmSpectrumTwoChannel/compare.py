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
nrgMax = 39
nrgMin = 4
Emin = -9.
Emax = 1
xr = np.logspace(Emin, Emax, 100)
bpar = 0.7
prefac = 1.0 / (bpar * np.sqrt(np.pi))
flist = "resultsSIAM.h5 resultsTcRT0.000000.h5"
flist = flist.split()
labels = ["SIAM", "Two Channel"]
for ifl in flist:
    f = h5py.File(ifl, "r")
    ii = 0
    x = np.array(f["GreenFnEnergyPoints"])
    z = np.array(
        f["GreenFnNegativeWeight{}".format(ii)])
    y = np.array(
        f["GreenFnPositiveWeight{}".format(ii)])
    print("Sum : ", np.sum(y))
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
    plt.plot(xr, yr,  label="Positive : {}".format(ifl))
    # plt.plot(xr, zr,  label="Negative: {}".format(ii))
    # plt.plot(rx[:-1], ry, label="Positive : {}".format(i))
    # plt.plot(rx[:-1], ry, label="Negative : {}".format(i))
# plt.plot(xr, yr*0 + 1, '--')
#
#
# plt.xlim(-15, 0)
plt.xscale('log')
plt.legend(labels)
pp.savefig()
#
plt.yscale('log')
plt.ylim(1e-3,)
pp.savefig()
plt.clf()
# ----------------------------------------------------------------
pp.close()
