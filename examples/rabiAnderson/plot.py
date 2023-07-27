"""
#!/usr/bin/env python3
DOC STRING
"""

# -*- coding: utf-8 -*-
import h5py
import matplotlib.pyplot as plt
# plt.figure(figsize=[10, 4.6])
import numpy as np
gamma = 0.01
cm = plt.get_cmap('jet')
# from matplotlib.backends.backend_pdf import PdfPages
nr_num = 20  # half of the nr iterations
nr_energy = 100
# ----------------------------------------------------------------
xr = np.full((nr_num*2, 21000), np.inf)
yr = np.full((nr_num*2, 21000), np.inf)
ifl = "results-test.h5"
f = h5py.File(ifl, "r")
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
# ----------------------------------------------------------------
xr = np.full((nr_num*2, 21000), np.inf)
yr = np.full((nr_num*2, 21000), np.inf)
ifl = "results.h5"
f = h5py.File(ifl, "r")
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
    ax1.plot(2*y, xr[enum][0:nr_num], "-", color='orange', alpha=0.5)
    ax2.plot(2*y+1, yr[enum][0:nr_num], "-", color='orange', alpha=0.5)


plt.savefig("out.pdf")
