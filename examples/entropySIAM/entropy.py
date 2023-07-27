"""
#!/usr/bin/env python3
DOC STRING
"""

# -*- coding: utf-8 -*-
import h5py
import numpy as np
import matplotlib.pyplot as plt

# import sys
# plt.figure(figsize=[10, 4.6])


ifl = "resultSIAMZero.h5"
f = h5py.File(ifl, "r")
x20 = np.array(f["entropy"])
kbt = np.array(f["temperatureArray"])


ifl = "resultSIAMU0.200000.h5"
f = h5py.File(ifl, "r")
x2 = np.array(f["entropy"])


# plt.plot(kbt, x1 - x10, "o", label="Kondo")
plt.plot(kbt, x2 - x20, "o", label="SIAM")
plt.plot(kbt, x2 * 0 + np.log(2), "--", label="log(2)")
plt.plot(kbt, x2 * 0 + np.log(4), "--", label="log(4)")
# set x-axis as log
plt.xscale("log")
plt.xlim(1e-8, 0.250)
plt.legend()
#
plt.savefig("out.pdf")
