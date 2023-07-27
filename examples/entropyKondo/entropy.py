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

ifl = "resultsZero.h5"
f = h5py.File(ifl, "r")
x10 = np.array(f["entropy"])

ifl = "resultsJ0.250000.h5"
f = h5py.File(ifl, "r")
x1 = np.array(f["entropy"])
kbt = np.array(f["temperatureArray"])

plt.plot(kbt, x1 - x10, "o", label="Kondo")
plt.plot(kbt, x1 * 0 + np.log(2), "--", label="log(2)")
# plt.plot(kbt, x1 * 0 + np.log(4), "--", label="log(4)")
# set x-axis as log
plt.xscale("log")
plt.xlim(1e-6, 0.10)
plt.ylabel("Entropy")
plt.xlabel("Temperature")
# plt.legend()
#
plt.savefig("out.pdf")
