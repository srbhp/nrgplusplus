#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np 
import matplotlib.pyplot as pl 
from matplotlib.backends.backend_pdf import PdfPages
pp=PdfPages("out.pdf") 

flist="spectr-U0.001.dat  spectr-U0.01.dat  spectr-U0.1.dat"

flist=flist.split()
for ifl in flist:
    x,y,z=np.loadtxt(ifl,unpack=True)
    xx,yy=np.histogram(x,weights=y,bins=100)
    print len(yy),len(xx)
    #pl.xscale("log")
    #pl.yscale("log")
    pl.plot((yy[:-1]+yy[1:])/2.0,xx,label='U='+ifl[8:13])

#pl.yscale("log")
#pl.xlim(-0.01,0.01)
pl.ylabel(r"DoS")
pl.xlabel(r"$\omega$")
pl.legend(loc="best")
pp.savefig()
pl.clf()
#fig = pl.gcf()
#fig.set_size_inches(10.0,4.0)

pp.close()

