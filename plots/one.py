#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np 
import matplotlib.pyplot as pl 
from matplotlib.backends.backend_pdf import PdfPages
pp=PdfPages("out.pdf") 

flist=" "
flist=flist.split()
x,y,z=np.loadtxt("spectr-U0.001.dat",unpack=True)
xx,yy=np.histogram(x,weights=y,bins=100)
#xx,yy=np.histogram(abs(x),weights=y+z,bins=np.logspace(-5,-1,1000))
print len(yy),len(xx)
#pl.xscale("log")
pl.plot((yy[:-1]+yy[1:])/2.0,abs(xx),label='Up')
#pl.xscale("log")
"""
xx,yy=np.histogram(x,weights=z,bins= 100 )
print len(yy),len(xx)
#pl.xscale("log")
#pl.yscale("log")
pl.plot((yy[:-1]+yy[1:])/2.0,xx,label='Down')
"""

#pl.xlim(-0.5,0.5)
pl.ylabel(r"DoS")
pl.xlabel(r"$\omega$")
pl.legend(loc="best")
pp.savefig()
pl.clf()
#fig = pl.gcf()
#fig.set_size_inches(10.0,4.0)

pp.close()

