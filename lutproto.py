import numpy as np
a0=1.42
delta = 0.00
nu = 0.85

ax = (np.sqrt(3.0)/8.0)*a0*(4.0+delta-(3.0*delta*nu))
ay = (3.0/8.0)*a0*(4.0+(3.0*delta)-(delta*nu))

stretch = 3
xb = stretch*ax
yb = stretch*ay
zb = 20.0
resolution = 3
zres = 3
x = np.linspace(-xb,xb,2*resolution*stretch+1)
x2 = np.linspace(-xb,xb,2*resolution*stretch+1)
y = np.linspace(-yb,yb,2*resolution*stretch+1)
z = np.linspace(-zb,zb,2*zres+1)

xgrid = np.linspace(0.0,ax,resolution+1)
ygrid = np.linspace(0.0,ay,resolution+1)
zgrid = np.linspace(0.0,zb,zres+1)
x[np.where(x < 0.0)] = -x[np.where(x < 0.0)]
x[np.where(x > 2.0*ax)] = x[np.where(x > 2.0*ax)] - (2.0*ax*(np.floor(x[np.where(x > 2.0*ax)]/(2.0*ax))))
x[np.where(x > ax)] = 2.0*ax - x[np.where(x > ax)]
#for i in range(len(x)):
#	print np.where(xgrid == x[i])
#print ax
#print x
#print xgrid
#print len(xgrid)
#print x[0]!=xgrid[-1]

xind = np.arange(len(x))
if stretch%2:
	xind = np.abs(xind%((resolution)*2)-(resolution))
else:
	xind = np.abs(np.abs(xind%((resolution)*2)-(resolution))-resolution)
#print xind

#print xgrid[xind]
#print xgrid
#print ygrid
#print zgrid
print x2
print y
print z
