from coeio import *
#from coeplottk import *
#from cosmocoe import *
#from scipy.interpolate import interp2d

field = 'omegaCen'
#field = 'simulated'

RGB = loadrgb('starrybackground-%s.jpg' % field)
three, ny, nx = RGB.shape

trim = (nx - ny) / 2
RGB = RGB[:, :, trim:-trim]
three, ny, nx = RGB.shape

yc = ny / 2
xc = nx / 2

xo = 0
yo = 0

yy, xx = mgrid[0:ny,0:nx]
dx = xx - xc - xo
dy = yy - yc - yo
rr = hypot(dx, dy)

M = 2e5
#M = M * (1/2.)**2
M = M * (2/3.)**2
#M = 1e5
# RE ~ sqrt(M)

ax = M * dx / rr**2
ay = M * dy / rr**2

ys = yy - ay
xs = xx - ax

ys = ys.astype(int)
xs = xs.astype(int)

edge = 800
ys2 = (abs(ys) - ny/2) % edge
xs2 = (abs(xs) - nx/2) % edge

goody = between(0, ys, ny-1)
goodx = between(0, xs, nx-1)
#good = goodx * goody

ys = where(goody, ys, ys2)
xs = where(goodx, xs, xs2)

#ys = ys % nyc  # put image on loop
#xs = xs % nxc  # put image on loop

lensed = RGB[:,ys,xs]

savergb(lensed, 'lensed-%s.jpg' % field)
