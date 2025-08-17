from coeplottk import *
from cosmocoe import *

zl = 0.584  # MACS0647

def loadmodel(i=None, fac=1):
    global ax, ay, nx, ny, X, Y, XX, YY, model
    
    ax = loadfits('hlsp_relics_model_model_abell2537_lenstool_v1_x-pixels-deflect')
    ay = loadfits('hlsp_relics_model_model_abell2537_lenstool_v1_y-pixels-deflect')
            
    ny, nx = ax.shape
    X = arange(nx)
    Y = arange(ny)
    XX, YY = meshgrid(X, Y)

# ~/A2261/SL/Adi/4/color/imAdi.py
#xo, yo = 1437, 1123
xo, yo = 1756, 1750
imdx = -xo  # 1437
imdy = -yo  # 1878#
#xfac = yfac = 1
xfac = yfac = 0.065 / 0.05
#imdx = imdy = 0

#ny = nx = 2496
#yo = ny - yo

# kappa BCG:    1250.5, 1251.5
# mosdriz BCG: ~2687,   2373
# drex BCG:    ~3150,   2907
# difference is xo, yo

#imdx = -2148
#imdy = -1896

# drex -> mosdriz
#imdx = imdx - 462
#imdy = imdy - 543

def addcrosshair(R, G, B, x, y, color=(255,255,255), d=15, d2=2):
    #print x, y, c
    R[y-d : y+d+1, x-d2 : x+d2+1] = color[0]
    G[y-d : y+d+1, x-d2 : x+d2+1] = color[1]
    B[y-d : y+d+1, x-d2 : x+d2+1] = color[2]
    
    if d <> d2:
        R[y-d2 : y+d2+1, x-d : x+d+1] = color[0]
        G[y-d2 : y+d2+1, x-d : x+d+1] = color[1]
        B[y-d2 : y+d2+1, x-d : x+d+1] = color[2]
    
    return R, G, B

#def relenscolor((x, y, zs)):
def relenscolor(x, y, zs, xc=[], yc=[], sh=1):  # imodel=4, 
    #loadmodel(imodel)

    x = xfac * x + imdx
    y = yfac * y + imdy
    if len(xc):
        xc = array(xc) + imdx
        yc = array(yc) + imdy
    
    ax1 = ax[y,x]
    ay1 = ay[y,x]
    
    Dds_Ds1 = Dds_Ds(zl, zs)
    
    xs1 = x - Dds_Ds1 * ax1
    ys1 = y - Dds_Ds1 * ay1
    
    xs0 = XX - Dds_Ds1 * ax
    ys0 = YY - Dds_Ds1 * ay
    
    dist = hypot(xs0-xs1, ys0-ys1)
    dd = clip(dist, 0, dmax)
    dd = (dmax - dd) / dmax * 255
    #showarr(dd)
    dd = dd.astype(int)

    #im1, RGBin = loadimage(inimage)
    im1, RGBin = loadimage()
    #R, G, B = RGBin
    R, G, B = RGBin[:,1:-2,1:-2]
    #R, G, B = RGBin[:]
    #R, G, B = RGBin = loadrgb(inimage)

    #print G.shape
    print dd.shape
    
    if acsir:
        dd = where(greater(dd, 70), dd, 0)

        G = where(dd, 0, G)
        R = where(dd, dd,  R)
        B = where(dd, 0,  B)

        # green blobs:
        #G = where(dd, dd, G)
        #R = where(dd, 0,  R)
        #B = where(dd, 0,  B)
        
        #G = where(dd, dd, G)
        #G = where(greater(dd, 50), dd, G)
        #G = max((G, dd))  # BRIGHTEN
    else:
        # green overwrite
        #dd = where(greater(dd, 70), dd, 0)
        #G = where(dd, dd, G)
        #R = where(dd, 0,  R)
        #B = where(dd, 0,  B)

        # green:
        B = min((B, 255-dd))
        R = min((R, 255-dd))
        # magenta:
        # G = min((G, 255-dd))
        # cyan:
        # R = min((R, 255-dd))
        # yellow:
        # B = min((B, 255-dd))

    for i in range(len(xc)):
        R, G, B = addcrosshair(R, G, B, xc[i], yc[i], color=(255,70,255))
        #R, G, B = addcrosshair(R, G, B, xc[i], yc[i], color=(0,0,0), d2=1)
    
    if acsir:
        crosshaircolor = (255,255,255)
    else:
        crosshaircolor = (0,0,0)
    
    R, G, B = addcrosshair(R, G, B, x, y, color=crosshaircolor)
    
    RGB = R, G, B
    
    im = rgb2im(RGB)
    if sh:
        im.show()
    return im

acsir = 1

#ny, nx = ax.shape

#dmax = 10.
#dmax = 30.
dmax = 12.

#def loadimage(inimage):
def loadimage():
    inimage = 'imAdilab.png'
    im = Image.open(inimage)
    RGBin = loadrgb(inimage)
    return im, RGBin



#################################

# grab coords from glens/arcs.py

# im = relenscolor((x, y, zs))

#dmax = 100

#x, y, zs = 3532, 2760, 1.13  #  1.1  big arc: drex coords
#x, y, zs = 3071, 2221, 1.13  #  1.1  big arc: drex coords

#x, y, zs = 3029, 2797, 4.0  #  claw.1
#x, y, zs = 2669, 2144, 4.0  #  claw.1
#x, y, zs = 2901, 2701, 4.0  # ??

#x, y, zs = 3029, 2797, 1.0  #  claw.1

#im = relenscolor((x, y, zs))
