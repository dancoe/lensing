# ~/p/shearconv.py -- shear2xy, shear2etheta
# ~/p/shear2xy.py -- NOW INPUTS ARRAYS
# PREVIOUSLY conv2()

# FROM ~/p/whiskerplot.py
# FROM /im2shape/test/2c/testetheta/conv.py


# a, b = Major, Minor axes of ellipse

# g = (a - b) / (a + b)
#  = (1 - b / a) / (1 + b / a)
#  = e / (2 - e)
# (a / b) = (1 + g) / (1 - g)

# e = 1 - (b / a) = 2 * g / (g + 1)


# FORMULATED IN THIS WAY, g SHOULD MATCH THE REDUCED SHEAR:
# g = gamma / (1 - kappa)
# (or 1 / g* if |g| > 1)
# where gamma1 = 1/2 (phi_11 - phi_22) ; gamma2 = phi_12 = phi_21
# and   kappa  = 1/2 (phi_11 + phi_22)


# g1 + 0 - 0
# g2 0 + 0 -
#    - / | \  (KSB, me usually & above)
#    | \ - /  (Bridle+, mkstarssheared.py)

# conversion: Bridle theta = KSB theta + 90
#             g1 & g2 reverse their signs

# g = sqrt(g1^2 + g2^2)
# g1 = g cos(2*theta)
# g2 = g sin(2*theta)

from MLab_coe import *

def atanxy(x, y, degrees=0):
    """ANGLE CCW FROM x-axis"""
    theta = arctan(divsafe(y, x, inf=1e30, nan=0))
    theta = where(less(x, 0), theta + pi, theta)
    theta = where(logical_and(greater(x, 0), less(y, 0)), theta + 2*pi, theta)
    if degrees:
        theta = theta * 180. / pi
    return theta

def shear2(g, theta, vert=0):
    theta = theta / 180. * pi
    g1 = g * cos(2*theta)
    g2 = g * sin(2*theta)
    if vert:
        g1 = -g1
        g2 = -g2
    return g1, g2

def shear2xy(g1, g2, vert=0):
    """(g1, g2) -> (e, theta) -> (ex, ey)"""
    theta = arctan(g2 / g1) / 2.  # CCW from x-axis
    #
    g1pos = greater(g1,0)
    g2pos = greater(g2,0)
    #
    thetanegg1 = where(g2pos, theta+pi/2., theta-pi/2.)
    theta = where(g1pos, theta, thetanegg1)
    theta = where(g1, theta, where(g2pos, pi/4., -pi/4.))  # / \
    
    # g1 + 0 - 0
    # g2 0 + 0 -
    #    - / | \  (KSB, me usually & above)
    #    | \ - /  (Bridle+, mkstarssheared.py)
    # conversion: Bridle theta = KSB theta + 90
    if vert:
        theta += pi/2.  # | \ - /

    g = hypot(g1, g2)
    e = 2 * g / (g + 1)
    ex = e * cos(theta)
    ey = e * sin(theta)

    return ex, ey

def shear2gtheta(g1, g2, vert=0):
    """(g1, g2) -> (g, theta)"""
    # g1 + 0 - 0
    # g2 0 + 0 -
    #    - / | \  (KSB, me usually & above)
    #    | \ - /  (Bridle+, mkstarssheared.py)
    # conversion: Bridle theta = KSB theta + 90
    theta = atanxy(g1, g2) / 2.    
    if vert:
        theta += pi/2.  # | \ - /
    
    g = hypot(g1, g2)
    
    return g, theta

def shear2etheta(g1, g2, vert=0):
    """(g1, g2) -> (e, theta) -> (ex, ey)"""
    g, theta = shear2gtheta(g1, g2, vert)
    e = 2 * g / (g + 1)
    return e, theta
    
def exy2shear(ex, ey, vert=0):
    """(ex, ey) -> (e, theta) -> (g1, g2)"""
    e = hypot(ex, ey)
    g = e / (2 - e)
    ang = atanxy(ex, ey)
    if vert:
        ang -= pi/2.  # | \ - /
        
    g1 = g * cos(2*ang)
    g2 = g * sin(2*ang)
    return g1, g2


class NoClass:
    pass

# ~/im2shape/test/MS1358/backpsfcor/etprofile4.py
def radial(x, y, g1, g2, xc=0, yc=0):
    TINY = 1.e-8
    dx = x - xc
    dy = y - yc
    dd = dx**2 + dy**2 + TINY
    cos2 = (dx ** 2 - dy ** 2) / dd
    sin2 = 2 * dx * dy / dd
    e0, e1 = g1, g2
    out = NoClass()
    out.r = sqrt(dd)
    out.et = - (cos2 * e0 + sin2 * e1)
    out.er =    sin2 * e0 - cos2 * e1
    out.ee =   (sin2 * e0 + cos2 * e1) ** 2
    return out


def printconvention():
    print r'MUST SPECIFY SHEAR CONVENTION:'
    print
    print r'"h" (KSB)      "v" (Bridle)'
    print r'----------     ------------'
    print r'  g2              g2'
    print
    print r'   /               \ '
    print
    print
    print r'|  +  - g1      -  +  | g1'
    print
    print
    print r'   \               / '
    print
    #print 'g1 + 0 - 0'
    #print 'g2 0 + 0 -'
    #print '   - / | \  ("h": KSB)'
    #print '   | \ - /  ("v": Bridle)'


# shear2gtheta:
    #g1pos = greater(g1,0)
    #g2pos = greater(g2,0)
    #    
    #theta = arctan(divsafe(g2, g1, inf=1e30)) / 2.  # CCW from x-axis
    #thetanegg1 = where(g2pos, theta+pi/2., theta-pi/2.)
    #theta = where(g1pos, theta, thetanegg1)
    #theta = where(g1, theta, where(g2pos, pi/4., -pi/4.))  # / \
