# ~/LensPerfect/A1689/2009/DMprofiles/plots/NFWfit.py
# ~/LensPerfect/A1689/analysis/NFWfit.py
# ~/LensPerfect/A1689/analysis/NFWconst3.py
# ~/glens/A1689/arcs/take4/
# NFWconst3.py
# NFWconst2.py
# profile3.py

from units import *

pi = 3.14159

from coeio import *
from cosmocoe import DA, H, Dds_Ds
from scipy.optimize import golden

Om = 0.27
OL = 1 - Om
def dccalc(z=None):  # Virial overdensity (instead of 200)
    if z == None:
        z = zL
    
    Omz = 1 / (1 + (OL/Om)*(1+z)**-3)
    OLz = 1 - Omz
    dc = 18 * pi**2 - 82*OLz - 39*OLz**2
    return dc

def ddc(c, dc, dvir):
    #print c, dc, dvir
    dc1 = (dvir / 3.) * c**3 / (log(1.+c) - c/(1.+c))
    return abs(dc1 - dc)

def cconv(c, dvir, dvirnew):
    if singlevalue(c):
	dc = (dvir / 3.) * c**3 / (log(1.+c) - c/(1.+c))
	cnew = golden(ddc, args=(dc, dvirnew), brack=[1, 30])
    else:
	cnew = []
	for c1 in c:
	    cnew.append(cconv(c1, dvir, dvirnew))
	cnew = array(cnew)
    return cnew

def dMvirrs(rsKPC, c, dvir, Mvir):
    NFW1 = NFW(zL, H0/H100)  # zL,h added 4/20/10!  Otherwise defaulted to A1689, h=0.7
    NFW1.set(rsKPC * KPC, c = c, dvir = dvir)
    Mvir1 = NFW1.Msph(NFW1.rvir).asNumber(MSUN)
    return abs(Mvir.asNumber(MSUN) - Mvir1)

class NFW:
    def __init__(self, z='A1689', h=0.7):
	global zL, DL, Ec, H0, pc0, pc, p_E
	#zL = loaddict('~/lp/zlens.dict', silent=1)[clusname]
	zL = loaddict('~/lp/zlens.dict', silent=1).get(z, z)
	H0 = h * H100
	pc0 = (3 * H0**2) / (8 * pi * G)
	#pc = pc0 * (1 + zL) ** 3  # 227 Msun / Kpc**3
	pc  = (3 * H(zL)**2) / (8 * pi * G)  # ERROR CORRECTED IN profile5.py

        if z > 0:
            DL = DA(0, zL) # * MPC
            Ec = c**2 / (4 * pi * G * DL)
            p_E = pc / Ec
    
    def test(self):
	print 'DL', DL.asunit(MPC)
	# 641.516168049 [Mpc]
	print 'c', c
	print 'pi', pi
	print 'G', G
	print 'DL', DL
	print 'Ec', Ec.asunit(KG / M**2)
	# 5.3848648036 [kg/m2]
	print 'H0', H0.asunit(KM / S / MPC)
	# 2.26854372456e-18 [1/s]
	print 'pc0', pc0.asunit(MSUN / KPC**3)
	# 135.993372864 [Msun/Kpc3]
	print 'p_E', p_E.asunit(1/KPC)
	# 8.80291134071e-08 [1/Kpc]
	print 'pc', pc.asunit(MSUN / KPC**3)
	# 226.982044872 [Msun/Kpc3]
	
    # ################################
    # NFW FIT
    # ~/glens/bestfit/LensPerfect/Eddie/NFWconst.py
    # SEE ALSO NFWfit3.py
    
    def set(self, rs=0*MPC, c=None, dvir=200, rvir=0*MPC, dc=None, ks=None, ps=0*KG/M**3,
	    Mvir=0*MSUN):

	if Mvir > 0*MSUN and c and dvir:
	    rs = golden(dMvirrs, args=(c, dvir, Mvir), brack=[1, 10000]) * KPC
	    #self.Mvir = Mvir
	
	if ks:
	    if rs == 0*MPC:
		rs = asNumber(ks / p_E / dc)
	    elif dc == None:
		dc = asNumber(ks / p_E / rs)
	
	if ps > 0*KG/M**3:
	    dc = asNumber(ps / pc)

	if c == None:
	    if rs > 0*MPC and rvir > 0*MPC:
		c = rvir / rs
	    elif dc:
		c = golden(ddc, args=(dc, dvir), brack=[1, 30])
	
	if dc == None:
	    dc = (dvir / 3.) * c**3 / (log(1.+c) - c/(1.+c))
	
	if rs == 0*MPC:
	    if rvir > 0*MPC and c:
		rs = rvir / c
	if rvir == 0*MPC:
	    rvir = rs * c
	
	if ks == None:
	    ks = asNumber(p_E * dc * rs)

	if ps == 0*KG/M**3:
	    ps = dc * pc

	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('self.%s = %s' % (label, label))
	
    def p(self, R):
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))
	
	"""NFW 3-D profile GIVEN R, rs in kpc"""
	x  = asNumber(R / rs)
	p = ps / x / (1 + x)**2
	return p

    def F(self, x):
	Flo = (1 - arccosh(1/x) / sqrt(1 - x**2)) / (x**2 - 1)
	F1 = 1 / 3.
	Fhi = (1 - arccos(1/x)  / sqrt(x**2 - 1)) / (x**2 - 1)
	Fall = where(greater(x, 1), Fhi, Flo)
	Fall = where(equal(x, 1), F1, Fall)
	return Fall
    
    def k(self, R, zs=None):
	"""NFW profile GIVEN R, rs in kpc"""
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))
	
	x  = asNumber(R / rs)
	k  = 2 * ks * self.F(x)
        if zs <> None:
            k = k * Dds_Ds(zL, zs)
        
	return k
    
    def y(self, R, zs=None):  # SHEAR
	"""NFW shear GIVEN R, rs in kpc"""
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))
	
	x  = asNumber(R / rs)
	y  = 2 * ks * (2 * self.Gcyl(x) / x**2 - self.F(x))
        if zs <> None:
            y = y * Dds_Ds(zL, zs)
        
	return y

    def g(self, R, zs=None):  # REDUCED SHEAR!!
	"""NFW reduced shear GIVEN R, rs in kpc"""
	g = self.y(R, zs) / (1 - self.k(R, zs))
	return g

    # ~/glens/NFW/NFW.py
    # ~/glens/A1689/arcs/take4/profile5.py
    def Gcyl(self, x):
	Glo = log(x / 2.) + arccosh(1/x) / sqrt(1 - x**2)
	G1  = log(x / 2.) + 1
	Ghi = log(x / 2.) + arccos(1/x)  / sqrt(x**2 - 1)
	Gall = where(greater(x, 1), Ghi, Glo)
	Gall = where(equal(x, 1), G1, Gall)
	return Gall

    def Mcyl(self, R):
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))

	"""MASS INSIDE RADIUS"""
	x  = asNumber(R / rs)
	return 4 * pi * ks * self.Gcyl(x) * Ec * rs**2

    def kavgcyl(self, R, zs=1e30):
	"""AVERAGE kappa INSIDE CYLINDER"""
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))

	# Mcyl / (pi * R**2 * Ec)
	x  = asNumber(R / rs)
	Dds_Ds1 = Dds_Ds(zL, zs)
	return 4 * ks * self.Gcyl(x) / x**2 * Dds_Ds1

    def RE(self, xunit=KPC, zs=1e30):
	"""EINSTEIN RADIUS: THAT WHICH CONTAINS kavgcyl = 1"""
	#RE = golden(lambda Rx: abs(self.kavgcyl(Rx*xunit) - 1),brack=[30, 70])
	Dds_Ds1 = Dds_Ds(zL, zs)
	def dk(Rx):
	    x  = asNumber(Rx * xunit / self.rs)
	    kavg = 4 * self.ks * self.Gcyl(x) / x**2 * Dds_Ds1
	    return abs(kavg - 1) 
	    # SLOW BECAUSE CALCULATES Dds_Ds
	    #return abs(self.kavgcyl(Rx*xunit, zs=zs) - 1) 

	RE = golden(dk, brack=[30, 70]) * xunit
	#print RE * xunit
	#print RE * xunit.asunit(KPC)
	#print self.kavgcyl(RE*xunit, zs=zs)
	return RE

    def Gsph(self, x):
	return log(1 + x) - x / (1 + x)

    def Msph(self, R):
	for label in string.split('rs c dvir rvir dc ks ps'):
	    exec('%s = self.%s' % (label, label))

	"""MASS INSIDE RADIUS"""
	x  = asNumber(R / rs)
	return 4 * pi * ps * rs**3 * self.Gsph(x)

    def Mvir(self):
        return self.Msph(self.rvir).asunit(MSUN)
