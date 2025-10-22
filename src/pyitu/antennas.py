import numpy as np
from .tooling import c0

__desc__="""contains several ITU antenna definitions as functions with the same name code as the ITU document
def p465( freq , phi , D ) 
        - freq in GHz
        - phi in degrees
        - D in meters
        - bVerbose boolean prints more information if true
        - bLegacy uses model in use prior to 1993

        Recommendation ITU-R S.465-6
(01/2010)
Reference radiation pattern for earth station
antennas in the fixed-satellite service for
use in coordination and interference
assessment in the frequency range
from 2 to 31 GHz

NOTE 3 – For the purpose of determining the maximum permissible levels of interference in
Recommendations ITU-R S.466, ITU-R S.483, ITU-R S.523 and ITU-R S.735, receiving earthstation
antenna reference patterns no worse than stated in those Recommendations should apply


def f699( freq , phi , D , Gmax = None , GI = None, phi_r = None , phi_m = None , bF1245 = False )
        - freq in GHz
        - phi in degrees
        - D in meters
        | Gmax maximum gain
        | GI gain of the first sidelobe
        | phi_r sidelobe position (non gain dependent)
        | phi_m sidelobe position (gain dependent shift)
        | bF1245 boolean sets calculation to corresponding ITU-F.1245 pattern for comparison
 
Recommendation ITU-R F.699
(01/2018)
Reference radiation patterns for fixed
wireless system antennas for use in
coordination studies and interference
assessment in the frequency range
from 100 MHz to 86 GHz

phi_r is assumed to correspond to the off-axis angle of the peak of the first side-lobe and the phase at phi_r is assumed to be 1.5*pi
phi_m also corresponds to the off-axis angle of the first sidelobe


 def ap8( freq , phi , D , bVerbose = False )
        - freq in GHz
        - phi in degrees
        - D in meters
        | bVerbose boolean prints extra information
         
Appendix 8 Earth station antenna pattern for GSO networks. Only
for maximum antenna gain greater than 9.3 dB.

    Name \t: APERR_001V01
    Type \t: Earth station, Receiving and Transmitting
    AP8  \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APERR_001V01.pdf
"""

def help() :
    print(__desc__)

def ap8( freq , phi , D , bVerbose = False ) :
    #
    # AP8 : https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APERR_001V01.pdf
    #
    desc_ = """Appendix 8 Earth station antenna pattern for GSO networks. Only
for maximum antenna gain greater than 9.3 dB.

    Name \t: APERR_001V01
    Type \t: Earth station, Receiving and Transmitting
    AP8  \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APERR_001V01.pdf
    """
    if bVerbose :
        print ( desc_ )
    #
    # CONSTANTS
    # lambda is a reserved word in python so we use k for wavelength (not wavenumber)
    k     = c0/(freq*10**9)
    kD    = k / D
    Dk    = 1 / kD
    bCase = Dk >= 100
    phi   = 0 if phi<0 else 180 if phi>180 else phi

    if phi_b < phi_r : 
        print("Warning:: Phib () is less than Phir ().")

    Gmax  = 20 * np.log10(Dk) + 7.7
    G1    = 2 + 15 * np.log10(Dk)
    phi_m = 20 * kD * np.sqrt(Gmax-G1)
    phi_b = 48
    phi_r = 15.85*kD**0.6 if bCase else 100*kD

    constants = {True:tuple((32,0,-10)),False:tuple((52,-10,10))}[bCase]

    phicase = 1*int(phi < phi_m) + 2*int(phi >= phi_m and phi < phi_r) + 3*int(phi_r <= phi and phi < phi_b) + 4*int(phi_b<phi)

    match  phicase :
        case 1:
            G = Gmax - 2.5 * 10**(-3) * (Dk * phi)**2
        case 2:
            G = G1
        case 3:
            G = constants[0] + constants[1] * np.log10(Dk) - 25 * np.log10(phi)
        case 4:
            G = constants[2] + constants[1] * np.log10(Dk)
        case _ :
            print("ERROR: Bad case value")
    return (G)


def f699( freq , phi , D , Gmax = None , GI = None, phi_r = None , phi_m = None , bF1245 = False )
    k = c0/(freq*10**9)

    if GI is None :
        # Implements ITU-R F.1245 / F.699 Note
        GI = 2 + 15*np.log10(D/k)
    if Gmax is None :
        # Implements ITU-R F.1245 / F.699 NOTE
        Gmax = 20.0*np.log10(D/k) + 7.7

    if phi_r is None or phi_m is None :
        # Implements ITU-R F.1245 / F.699
        factor = 12.02 if bF1245 else 15.85
        phi_r = factor * (k/D)**0.6
        phi_m = 20*(k/D) * np.sqrt(Gmax-GI)

    bLow = freq>=1 and freq<=70

    if D/k > 100 :
        phic = {True:48,False:120}
        constants = { tuple((True,True))  : [29,-13] , tuple((False,True)) : [29,-23] ,
                      tuple((True,False)) : [32,-10] , tuple((False,False)): [32,-20] }

        if phi>0 and phi<phi_m :
            G = Gmax - 2.5*10**(-3) * (D/k * phi)**2
        if phi >= phi_m and phi < np.max([phi_m,phi_r]) :
            G = GI
        if np.max([phi_m,phi_r]) <= phi  and phi < phic[bLow] :
            G = constants[tuple((bLow,bF1245))][0] - 25 * np.log10(phi)
        if phi >= phic[bLow] and phi <= 180 :
            G = constants[tuple((bLow,bF1245))][1]

    if D/k<=100 :
        phic = {True:48,False:120}
        constants = { tuple((True,True))  : [39,-3,0.5] , tuple((False,True)) : [39,-13,0.5] ,
                      tuple((True,False)) : [52,-10,1.0] , tuple((False,False)): [52,0,1.0] }

        phi_L = 100*k/D
        if phi>0 and phi<phi_m :
            G = Gmax - 2.5*10**(-3) * (D/k * phi)**2
        if phi >= phi_m and phi < phi_L :
            G = GI
        if phi_L <= phi  and phi < phic[bLow] :
            G = constants[tuple((bLow,bF1245))][0] - 25 * np.log10(phi) - 10 * np.log10(D/k) * constants[tuple((bLow,bF1245))][2]
        if phi >= phic[bLow] and phi <= 180 :
            G = constants[tuple((bLow,bF1245))][1] - 10 * np.log10(D/k) * constants[tuple((bLow,bF1245))][2]
    return G

def s465(freq , phi , D ,
          bVerbose = False , bLegacy = False) :
    return ( p465( freq , phi , D , bVerbose , bLegacy ) )

def p465( freq , phi , D ,
          bVerbose = False , bLegacy = False) :

    k = c0/(freq*10**9)
    if bLegacy :
        phimin = 100*k/D
        G = 52 - 10*np.log10(D/k) - 25*np.log10(phi) if phi >=phimin and phi<48 else 10 - 10 *np.log10(D/k)
        return G

    if D/k>=50 :
        phimin = np.max([1,100*k/D])
    else :
        phimin = np.max([2,114*(k/D)**1.09])
    if bVerbose :
        print ( phimin , D/k )
    if freq >= 2 and freq < 31 :
        if phi >= phimin and freq < 48 :
            G = 32 - 25*np.log10(phi) # ϕmin ≤ ϕ < 48°
        if phi > 48 and phi<=180 :
            G = -10
    return G

if __name__=='__main__':

    print("ANTENNAS")
    help()
    print ( p465( freq=2 , phi=15 , D=1.0 , bVerbose=True ) )
