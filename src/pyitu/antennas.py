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

Mode        Pattern id
Receiving   	33
Transmitting	76


def ap30B( freq , phi , D ,
	  bVerbose = False , antenna_efficiency = 0.7 ) :
    #
    desc_ = Appendix 30B reference Earth station antenna pattern.
Recommendation ITU-R S.580-6 reference Earth station antenna pattern.

    Name \t: APEREC015V01
    Type \t: Earth station, Receiving and Transmitting
    S580 \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APEREC015V01.pdf

Appendix 30B Earth station antenna pattern since WRC-03 applicable for D/lambda > 50.
Pattern is extended for D/lambda < 50 as in Appendix 8.
Pattern is extended for angles greater than 20 degrees as in Recommendation ITU-R S.465-5.
Pattern is extended in the main-lobe range as in Appendix 7 to produce continuous curves.
BR software sets antenna efficiency to 0.7 for technical examination.

Receiving    \t 605
Transmitting \t	606


def ap30Baeq29( freq , phi , D , coefA = 29 ,
		bVerbose = False , antenna_efficiency = 0.7 ) :
		- freq in GHz
        - phi in degrees
        - D in meters
        | bVerbose boolean prints extra information
        | antenna_efficiency = 0.7
    
    Appendix 30B reference Earth station pattern with the improved side-lobe for coefficient A = 29.

    Name \t: APERR_002V01
    Type \t: Earth station, Receiving and Transmitting
    AP30B A = 29  \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APERR_002V01.pdf

Receiving   	30, 31
Transmitting	73, 74
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

    Gmax  = 20 * np.log10(Dk) + 7.7
    G1    = 2 + 15 * np.log10(Dk)
    phi_m = 20 * kD * np.sqrt(Gmax-G1)
    phi_b = 48
    phi_r = 15.85*kD**0.6 if bCase else 100*kD


    if phi_b < phi_r : 
        print("Warning:: Phib () is less than Phir ().")

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

def ap30B( freq , phi , D ,
	  bVerbose = False , antenna_efficiency = 0.7 ) :
    #
    desc_ = """Appendix 30B reference Earth station antenna pattern.
Recommendation ITU-R S.580-6 reference Earth station antenna pattern.

    Name \t: APEREC015V01
    Type \t: Earth station, Receiving and Transmitting
    S580 \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APEREC015V01.pdf

Appendix 30B Earth station antenna pattern since WRC-03 applicable for D/lambda > 50.
Pattern is extended for D/lambda < 50 as in Appendix 8.
Pattern is extended for angles greater than 20 degrees as in Recommendation ITU-R S.465-5.
Pattern is extended in the main-lobe range as in Appendix 7 to produce continuous curves.
BR software sets antenna efficiency to 0.7 for technical examination.

Receiving    \t 605
Transmitting \t	606
    """
    if bVerbose :
        print ( desc_ )
    #
    # CONSTANTS
    # lambda is a reserved word in python so we use k for wavelength (not wavenumber)
    k     = c0/(freq*10**9)
    kD    = k / D
    Dk    = 1 / kD
    bCase = Dk >= 50
    phi   = 0 if phi<0 else 180 if phi>180 else phi

    Gmax  = 10 * np.log10( Dk**2 * np.pi**2 * antenna_efficiency )
    G1    = 2 + 15 * np.log10( Dk ) if Dk<50 else -21 + 25 * np.log10(Dk) if Dk<100 else -1 + 15 * np.log10(Dk)

    phi_m = 20 * kD * np.sqrt( Gmax - G1 )
    phi_b = 10 ** ( 42/25 )
    phi_r = 15.85*kD**0.6 if Dk>=100 else 100*kD

    if Gmax < G1 :
        print ("ERROR: UNSUPPORTED Gmax<G1 COMPLEX")
        coefA = 29
    if phi_b < phi_r : 
        print ("ERROR: Phib () is less than Phir ().")

    constants = {True:tuple((29,0,-10)),False:tuple((52,-10,10))}[bCase]

    phicase =	1*int(phi < phi_m) + 2*int(phi >= phi_m and phi < phi_r) +\
		3*int(phi_r <= phi and phi < phi_b) * int(not bCase) +\
                ( 3*int(phi_r <= phi and phi < 19.95) + 5*int(19.95 <= phi and phi < phi_b)  ) * int( bCase ) +\
		4*int(phi_b<phi)

    match  phicase :
        case 1:
            G = Gmax - 2.5 * 10**(-3) * (Dk * phi)**2
        case 2:
            G = G1
        case 3:
            G = constants[0] + constants[1] * np.log10(Dk) - 25 * np.log10(phi)
        case 4:
            G = constants[2] + constants[1] * np.log10(Dk)
        case 5:
            G = np.min([ -3.5 , 32 - 25*np.log10(phi) ])
        case _ :
            print("ERROR: Bad case value")
    return (G)

def ap30Baeq29( freq , phi , D , coefA = 29 ,
		bVerbose = False , antenna_efficiency = 0.7 ) :
    #
    desc_ = """Appendix 30B reference Earth station pattern with the improved
side-lobe for coefficient A = 29.

    Name \t: APERR_002V01
    Type \t: Earth station, Receiving and Transmitting
    AP30B A = 29  \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APERR_002V01.pdf

Appendix 30B Earth station antenna reference pattern applicable for D/lambda > 100. It is used for the determination of
coordination requirements and interference assessment in FSS Plan.

Pattern contains an optional improved near side-lobe (coefA=29) which may be used if so desired by administrations,
particularly in the cases where an aggregate C/I ratio of 26 dB cannot be obtained.

Pattern is extended for D/lambda < 100 as in Appendix 8.

Original Plan was based on the antennas having diameter 7 m for the 6/4 GHz band and 3 m for the 13/10-11 GHz band
and the antenna efficiency of 0.7.

WRC-03 replaced this Appendix 30B reference antenna pattern for coefA=32 by pattern APEREC015V01 (RR-2003). This
pattern (APERR_002V01) is still used as improved side-lobe Appendix 30B reference antenna pattern with coefA=29 for
D/lambda > 100.

BR software sets antenna efficiency to 0.7 for technical examination.
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

    if Dk > 100 :
        coefA = 29

    Gmax  = 10 * np.log10( Dk**2 * np.pi**2 * antenna_efficiency )
    G1    = 15 * np.log10( Dk ) - 30 + coefA
    phi_m = 20 * kD * np.sqrt( Gmax - G1 )
    phi_b = 10 ** ( (coefA + 10)/25 )
    phi_r = 15.85*kD**0.6 if bCase else 100*kD

    if coefA != 29 and coefA != 32 :
        print("WARNING: UNSUPPORTED coefA SETTING TO 29")
        coefA = 29
    if phi_b < phi_r : 
        print("Warning:: Phib () is less than Phir ().")

    constants = {True:tuple((0,0,-10)),False:tuple((20,-10,10))}[bCase]

    phicase = 1*int(phi < phi_m) + 2*int(phi >= phi_m and phi < phi_r) + 3*int(phi_r <= phi and phi < phi_b) + 4*int(phi_b<phi)

    match  phicase :
        case 1:
            G = Gmax - 2.5 * 10**(-3) * (Dk * phi)**2
        case 2:
            G = G1
        case 3:
            G = coefA + constants[0] + constants[1] * np.log10(Dk) - 25 * np.log10(phi)
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

def s672( Gmax , phi , phi0 , Ls = –20 ):
    desc_ = """Recommendation ITU-R S.672-4 space station antenna pattern for
use as a designe objective in the FSS employing geostationary
satellites. For Ls = -20 dB. For Ls = -10 dB.

    Name \t: APSREC408V01 APSREC407V01
    Type \t: Single feed circular beams. 
	    Ls = -10 dB for use in article 22, table 22-2.
		Ls = -20 dB for use in article 22, tables 22-2 and 22-3.

     S.672 \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APSREC408V01.pdf
           \t: https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/APSREC407V01.pdf

    https://www.itu.int/dms_pubrec/itu-r/rec/s/R-REC-S.672-4-199709-I!!PDF-E.pdf
    """
    a = {-20:2.58,-10:1.83}[Ls]
    b = 6.32
	pp0 = phi/phi0
    if pp0 < 0 :
        print ("ERROR: Not recognized value")
    if  pp0 <= a/2 :
        G = Gmax - 12 * pp0**2
    elif pp0 <= b/2 :
        G = Gmax + Ls
    elif b/2 < pp0 :
        G = Gmax + Ls + 20 - 25 * np.log10( 2*pp0 )
    if G < 0 :
        G = 0
    return G


def s1528( Gmax , phi ):
    desc_ = """ Recommendation ITU-R S.1528-0 space station antenna pattern for
non-GSO satellites operating in FSS below 30 GHz. Recommends 1.2

    Name \t: APSREC409V01
    Type \t: Space station, Receiving and Transmitting
    S1528\t:  https://www.itu.int/en/ITU-R/software/Documents/ant-pattern/APL_DOC_BY_PATTERN_NAME/REC-1528.pdf

    https://www.itu.int/dms_pubrec/itu-r/rec/s/R-REC-S.1528-0-200106-I!!PDF-E.pdf
    """
    Ln	= -15 dB
    z	= 1

    a	= 2.58
    b	= 6.32
    Dk		= 10**( (Gmax-7.7)/20 )
    psib	= np.sqrt(1200) * Dk
    Y		= b*psib*10**( 0.04*(Gmax-15) )

    G1,G2,G3,G4,G5 = 0,0,0,0,0
    if phi>0 and phi <= a*psib :
        G1 = Gmax - 3 * (phi/psib)**1.5
    elif phi <= b * psib :
        G2 = Gmax - 15
    elif  phi <= Y :
        G3 = Gmax - 15 - 25 * np.log10 (phi/b/psib)
    elif phi <= 90 :
        G4 = 0
    elif phi<=180 :
        G5 = 0.25 * Gmax
    G = np.max([G1, G2, G3, G4, G5])
    return ( G )


if __name__=='__main__':

    print("ANTENNAS")
    help()
    print ( p465( freq=2 , phi=15 , D=1.0 , bVerbose=True ) )
