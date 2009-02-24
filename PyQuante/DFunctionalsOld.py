"""\
 DFunctionals.py Exchange and Correlation Density Functionals

 This is the old version, that contains code prior to Ann Mattsson's
 modifications, and is being kept around for archival purposes only.
 Putting this comment in 2007-02; should delete the file after 2007-08.
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution.

"""

from math import sqrt,pi,exp,log,atan
from NumWrap import zeros,arcsinh,where

def XC(dens,gamma,**opts):
    """\
    Simplified top level routine for all XC functionals.
    """
    functional = opts.get('functional','SVWN')
    assert functional in xfuncs.keys() and functional in cfuncs.keys()
    npts = len(dens)
    exc = zeros(npts,'d')
    vxc = zeros(npts,'d')
    if xfuncs[functional]:
        ex,vx = xfuncs[functional](dens,gamma)
        vxc = vxc + vx
        exc = exc + ex
    if cfuncs[functional]:
        ec,vc = cfuncs[functional](dens,gamma)
        vxc = vxc + vc
        exc = exc + ec
    return exc,vxc

def XCNEW(dens,gamma,**opts):
    """\
    New top level routine for all XC functionals.
    dens is a two x npts component matrix with spin-up and spin-down 
    densities.
    gamma is a three x npts component matrix with grad(n1)*grad(n2) 
    with n1 and n2 being (spin-up, spin-up), (spin-up, spin-down), 
    and (spin-down,spin-down) densities. Should contain the right 
    number of elements (for example zeros) even if gradients are not needed.
    For non-spin calculations set 
    spin-up density = spin-down density = density/2 and
    gammaupup=gammaupdown=gammadowndown = |grad(rho)|**2/4.
    For non-spin we should get dfxcdna=dfxcdnb and dfxcdgaa=dfxcdgbb.
    dfxdgab = 0 always, it's included for consistency with dfcdgab.
    AEM June 2006.
    """
    functional = opts.get('functional','SVWN')
    derivative = opts.get('derivative','analyt')
    assert functional in xfuncsn.keys() and functional in cfuncsn.keys()
    #that npts is the same for all 5 vectors should be checked elsewhere    
    npts = len(dens[0]) 
    fxc = zeros(npts,'d')
    dfxcdna = zeros(npts,'d')
    dfxcdnb = zeros(npts,'d')
    dfxcdgaa = zeros(npts,'d')
    dfxcdgab = zeros(npts,'d')
    dfxcdgbb = zeros(npts,'d')
    if xfuncsn[functional]:
        if derivative == 'analyt' and analyt[functional]:
            fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb = \
	                      xfuncsn[functional](dens,gamma)
        else:
	    fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb = \
	                      xfuncsn[functional](dens,gamma)	
	    #print fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb
	    dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb = \
	                      numder('x',functional,dens,gamma)
	    #print fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb
        fxc = fxc + fx
	dfxcdna = dfxcdna + dfxdna
        dfxcdnb = dfxcdnb + dfxdnb
        dfxcdgaa = dfxcdgaa + dfxdgaa
        dfxcdgab = dfxcdgab + dfxdgab
        dfxcdgbb = dfxcdgbb + dfxdgbb
    if cfuncsn[functional]:
        if derivative == 'analyt' and analyt[functional]:
            fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb = \
	                      cfuncsn[functional](dens,gamma)
        else:
            fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb = \
	                      cfuncsn[functional](dens,gamma)
	    #print fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb
	    dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb = \
	                      numder('c',functional,dens,gamma)
	    #print fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb
        fxc = fxc + fc
	dfxcdna = dfxcdna + dfcdna
        dfxcdnb = dfxcdnb + dfcdnb
        dfxcdgaa = dfxcdgaa + dfcdgaa
        dfxcdgab = dfxcdgab + dfcdgab
        dfxcdgbb = dfxcdgbb + dfcdgbb
    return fxc,dfxcdna,dfxcdnb,dfxcdgaa,dfxcdgab,dfxcdgbb

def S(dens,gamma=None):
    """
    Slater exchange functional for a vector of densities. Based upon
    the classic Slater exchange (which actually came from Dirac). See
    JC Slater 'The self consistent field for molecules and solids',
    McGraw Hill, New York, 1974.

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.
    """
    npts = len(dens)
    ex = zeros(npts,'d')
    vx = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i]) # Density of the alpha spin
        exa,vxa = xs(rho)
        ex[i] = 2*exa
        vx[i] = vxa
    return ex,vx

def SN(dens,gamma=None):
    """
    Slater exchange functional for vectors of spin-up and spin-down densities. 
    Based upon the classic Slater exchange (which actually came from Dirac). 
    See JC Slater 'The self consistent field for molecules and solids',
    McGraw Hill, New York, 1974.

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.

    Also note that this exchange is written as a sum but doesn't use the
    spin-scaling relationship used for 'physics' exchange (LDA, PBE, PW91).

    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    fx = zeros(npts,'d')
    dfxdna = zeros(npts,'d')
    dfxdnb = zeros(npts,'d')
    dfxdgaa = zeros(npts,'d')
    dfxdgab = zeros(npts,'d')
    dfxdgbb = zeros(npts,'d')
    for i in range(npts):
        na = float(dens[0][i]) # Density of the alpha spin
	nb = float(dens[1][i]) # Density of the beta spin
        fxa,vxa = xs(na)
	fxb,vxb = xs(nb)
        fx[i] = fxa + fxb
        dfxdna[i] = vxa
	dfxdnb[i] = vxb
    return fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb

def VWN(dens,gamma=None):
    """
    Vosko-Wilk-Nusair correlation functional for a vector of
    densities. From 'Accurate spin-dependent electron liquid
    correlation energies for local spin density calculations:
    a critical analysis.' SH Vosko, L Wilk, M Nusair, Can J
    Phys, 58, 1200 (1980).

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.
    """
    npts = len(dens)
    ec = zeros(npts,'d')
    vc = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i]) # Density of the alpha spin
        ecab,vca,vcb = cvwn(rho,rho)
        ec[i] = ecab
        vc[i] = vca
    return ec,vc

def VWNN(dens,gamma=None):
    """
    Vosko-Wilk-Nusair correlation functional for vectors of
    spin up and spin down densities. From 'Accurate spin-dependent 
    electron liquid correlation energies for local spin density 
    calculations: a critical analysis.' SH Vosko, L Wilk, M Nusair, 
    Can J Phys, 58, 1200 (1980).

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.

    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    fc = zeros(npts,'d')
    dfcdna = zeros(npts,'d')
    dfcdnb = zeros(npts,'d')
    dfcdgaa = zeros(npts,'d')
    dfcdgab = zeros(npts,'d')
    dfcdgbb = zeros(npts,'d')
    for i in range(npts):
        na = float(dens[0][i]) # Density of the alpha spin
	nb = float(dens[1][i]) # Density of the beta spin
        fcab,vca,vcb = cvwn(na,nb)
        fc[i] = fcab
	dfcdna[i] = vca
        dfcdnb[i] = vcb
    return fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb

def PW(dens,gamma=None):
    """
    Perdew-Wang correlation functional for a vector of
    densities. From 'Accurate and simple analytical representation
    of the electron-gas correlation energy' 
    John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.
    
    AEM June 2006.
    """
    npts = len(dens)
    ec = zeros(npts,'d')
    vc = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i]) # Density of the alpha spin
        ecab,vca,vcb = pw(rho,rho)
        ec[i] = ecab
        vc[i] = vca
    return ec,vc

def PWN(dens,gamma=None):
    """
    Perdew-Wang correlation functional for vectors of
    spin up and spin down densities. From 'Accurate and simple 
    analytical representation of the electron-gas correlation energy' 
    John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).

    Note that gamma is ignored, and is only included for consistency
    with the other XC functionals.

    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    fc = zeros(npts,'d')
    dfcdna = zeros(npts,'d')
    dfcdnb = zeros(npts,'d')
    dfcdgaa = zeros(npts,'d')
    dfcdgab = zeros(npts,'d')
    dfcdgbb = zeros(npts,'d')
    for i in range(npts):
        na = float(dens[0][i]) # Density of the alpha spin
	nb = float(dens[1][i]) # Density of the beta spin
        fcab,vca,vcb = pw(na,nb)
        fc[i] = fcab
	dfcdna[i] = vca
        dfcdnb[i] = vcb
    return fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb

def B(dens,gamma):
    """
    Becke 1988 Exchange Functional. From 'Density-functional exchange
    energy approximation with correct asymptotic behavior.' AD Becke,
    PRA 38, 3098 (1988).
    """
    npts = len(dens)
    assert len(gamma) == npts
    ex = zeros(npts,'d')
    vx = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i])  # Density of the alpha spin
        gam = 0.25*gamma[i]
        exa,vxa = xb(rho,gam)
        ex[i] = 2*exa
        vx[i] = vxa
    return ex,vx

def BN(dens,gamma):
    """
    Becke 1988 Exchange Functional. From 'Density-functional exchange
    energy approximation with correct asymptotic behavior.' AD Becke,
    PRA 38, 3098 (1988).
    Chemistry way.
    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    assert len(gamma[0]) == npts
    assert len(gamma[1]) == npts
    assert len(gamma[2]) == npts
    fx = zeros(npts,'d')
    dfxdna = zeros(npts,'d')
    dfxdnb = zeros(npts,'d')
    dfxdgaa = zeros(npts,'d')
    dfxdgab = zeros(npts,'d')
    dfxdgbb = zeros(npts,'d')
    for i in range(npts):
        na = float(dens[0][i])  # Density of the alpha spin
	nb = float(dens[1][i])  # Density of the beta spin
        gamaa = gamma[0][i]
	gambb = gamma[2][i]
        fxa,fxna,fxgaa = xb(na,gamaa,return_flag = 1)
	fxb,fxnb,fxgbb = xb(nb,gambb,return_flag = 1)
        fx[i] = fxa + fxb
        dfxdna[i] = fxna
	dfxdnb[i] = fxnb
        dfxdgaa[i] = fxgaa
        dfxdgbb[i] = fxgbb
    return fx,dfxdna,dfxdnb,dfxdgaa,dfxdgab,dfxdgbb


def LYP(dens,gamma):
    """Transformed version of LYP. See 'Results obtained with correlation
    energy density functionals of Becke and Lee, Yang, and Parr.' Miehlich,
    Savin, Stoll and Preuss. CPL 157, 200 (1989).
    """
    npts = len(dens)
    ec = zeros(npts,'d')
    vc = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i]) # Density of the alpha spin
        gam = 0.25*gamma[i]
        ecab,vca,vcb = clyp(rho,rho,gam,gam,gam)
        ec[i] = ecab
        vc[i] = vca
    return ec,vc

def LYPN(dens,gamma):
    """Transformed version of LYP. See 'Results obtained with correlation
    energy density functionals of Becke and Lee, Yang, and Parr.' Miehlich,
    Savin, Stoll and Preuss. CPL 157, 200 (1989).
    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    assert len(gamma[0]) == npts
    assert len(gamma[1]) == npts
    assert len(gamma[2]) == npts
    fc = zeros(npts,'d')
    dfcdna = zeros(npts,'d')
    dfcdnb = zeros(npts,'d')
    dfcdgaa = zeros(npts,'d')
    dfcdgab = zeros(npts,'d')
    dfcdgbb = zeros(npts,'d')
    for i in range(npts):
        na = float(dens[0][i]) # Density of the alpha spin
	nb = float(dens[1][i]) # Density of the beta spin
        gamaa = gamma[0][i]
	gamab = gamma[1][i]
        gambb = gamma[2][i]
        fcab,fcna,fcnb,fcgaa,fcgab,fcgbb = clyp(na,nb,gamaa,gamab,gambb,return_flag=1)
        fc[i] = fcab
	dfcdna[i] = fcna
        dfcdnb[i] = fcnb
        dfcdgaa[i] = fcgaa
	dfcdgab[i] = fcgab
        dfcdgbb[i] = fcgbb
    return fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb

def XPBE(dens,gamma):
    "PBE Exchange Functional"
    npts = len(dens)
    assert len(gamma) == npts
    ex = zeros(npts,'d')
    vx = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i])  # Density of the alpha spin
        gam = 0.25*gamma[i]
        exa,vxa = xpbe(rho,gam)
        ex[i] = 2*exa
        vx[i] = vxa
    return ex,vx

def CPBE(dens,gamma):
    "PBE Correlation Functional"
    npts = len(dens)
    ec = zeros(npts,'d')
    vc = zeros(npts,'d')
    for i in range(npts):
        rho = 0.5*float(dens[i]) # Density of the alpha spin
        gam = 0.25*gamma[i]
        ecab,vca,vcb = cpbe(rho,rho,gam,gam,gam)
        ec[i] = ecab
        vc[i] = vca
    return ec,vc

def EXXC1(dens,gamma):
    "AEM's EXX compatible correlation #1 (note: no spin). AEM June 2006."
    npts = len(dens)
    ec = zeros(npts,'d')
    vc = zeros(npts,'d')
    for i in range(npts):
        rho = float(dens[i]) # Density 
        gam = gamma[i]
        ecpnt,dnedrho,dnedgamma = c1(rho,gam)
        ec[i] = ecpnt
        vc[i] = dnedrho   # more derivatives need to be added
    return ec,vc

def EXXC1N(dens,gamma):
    "AEM's EXX compatible correlation #1 (note: no spin). AEM June 2006."
    npts = len(dens[0])
    assert len(dens[1]) == npts
    assert len(gamma[0]) == npts
    assert len(gamma[1]) == npts
    assert len(gamma[2]) == npts
    fc = zeros(npts,'d')
    dfcdna = zeros(npts,'d')
    dfcdnb = zeros(npts,'d')
    dfcdgaa = zeros(npts,'d')
    dfcdgab = zeros(npts,'d')
    dfcdgbb = zeros(npts,'d')
    for i in range(npts):
        rho = float(dens[0][i]+dens[1][i]) # Total density 
        gam = gamma[0][i]+gamma[2][i] + 2.0*gamma[1][i] # Total gamma
        fcpnt,dfdrho,dfdgamma = c1(rho,gam)
        fc[i] = fcpnt
        dfcdna[i] = dfdrho 
	dfcdnb[i] = dfdrho
	dfcdgaa[i] = dfdgamma
	dfcdgab[i] = 2.0*dfdgamma
        dfcdgbb[i] = dfdgamma
    return fc,dfcdna,dfcdnb,dfcdgaa,dfcdgab,dfcdgbb

def AM05(dens,gamma):
    """Armiento and Mattsson functional from 2005. (note: no spin)
    Rickard Armiento and Ann E Mattsson, PRB 72, 085108 (2005).
    AEM June 2006.
    """
    npts = len(dens[0])
    assert len(dens[1]) == npts
    assert len(gamma[0]) == npts
    assert len(gamma[1]) == npts
    assert len(gamma[2]) == npts
    fxc = zeros(npts,'d')
    dfxcdna = zeros(npts,'d')
    dfxcdnb = zeros(npts,'d')
    dfxcdgaa = zeros(npts,'d')
    dfxcdgab = zeros(npts,'d')
    dfxcdgbb = zeros(npts,'d')
    for i in range(npts):
        rho = float(dens[0][i]+dens[1][i]) # Total density 
        gam = gamma[0][i]+gamma[2][i] + 2.0*gamma[1][i] # Total gamma
        fpnt,dfdrho,dfdgamma = am05xc(rho,gam)
        fxc[i] = fpnt
        dfxcdna[i] = dfdrho 
	dfxcdnb[i] = dfdrho
	dfxcdgaa[i] = dfdgamma
	dfxcdgab[i] = 2.0*dfdgamma
        dfxcdgbb[i] = dfdgamma
    return fxc,dfxcdna,dfxcdnb,dfxcdgaa,dfxcdgab,dfxcdgbb    

# Functional terms themselves.
# Functionals are defined in their spin-polarized versions. However,
# since the exchange functionals factor (e.g. fx(rhoa,rhob) = fx(rhoa)+fx(rhob))
# we define them in the separate form, understanding that we can always
# use the above relationship to refactor them
def xs(rho,**opts):
    "Xalpha X functional. Can pass in 'alpha' as a kwarg"
    tol = opts.get('tol',1e-10)
    Xalpha = opts.get('Xalpha',2./3.)
    fac=-2.25*Xalpha*pow(3./4./pi,1./3.)
    if rho < tol: rho=0
    rho3 = pow(rho,1./3.)
    ex = fac*rho*rho3
    vx = (4./3.)*fac*rho3
    return ex,vx

def xb(rho,gam,**opts):
    # Corrected by AEM June 2006.
    tol = opts.get('tol',1e-10)
    return_flag = opts.get('return_flag',0)
    fx=dfxdrho=dfxdgam=0
    if rho > tol:
        rho13 = pow(rho,1./3.)
        x = sqrt(gam)/rho13/rho
        g = b88_g(x)
        dg = b88_dg(x)
        dfxdrho = (4./3.)*rho13*(g-x*dg)
        if gam > tol: dfxdgam = 0.5*dg/sqrt(gam)
        fx = rho*rho13*g
    if return_flag == 1: return fx,dfxdrho,dfxdgam
    return fx,dfxdrho

def xpbe(rho,gam,**opts):
    tol = opts.get('tol',1e-10)
    return_flag = opts.get('return_flag',0)
    kap = 0.804
    mu = 0.449276922095889E-2
    ex=vxrho=vxgam=0
    if rho > tol:
        ex0,vx0 = xs(rho)
        rho13 = rho**(1.E0/3.E0)
        rho43 = rho13*rho
        den = 1.E0+mu*gam/rho43/rho43
        F = 1+kap-kap/den
        ex = ex0*F
        dFdr = -(8./3.)*kap*mu*gam/den/den*rho**(-11./3.)
        vxrho = vx0*F+ex0*dFdr
        dFdg = -kap*mu/rho43/rho43/den/den
        vxgam = ex0*dFdg
    if return_flag == 1: return ex,vxrho,vxgam
    return ex,vxrho


def cvwn(rhoa,rhob,**opts):
    tol = opts.get('tol',1e-10)
    rho = rhoa+rhob
    ec = vcrhoa = vcrhob = 0
    if rho < tol: return ec,vcrhoa,vcrhob
    zeta=(rhoa-rhob)/rho
    x = pow(3./4./pi/rho,1./6.)
    epsp = vwn_epsp(x)
    depsp = vwn_depsp(x)
    g = vwn_g(zeta)
    # Can uncomment for a little more speed:
    if g < tol: return epsp*rho,epsp-(x/6.)*depsp,epsp-(x/6.)*depsp
    epsf = vwn_epsf(x)
    eps = epsp + g*(epsf-epsp)
    ec = eps*rho
    depsf = vwn_depsf(x)
    dg = vwn_dg(zeta)
    deps_dx = depsp + g*(depsf-depsp)
    deps_dg = (epsf-epsp)*dg
    vcrhoa = eps - (x/6.)*deps_dx + deps_dg*(1-zeta)
    vcrhob = eps - (x/6.)*deps_dx - deps_dg*(1+zeta)
    return ec,vcrhoa,vcrhob

def clyp(rhoa,rhob,gamaa,gamab,gambb,**opts):
    # Modified and corrected by AEM in June 2006.
    tol = opts.get('tol',1e-10)
    return_flag = opts.get('return_flag',0)
    a = 0.04918  # Parameters from the LYP papers
    b = 0.132
    c = 0.2533
    d = 0.349
    rho = rhoa+rhob
    fc=fcrhoa=fcrhob=fcgamaa=fcgamab=fcgambb=0
    assert rhoa >= 0.0
    assert rhob >= 0.0
    if rho > tol:
        rhom3 = pow(rho,-1./3.)
        w = exp(-c*rhom3)/(1+d*rhom3)*pow(rho,-11./3.)
        dl = c*rhom3+d*rhom3/(1+d*rhom3)

        fcgamaa = -a*b*w*((1./9.)*rhoa*rhob*(1-3*dl-(dl-11)*rhoa/rho)-rhob*rhob)
        fcgamab = -a*b*w*((1./9.)*rhoa*rhob*(47-7*dl)-(4./3.)*rho*rho)
	fcgambb = -a*b*w*((1./9.)*rhoa*rhob*(1-3*dl-(dl-11)*rhob/rho)-rhoa*rhoa)

        fc = -4*a/(1+d*rhom3)*rhoa*rhob/rho \
             -pow(2,11./3.)*0.3*pow(3*pi*pi,2./3.)*a*b*w \
             *rhoa*rhob*(pow(rhoa,8./3.)+pow(rhob,8./3.)) \
             + fcgamaa*gamaa + fcgamab*gamab + fcgambb*gambb

        dw = -(1./3.)*pow(rho,-4./3.)*w*(11*pow(rho,1./3.)-c-d/(1+d*rhom3))
        ddl = (1./3.)*(d*d*pow(rho,-5./3.)/pow(1+d*rhom3,2)-dl/rho)

        d2f_dradgaa = dw/w*fcgamaa - a*b*w*(
            (1./9.)*rhob*(1-3*dl-(dl-11)*rhoa/rho)
            -(1./9.)*rhoa*rhob*((3+rhoa/rho)*ddl+(dl-11)*rhob/rho/rho))
        d2f_dradgbb = dw/w*fcgambb - a*b*w*(
            (1./9.)*rhob*(1-3*dl-(dl-11)*rhob/rho)
            -(1./9.)*rhoa*rhob*((3+rhob/rho)*ddl-(dl-11)*rhob/rho/rho)
            -2*rhoa)
        d2f_dradgab = dw/w*fcgamab-a*b*w*(
            (1./9)*rhob*(47-7*dl)-(7./9.)*rhoa*rhob*ddl-(8./3.)*rho)

	d2f_drbdgaa = dw/w*fcgamaa - a*b*w*(
            (1./9.)*rhoa*(1-3*dl-(dl-11)*rhoa/rho)
            -(1./9.)*rhoa*rhob*((3+rhoa/rho)*ddl-(dl-11)*rhoa/rho/rho)
            -2*rhob)
	d2f_drbdgbb = dw/w*fcgambb - a*b*w*(
            (1./9.)*rhoa*(1-3*dl-(dl-11)*rhob/rho)
            -(1./9.)*rhoa*rhob*((3+rhob/rho)*ddl+(dl-11)*rhoa/rho/rho))
        d2f_drbdgab = dw/w*fcgamab-a*b*w*(
            (1./9)*rhoa*(47-7*dl)-(7./9.)*rhoa*rhob*ddl-(8./3.)*rho)

        fcrhoa = -4*a/(1+d*rhom3)*rhoa*rhob/rho*(
            (1./3.)*d*pow(rho,-4./3.)/(1+d*rhom3)+1/rhoa-1/rho)\
            -pow(2,11./3.)*0.3*pow(3*pi*pi,2./3.)*a*b*(
            dw*rhoa*rhob*(pow(rhoa,8./3.)+pow(rhob,8./3.))
            +w*rhob*((11./3.)*pow(rhoa,8./3.)+pow(rhob,8./3.))) \
            +d2f_dradgaa*gamaa + d2f_dradgbb*gambb + d2f_dradgab*gamab
        fcrhob = -4*a/(1+d*rhom3)*rhoa*rhob/rho*(
            (1./3.)*d*pow(rho,-4./3.)/(1+d*rhom3)+1/rhob-1/rho)\
            -pow(2,11./3.)*0.3*pow(3*pi*pi,2./3.)*a*b*(
            dw*rhoa*rhob*(pow(rhob,8./3.)+pow(rhoa,8./3.))
            +w*rhoa*((11./3.)*pow(rhob,8./3.)+pow(rhoa,8./3.))) \
            +d2f_drbdgaa*gamaa + d2f_drbdgbb*gambb + d2f_drbdgab*gamab
    if return_flag == 1: return fc,fcrhoa,fcrhob,fcgamaa,fcgamab,fcgambb
    return fc,fcrhoa,fcrhob

def cpbe(rhoa,rhob,gama,gamb,gamab,**opts):
    # My attempt to rewrite the functional. Not finished yet.
    tol = opts.get('tol',1e-10)
    return_flag = opts.get('return_flag',0)
    rho = rhoa+rhob
    ec = vca = vcb = vcgama = vcgamb = vcgamab = 0
    gam = 0.031091
    ohm = 0.046644
    bet = 0.066725
    if rho > tol:
        Rs = pow(3./(4.*pi*rho),1./3.)
        Zeta = (rhoa-rhob)/rho
        Kf = pow(3*pi*pi*rho,1./3.)
        Ks = sqrt(4*Kf/pi)
        Phi = 0.5*(pow(1+Zeta,2./3.) + pow(1-Zeta,2./3.))
        Phi3 = Phi*Phi*Phi
        gradrho = sqrt(gama+gamb+2.*gamab)
        T = gradrho/(2*Phi*Ks*rho)
        T2 = T*T
        T4 = T2*T2

        eps,vc0a,vc0b = cpbe_lsd(rhoa,rhob)

        expo = (exp(-eps/(gam*Phi3))-1.)
        A = bet/gam/expo
        N = T2+A*T4
        D = 1.+A*T2+A*A*T4
        H = gam*Phi3*log(1.+(bet/gam)*N/D)
        ec = rho*(eps+H)

        # Derivative stuff
        dZ_drhoa = (1.-Zeta)/rho
        dZ_drhob = -(1.+Zeta)/rho

        dPhi_dZ = pow(1.+Zeta,-1./3.)/3.-pow(1.-Zeta,-1./3.)/3.
        dPhi_drhoa = dPhi_dZ*dZ_drhoa
        dPhi_drhob = dPhi_dZ*dZ_drhob
        
        dKs_drho = Ks/(6*rho)
        
        dT_dPhi = -T/Phi
        dT_dKs = -T/Ks
        dT_drhoa = -T/rho + dT_dPhi*dPhi_drhoa + dT_dKs*dKs_drho
        dT_drhob = -T/rho + dT_dPhi*dPhi_drhob + dT_dKs*dKs_drho

        dA_dPhi = -A/expo*exp(-eps/(gam*Phi3))*(3*eps/(gam*Phi3*Phi))
        dA_deps = -A/expo*exp(-eps/(gam*Phi3))*(-1/(gam*Phi3))
        deps_drhoa = (vc0a-eps)/rho
        deps_drhob = (vc0b-eps)/rho
        dA_drhoa = dA_dPhi*dPhi_drhoa + dA_deps*deps_drhoa
        dA_drhob = dA_dPhi*dPhi_drhob + dA_deps*deps_drhoa

        dN_dT = 2*T+4*A*T2*T
        dD_dT = 2*A*T + 4*A*A*T*T2
        dN_dA = T4
        dD_dA = T2+2*A*T4

        dH_dPhi = 3*H/Phi
        dH_dT = bet*Phi3/(1.+bet/gam*N/D)*(D*dN_dT-N*dD_dT)/D/D
            
        dH_dA = bet*Phi3/(1.+bet/gam*N/D)*(D*dN_dA-N*dD_dA)/D/D
        
        dH_drhoa = dH_dPhi*dPhi_drhoa + dH_dT*dT_drhoa + dH_dA*dA_drhoa
        dH_drhob = dH_dPhi*dPhi_drhob + dH_dT*dT_drhob + dH_dA*dA_drhob
        
        vca = vc0a + H + rho*dH_drhoa
        vcb = vc0b + H + rho*dH_drhob
    # Havent done the dE_dgamma derives yet
    return ec,vca,vcb

def c1(rho,gam,**opts):
    # AEM's EXX compatible correlation 1.
    # AEM June 2006.
    tol = opts.get('tol',1e-10)
    return_flag = opts.get('return_flag',0)
    fc = dfcdrho = dfcdgamma = 0
    g1 = 0.3060
    g2 = 0.04108
    t = 0.3123
    if rho > tol:
	rs = pow(3./(4.*pi*rho),1./3.)
	snorm2 = (4.*pow(3.*pi**2,2./3.)*pow(rho,8./3.))
	s2 = gam / snorm2
	X = 2.*exp(-s2/t)/(1.+exp(-s2/t))
    	eps,vc0a,vc0b = cpbe_lsd(0.5*rho,0.5*rho)
    	fc = rho*eps*(X + (1. - X)*(g1 + g2*rs))
	#derivatives
	dXds2 = -X/t/(1.+exp(-s2/t))
	dfcdrho = vc0a*(X + (1. - X)*(g1 + g2*rs)) - \
		eps/3.*(g2*(1. - X)*rs + 8*(1.-g1-g2*rs)*s2*dXds2)
	dfcdgamma = eps*(1.-g1-g2*rs)*rho/snorm2*dXds2
    return fc,dfcdrho,dfcdgamma

def pw(rhoa,rhob,**opts):
    # AEM June 2006.
    tol = opts.get('tol',1e-10)
    rho = rhoa+rhob
    ec = vca = vcb = 0
    if rho < tol: return ec,vca,vcb
    eps,vca,vcb = cpbe_lsd(rhoa,rhob)
    ec = rho*eps
    return ec,vca,vcb
    
def cpbe_lsd(rhoa,rhob):
    # Not quite VWN. AEM: It's usually called PW correlation
    # LSD terms
    # Note that this routine gives out ec, not fc.
    # If you rather have fc, use pw instead
    rho = rhoa+rhob
    Rs = pow(3./(4.*pi*rho),1./3.)
    Zeta = (rhoa-rhob)/rho
    thrd = 1./3.     # thrd*=various multiples of 1/3
    thrd4 = 4*thrd
    ggam=0.5198420997897463295344212145565 # gam= 2^(4/3)-2
    fzz=8./(9.*ggam) # fzz=f''(0)= 8/(9*gam)
    rtrs = sqrt(Rs)
    eu,eurs = pbe_gcor(0.0310907,0.21370,7.5957,
                       3.5876,1.6382,0.49294,rtrs)
    ep,eprs = pbe_gcor(0.01554535,0.20548,14.1189,
                       6.1977,3.3662,0.62517,rtrs)
    alfm,alfrsm = pbe_gcor(0.0168869,0.11125,10.357,
                           3.6231,0.88026,0.49671,rtrs)
    alfc = -alfm
    z4 = Zeta**4
    f=((1.+Zeta)**thrd4+(1.-Zeta)**thrd4-2.)/ggam
    eps = eu*(1.-f*z4)+ep*f*z4-alfm*f*(1.-z4)/fzz

    ecrs = eurs*(1.-f*z4)+eprs*f*z4-alfrsm*f*(1.-z4)/fzz
    fz = thrd4*((1.+Zeta)**thrd-(1.-Zeta)**thrd)/ggam
    eczet = 4.*(Zeta**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu-(1.-z4)*alfm/fzz)
    comm = eps -Rs*ecrs/3.-Zeta*eczet
    vca = comm + eczet
    vcb = comm - eczet
    return eps,vca,vcb
    
def pbe_gcor(a,a1,b1,b2,b3,b4,rtrs):
#      subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
# slimmed down version of gcor used in pw91 routines, to interpolate
# lsd correlation energy, as given by (10) of
# j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
# k. burke, may 11, 1996.
#      implicit real*8 (a-h,o-z)
      q0 = -2.*a*(1.+a1*rtrs*rtrs)
      q1 = 2.*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = log(1.+1./q1)
      gg = q0*q2
      q3 = a*(b1/rtrs+2.*b2+rtrs*(3.*b3+4.*b4*rtrs))
      ggrs = -2.*a*a1*q2-q0*q3/(q1*(1.+q1))
      return gg,ggrs

def vwn_xx(x,b,c): return x*x+b*x+c
def vwn_epsp(x): return vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352)
def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_eps(x,a,x0,b,c):
    q = sqrt(4*c-b*b)
    eps = a*(log(x*x/vwn_xx(x,b,c))
             - b*(x0/vwn_xx(x0,b,c))*log(pow(x-x0,2)/vwn_xx(x,b,c))
             + (2*b/q)*(1-(x0*(2*x0+b)/vwn_xx(x0,b,c))) * atan(q/(2*x+b)))
    return eps

def vwn_depsp(x): return vwn_deps(x,0.0310907,-0.10498,3.72744,12.9352)
def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_deps(x,a,x0,b,c):
    q = sqrt(4*c-b*b)
    deps = a*(2/x - (2*x+b)/vwn_xx(x,b,c)
              - 4*b/(pow(2*x+b,2)+q*q) - (b*x0/vwn_xx(x0,b,c))
              * (2/(x-x0)-(2*x+b)/vwn_xx(x,b,c)-4*(2*x0+b)/(pow(2*x+b,2)+q*q)))
    return deps

def vwn_g(z): return 1.125*(pow(1+z,4./3.)+pow(1-z,4./3.)-2)
def vwn_dg(z): return 1.5*(pow(1+z,1./3.)-pow(1-z,1./3.))

def b88_g(x,b=0.0042):
    return -1.5*pow(3./4./pi,1./3.)-b*x*x/(1.+6.*b*x*arcsinh(x))

def b88_dg(x,b=0.0042):
    num = 6*b*b*x*x*(x/sqrt(x*x+1)-arcsinh(x))-2*b*x
    denom = pow(1+6*b*x*arcsinh(x),2)
    return num/denom

def pbe_F(s,mu=0.2195149727645171,kap=0.8040):
    P0 = 1.+mu*s*s/kap
    return 1.+kap-kap/P0

def pbe_dF(s,mu=0.2195149727645171,kap=0.8040):
    P0 = 1.+mu*s*s/kap
    return 2*mu*s/(P0*P0)

def am05xc(rho,gam,**opts):
    # AM05, both exchange and correlation in one routine
    # AEM June 2006.
    tol = opts.get('tol',1e-16) # to conform with official am05 routine
    fxc = dfxcdrho = dfxcdgamma = 0
    g = 0.8098
    a = 2.804
    c = 0.7168
    if rho > tol:
	snorm2 = (4.*pow(3.*pi**2,2./3.)*pow(rho,8./3.))
	s2 = abs(gam) / snorm2
	s = pow(s2,1./2.)
	# LDAPW exchange and correlation
	fx0,vxlda = xs(0.5*rho) # xs(na) = fx(na) = na*ex(2*na)
	fxlda = 2.0*fx0
    	fclda,vclda,vc0b = pw(0.5*rho,0.5*rho)
	assert vclda == vc0b
	# Interpolation index
	X = 1.0/(1.0 + a*s2)
	w = am05_lambertw(pow(s,3./2.)/sqrt(24.0))
	if (s < 1.e-14):     # am05_lambertw give back argument if it is < 1.0e-20
              zosn = 1.0     # (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
      	else:
              zosn = 24.**(1./3.)*pow(w,2./3.)/s   # zosn = normalized z/s
      	zfac = s2*(zosn*27./32./pi**2)**2
	denom = 1.0 + c*s2*zosn*pow(1.0 + zfac,1./4.) # denom = denominator of Airy LAA refinement function
      	F = (c*s2 + 1.0)/denom   # Airy LAA refinement function
	# Refinement functions
	Hx = X + (1.0 - X)*F 
	Hc = X + g*(1.0 - X)
	# Exchange-correlation energy density, Exc = Integrate[fxc]
	fxc = fxlda*Hx + fclda*Hc
	# Derivatives
	# Interpolation index derivatives: 1/s dX/ds
      	Xsos = -2.0*a*X**2
      	szsoz = 1.0/(1.0 + w)  # szsoz = s*(dz/ds)/z 
	Fsos = c/denom**2*(2.0 - zosn*     # Airy LAA refinement function derivatives, 1/s dF/ds
                 ((1.0 - c*s2)*pow(1.0 + zfac,1./4.) +
                  (1.0 + c*s2)*(1.0 + 3./2.*zfac)/
                           pow(1.0 + zfac,3./4.)*szsoz)) 
	# Refinement function derivatives, 1/s dH{x,c}/ds
	Hxsos = (1.0 - X)*Fsos - (F - 1.0)*Xsos
      	Hcsos = Xsos*(1.0 - g)
	# Exchange-correlation energy density derivatives
	dfxcdrho = vxlda*Hx + vclda*Hc - 4./3.*s2/rho*(Hcsos*fclda +
			Hxsos*fxlda)
	dfxcdgamma = 1./2./snorm2*(Hcsos*fclda + Hxsos*fxlda)
    return fxc,dfxcdrho,dfxcdgamma

def am05_lambertw(z):
    # Used only in am05xc. AEM June, 2006.
    assert z >= 0.0
    # If z small, go with the first term of the power expansion, z
    if (z < 1.e-20):
       result = z
    else:
       e = exp(1.0)
       # Inital guess
       if (abs(z + 1.0/e) > 1.45 ):
          # Asymptotic expansion at 0 and Inf
          result = log(z)
          result = result - log(result)
       else:
          # Series expansion about -1/e to first order
          result = sqrt(2.0*e*z + 2.0) - 1.0
       # Find result through iteration
       for i in range(1,11):
          p = exp(result)
          t = result*p - z
          if (result != -1.0):
             t = t/(p*(result + 1.0) -
                  0.5*(result + 2.0)*t/(result + 1.0))
          else:
             t = 0.0
          result = result - t
          if (abs(t) < (2.48*1.e-14)*(1.0 + abs(result))): break
          assert i != 10
    return result
    
def numder(xc,functional,dens,gamma):
        """\
        Since numerical derivatives will not be used for production
	I do not care if this routine is fast or not. AEM June 2006.
        """
	maxdelta = 1.e-5
	mindelta = 1.e-14
	npts = len(dens[0]) #that npts is the same for all 5 vectors should be checked elsewhere
    	dfdna = zeros(npts,'d')
    	dfdnb = zeros(npts,'d')
    	dfdgaa = zeros(npts,'d')
    	dfdgab = zeros(npts,'d')
    	dfdgbb = zeros(npts,'d')
	delta1 = dens[0]/1000.0 #spin up
	delta = where(delta1 > maxdelta, maxdelta, delta1)
	dens1 = zeros(dens.shape,'d')
	dens1[0] = dens[0] + delta
	dens1[1] = dens[1]
	f1 = getf(xc,functional,dens1,gamma)
	dens1[0] = dens[0] - delta
	f2 = getf(xc,functional,dens1,gamma)
	delta1 = delta
	delta =  where(delta1 < mindelta, mindelta, delta1) # to avoid division by zero
	dfdna = (f1 - f2)/2.0/delta
	delta1 = dens[1]/1000.0 # spin down
	delta = where(delta1 > maxdelta, maxdelta, delta1)
	dens1[0] = dens[0] 
	dens1[1] = dens[1] + delta
	f1 = getf(xc,functional,dens1,gamma)
	dens1[1] = dens[1] - delta
	f2 = getf(xc,functional,dens1,gamma)
	delta1 = delta
	delta =  where(delta1 < mindelta, mindelta, delta1) # to avoid division by zero
	dfdnb = (f1 - f2)/2.0/delta
        delta1 = gamma[0]/1000.0 #grad spin up*grad spin up
	delta = where(delta1 > maxdelta, maxdelta, delta1)
	gamma1 = zeros(gamma.shape,'d')
	gamma1[0] = gamma[0] + delta
	gamma1[1] = gamma[1]
	gamma1[2] = gamma[2]
	f1 = getf(xc,functional,dens,gamma1)
	gamma1[0] = gamma[0] - delta
	f2 = getf(xc,functional,dens,gamma1)
	delta1 = delta
	delta =  where(delta1 < mindelta, mindelta, delta1) # to avoid division by zero
	dfdgaa = (f1 - f2)/2.0/delta
        delta1 = gamma[1]/1000.0 #grad spin up*grad spin down
	delta = where(delta1 > maxdelta, maxdelta, delta1)
	gamma1[0] = gamma[0] 
	gamma1[1] = gamma[1] + delta
	gamma1[2] = gamma[2]
	f1 = getf(xc,functional,dens,gamma1)
	gamma1[1] = gamma[1] - delta
	f2 = getf(xc,functional,dens,gamma1)
	delta1 = delta
	delta =  where(delta1 < mindelta, mindelta, delta1) # to avoid division by zero
	dfdgab = (f1 - f2)/2.0/delta
        delta1 = gamma[2]/1000.0 #grad spin down*grad spin down
	delta = where(delta1 > maxdelta, maxdelta, delta1)
	gamma1[0] = gamma[0] 
	gamma1[1] = gamma[1]
	gamma1[2] = gamma[2] + delta
	f1 = getf(xc,functional,dens,gamma1)
	gamma1[2] = gamma[2] - delta
	f2 = getf(xc,functional,dens,gamma1)
	delta1 = delta
	delta =  where(delta1 < mindelta, mindelta, delta1) # to avoid division by zero
	dfdgbb = (f1 - f2)/2.0/delta
	return dfdna,dfdnb,dfdgaa,dfdgab,dfdgbb

def getf(xc,functional,dens,gamma):
	# AEM June, 2006.
	if (xc == 'x'):
		f,t1,t2,t3,t4,t5 = xfuncsn[functional](dens,gamma)
	elif (xc == 'c'):
		f,t1,t2,t3,t4,t5 = cfuncsn[functional](dens,gamma)	
	return f
	
# Table mapping functional names into functions etc.
# LDA -> SVWN, GGA -> PBE
# xfuncsn, cfuncsn, and analyt added by AEM in June 2006.
# PW, AM05 and EXXC1 added by AEM in June 2006.
xfuncs = dict(LDA=S,S0=S,SVWN=S,BLYP=B,PBE=XPBE,LYP=None,CPBE=None,EXXC1=None,VWN=None,PW=None,LDAPW=S)
cfuncs = dict(LDA=None,S0=None,SVWN=VWN,BLYP=LYP,PBE=CPBE,LYP=LYP,CPBE=CPBE,EXXC1=EXXC1,VWN=VWN,PW=PW,LDAPW=PW)
xfuncsn = dict(LDA=SN,S0=SN,SVWN=SN,BLYP=BN,LYP=None,VWN=None,PW=None,LDAPW=SN,AM05=AM05)
cfuncsn = dict(LDA=None,S0=None,SVWN=VWNN,BLYP=LYPN,LYP=LYPN,VWN=VWNN,PW=PWN,LDAPW=PWN,AM05=None)
analyt = dict(LDA=True,S0=True,SVWN=True,BLYP=True,PBE=False,
              LYP=True,CPBE=False,EXXC1=False,VWN=True,PW=True,LDAPW=True,AM05=True)
need_gradients = dict(LDA=False,S0=False,SVWN=False,BLYP=True,PBE=True,
                      LYP=True,CPBE=True,EXXC1=True,VWN=False,PW=False,LDAPW=False,AM05=True)


