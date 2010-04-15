"""\

NAME
      qmmd.py

SYNOPSIS
      First principles molecular dynamics module
      
DESCRIPTION
      This is an experimental module for ab initio MD calculations which 
      makes use of numerical integrators to move nuclei in response to forces
      computed through analytic derivatives of the total energy expression.
	  
AUTHOR
      Hatem H. Helal, hhh23@cam.ac.uk

REPORT BUGS
      Report bugs to hhh23@cam.ac.uk

COPYRIGHT

"""
from leapfrog import leapfrog
from hartree_fock import rhf,uhf,uhf_fixed_occ
from Wavefunction import Wavefunction
from force import hf_force

from Element import symbol

def rhf_dyn(atoms,**opts):
    """\
    Uses RHF derived forces to compute dynamics.  
    
    Options:      Value   Description
    --------      -----   -----------
    job           pydyn   Descriptive job name
    nsteps        100     number of dynamics steps to take
    dt            0.1     time step size in picoseconds
                          units matter! we assume atom positions are 
                          stored in bohrs, velocities are in bohr/ps,
                          acceleration is in bohr/ps^2
                          and forces are in hartree/bohr
    
    
    Hartree-Fock options copied from hartree_fock.py -> rhf
    rhf(atoms,**opts) - Closed-shell HF driving routine
    
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    ETemp         False   Use ETemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    """
    #dynamics options
    job = opts.get('job','pydyn')
    nsteps = opts.get('nsteps',100)
    dt = opts.get('dt',0.1)
    
    #save any given RHF options
    cc = opts.get('ConvCriteria',1e-4)
    maxit = opts.get('MaxIter',20)
    doavg = opts.get('DoAveraging',False)
    temp = opts.get('ETemp',False)
    bfcns = opts.get('bfs',None)
    if not bfcns:
        bdat = opts.get('basis_data',None)
    ints = opts.get('integrals', None)
    init_orbs = opts.get('orbs',None)
    
    #open data file to store energy info 
    edat = open(job+'.edat', 'w')
    edat.write("#Step       Time    PE  KE  TE\n")
    
    #open trajectory file to store xyz info
    xyz = open(job+'.xyz', 'w')
    #xyz.write("#RHF molecular dynamics done by PyQuante\n")
    #xyz.write("#job: %s  nsteps: %d  dt:%f\n"%(job,nsteps,dt))
    xyz.write(xyz_str(atoms))
    t=0.
    for n in xrange(nsteps):
        t+=n*dt
        pe,orben,coefs = rhf(atoms,ConvCriteria=cc,MaxIter=maxit,\
                           DoAveraging=doavg,ETemp=temp,bfs=bfcns,\
                           basis_data=bdat,integrals=ints,orbs=init_orbs)

        ncl,nop = atoms.get_closedopen()

        wf = Wavefunction(orbs=coefs,orbe=orben,restricted=True,nclosed=ncl,nopen=nop)
        hf_force(atoms,wf,bdat)
        
        ke = leapfrog(atoms,t,dt)
        te = ke+pe
        bl = atoms[0].dist(atoms[1])
        edat.write('%d      %f	%f      %f      %f  %f\n' %(n,t,bl,pe,ke,te))
        xyz.write(xyz_str(atoms))

    edat.close()
    xyz.close()

    return 

def uhf_dyn(atoms,**opts):
    """\
    Uses UHF derived forces to compute dynamics.  
    
    Options:      Value   Description
    --------      -----   -----------
    job           pydyn   Descriptive job name
    nsteps        100     number of dynamics steps to take
    dt            0.1     time step size in picoseconds
                          units matter! we assume atom positions are 
                          stored in bohrs, velocities are in bohr/ps,
                          acceleration is in bohr/ps^2
                          and forces are in hartree/bohr
    
    
    Hartree-Fock options copied from hartree_fock.py -> uhf
    uhf(atoms,**opts) - Unrestriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS averaging for convergence acceleration
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not None, the guess orbitals
    """
    #dynamics options
    job = opts.get('job','pydyn')
    nsteps = opts.get('nsteps',100)
    dt = opts.get('dt',0.1)
    
    #save any given UHF options
    cc = opts.get('ConvCriteria',1e-4)
    maxit = opts.get('MaxIter',20)
    doavg = opts.get('DoAveraging',False)
    temp = opts.get('ETemp',False)
    bfcns = opts.get('bfs',None)
    if not bfcns:
        bdat = opts.get('basis_data',None)
    ints = opts.get('integrals', None)
    init_orbs = opts.get('orbs',None)
    
    #open data file to store energy info 
    edat = open(job+'.edat', 'w')
    edat.write("#Step       Time    PE  KE  TE\n")
    
    #open trajectory file to store xyz info
    xyz = open(job+'.xyz', 'w')
    #xyz.write("#RHF molecular dynamics done by PyQuante\n")
    #xyz.write("#job: %s  nsteps: %d  dt:%f\n"%(job,nsteps,dt))
    xyz.write(xyz_str(atoms))
    t=0.
    for n in xrange(nsteps):
        t+=n*dt
        pe,(orbea,orbeb),(coefsa,coefsb) = uhf(atoms,ConvCriteria=cc,MaxIter=maxit,\
                                               DoAveraging=doavg,ETemp=temp,bfs=bfcns,\
                                               basis_data=bdat,integrals=ints,orbs=init_orbs)

        na,nb = atoms.get_alphabeta()

        wf = Wavefunction(orbs_a=coefsa,orbs_b=coefsb,\
                          orbe_a=orbea,orbe_b=orbeb,\
                          unrestricted=True,nalpha=na,nbeta=nb)
        hf_force(atoms,wf,bdat)
        
        #bl = atoms[0].dist(atoms[1])  #testing with H2

        ke = leapfrog(atoms,t,dt)
        te = ke+pe
        
        edat.write('%d      %f	%f      %f      %f\n' %(n,t,pe,ke,te))
        xyz.write(xyz_str(atoms))

    edat.close()
    xyz.close()

    return 


def fixedocc_uhf_dyn(atoms,occa,occb,**opts):
    """\
    Uses Fixed occupation UHF derived forces to compute dynamics.  

    occa and occb represent the orbital occupation arrays for 
    the calculating spin orbitals with holes

    Options:      Value   Description
    --------      -----   -----------
    job           pydyn   Descriptive job name
    nsteps        100     number of dynamics steps to take
    dt            0.1     time step size in picoseconds
                          units matter! we assume atom positions are 
                          stored in bohrs, velocities are in bohr/ps,
                          acceleration is in bohr/ps^2
                          and forces are in hartree/bohr
    
    
    Hartree-Fock options copied from hartree_fock.py -> uhf
    uhf(atoms,**opts) - Unrestriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS averaging for convergence acceleration
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not None, the guess orbitals
    """
    #dynamics options
    job = opts.get('job','pydyn')
    nsteps = opts.get('nsteps',100)
    dt = opts.get('dt',0.1)
    
    #save any given UHF options
    cc = opts.get('ConvCriteria',1e-4)
    maxit = opts.get('MaxIter',20)
    doavg = opts.get('DoAveraging',False)
    temp = opts.get('ETemp',False)
    bfcns = opts.get('bfs',None)
    if not bfcns:
        bdat = opts.get('basis_data',None)
    ints = opts.get('integrals', None)
    init_orbs = opts.get('orbs',None)
    
    #open data file to store energy info 
    edat = open(job+'.edat', 'w')
    edat.write("#Step       Time    PE  KE  TE\n")
    
    #open trajectory file to store xyz info
    xyz = open(job+'.xyz', 'w')
    #xyz.write("#RHF molecular dynamics done by PyQuante\n")
    #xyz.write("#job: %s  nsteps: %d  dt:%f\n"%(job,nsteps,dt))
    xyz.write(xyz_str(atoms))
    t=0.
    for n in xrange(nsteps):
        t+=n*dt
        pe,(orbea,orbeb),(coefsa,coefsb) = uhf_fixed_occ(atoms,occa,occb, ConvCriteria=cc,MaxIter=maxit,\
                                               DoAveraging=doavg,ETemp=temp,bfs=bfcns,\
                                               basis_data=bdat,integrals=ints,orbs=init_orbs)

        na,nb = atoms.get_alphabeta()

        wf = Wavefunction(orbs_a=coefsa,orbs_b=coefsb,\
                          orbe_a=orbea,orbe_b=orbeb,\
                          occs_a=occa,occs_b=occb,\
                          fixedocc=True,nalpha=na,nbeta=nb)
        hf_force(atoms,wf,bdat)
        
        bl = atoms[0].dist(atoms[1])  #testing with H2
        
        ke = leapfrog(atoms,t,dt)
        te = ke+pe
    
        edat.write('%d      %f  %f	%f      %f      %f\n' %(n,t,bl,pe,ke,te))
        xyz.write(xyz_str(atoms))

    edat.close()
    xyz.close()

    return 


def xyz_str(mol):
    """
    Takes a mol and returns the xyz file as a string
    """
    natoms = len(mol.atoms)
    str='%d\n' %natoms
    str+='%s\n' %mol.name
    for atom in mol.atoms:
        sym = symbol[atom.atno]
        xyz = atom.r
        tmpstr = "%s    %12.6f  %12.6f  %12.6f\n" %(sym,xyz[0],xyz[1],xyz[2])
        str += tmpstr
    return str
