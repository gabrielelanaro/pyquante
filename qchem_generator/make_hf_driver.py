#!/usr/bin/env python
"Print out some basis information for the water molecule"

from PyQuante.Molecule import Molecule
from PyQuante.Ints import getbasis

# Routines for testing only: not required for normal operation
from Numeric import array2string as a2s
from PyQuante.Ints import getT,getS,getV,get2ints,getints,getJ,getK
from PyQuante.cints import ijkl2intindex
from PyQuante.LA2 import GHeigenvectors,mkdens,TraceProperty

def get_adat(atoms):
    "Break the molecule down into data arrays for printing"
    atnos = []
    xs = []
    ys = []
    zs = []
    for atom in atoms:
        atno,(x,y,z) = atom.atuple()
        atnos.append(atno)
        xs.append(x)
        ys.append(y)
        zs.append(z)
    return atnos,xs,ys,zs

def get_cdat(bfs):
    "Break the contracted basis data into arrays for printing"
    xcenter = []
    ycenter = []
    zcenter = []
    lpower = []
    mpower = []
    npower = []
    normc = []
    istart = []
    nprim = []
    i = 0
    nbf = len(bfs)
    for ibf in range(nbf):
        bf = bfs[ibf]
        x,y,z = bf.origin()
        xcenter.append(x)
        ycenter.append(y)
        zcenter.append(z)
        l,m,n = bf.powers()
        lpower.append(l)
        mpower.append(m)
        npower.append(n)
        normc.append(bf.norm())
        istart.append(i)
        nprim.append(len(bf._prims))
        i += len(bf._prims)
    return xcenter,ycenter,zcenter,lpower,mpower,npower,normc,istart,nprim

def get_pdat(bfs):
    "Break the primitive data into arrays for printing"
    alpha = []
    coef = []
    normp = []
    iprim = 0
    nbf = len(bfs)
    for ibf in range(nbf):
        bf = bfs[ibf]
        for prim in bf._prims:
            alpha.append(prim.exp())
            coef.append(prim.coef())
            normp.append(prim.norm())
            iprim += 1
    return alpha,coef,normp

def output_elements(array,tag,type):
    "Print all elements of array to file with appropriate decoration"
    lines = []
    for i in range(len(array)):
        element = array[i]
        if type == 'i':
            lines.append("  %s[%d] = %d;" % (tag,i,element))
        else:
            lines.append("  %s[%d] = %f;" % (tag,i,element))
    return "\n".join(lines)

def output_elements_table(array,tag,type):
    lines = ["  %s %s[] = {" % (type,tag)]
    for el in array:
        lines.append("    %s," % el)
    lines.append("  };")
    return "\n".join(lines)

def make_hf_driver(atoms,**opts):
    basis_data = opts.get('basis_data',None)
    do_print_one = opts.get('do_print_one',False)
    do_print_two = opts.get('do_print_two',False)
    maxiter = opts.get('maxiter',5)
    
    nat = len(atoms)
    bfs = getbasis(atoms,basis_data)
    nbf = len(bfs)
    nbf2 = nbf*nbf
    nint = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
    nclosed,nopen = atoms.get_closedopen()
    enuke = atoms.get_enuke()
    assert nopen==0, "code only works for closed-shell systems"

    nprim = 0
    for bf in bfs: nprim += len(bf._prims)

    atnos,xs,ys,zs = get_adat(atoms)
    atno_array = output_elements_table(atnos,'atno','int')
    x_array = output_elements_table(xs,'x','double')
    y_array = output_elements_table(ys,'y','double')
    z_array = output_elements_table(zs,'z','double')

    # Print out data over contracted bfs:
    xcenter,ycenter,zcenter,lpower,mpower,npower,normcs,istarts,nprims = get_cdat(bfs)
    xcenter_array = output_elements_table(xcenter,'xcenter','double')
    ycenter_array = output_elements_table(ycenter,'ycenter','double')
    zcenter_array = output_elements_table(zcenter,'zcenter','double')
    lpower_array = output_elements_table(lpower,'lpower','int')
    mpower_array = output_elements_table(mpower,'mpower','int')
    npower_array = output_elements_table(npower,'npower','int')
    normc_array = output_elements_table(normcs,'normc','double')
    istart_array = output_elements_table(istarts,'istart','int')
    nprim_array = output_elements_table(nprims,'nprim','int')

    # print out data over primitive bfs:
    alpha,coef,normp = get_pdat(bfs)
    alpha_array = output_elements_table(alpha,'alpha','double')
    coef_array = output_elements_table(coef,'coef','double')
    normp_array = output_elements_table(normp,'normp','double')

    template = open("main_template.c").read()

    open("%s.c" % atoms.name,'w').write(template % locals())

    return

def test_mol(mol,**opts):
    basis_data = opts.get('basis_data',None)
    do_python_tests = opts.get('do_python_tests',True)
    
    make_hf_driver(mol,basis_data=basis_data)
    if do_python_tests:
        bfs = getbasis(mol,basis_data)
        S= getS(bfs)
        T = getT(bfs)
        V = getV(bfs,mol)
        Ints = get2ints(bfs)
        nclosed,nopen = mol.get_closedopen()
        enuke = mol.get_enuke()
        assert nopen==0
  
        h = T+V

        orbe,orbs = GHeigenvectors(h,S)
        print "Eval of h: ",
        print orbe
        for i in range(10):
            D = mkdens(orbs,0,nclosed)
            J = getJ(Ints,D)
            K = getK(Ints,D)
            orbe,orbs = GHeigenvectors(h+2*J-K,S)
            eone = TraceProperty(D,h)
            ej = TraceProperty(D,J)
            ek = TraceProperty(D,K)
            energy = enuke+2*eone+2*ej-ek
            print i,energy,enuke,eone,ej,ek
    return
    

def test_h2o():
    h2o = Molecule('h2o',
                   [(8,(0.,0.,0.)),
                    (1,(1.,0.,0.)),
                    (1,(0.,1.,0.))],
                   units='Angs')
    from PyQuante.basis_sto3g import basis_data
    test_mol(h2o,basis_data=basis_data)
    return

def test_h2():
    h2 = Molecule('h2',[(1,(0.,0.,-0.35)),(1,(0.,0.,0.35))])
    test_mol(h2)
    return

def test_h2_sto():
    from PyQuante.basis_sto3g import basis_data
    h2 = Molecule('h2_sto',[(1,(0.,0.,-0.35)),(1,(0.,0.,0.35))])
    test_mol(h2,basis_data=basis_data)
    return

if __name__ == '__main__':
    test_h2o()
    test_h2()
    test_h2_sto()
    
