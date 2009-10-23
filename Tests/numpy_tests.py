#!/usr/bin/env python
def numpy_test():
    from numpy import array,identity
    from numpy import dot as matrixmultiply
    from numpy.linalg import eigh as Heigenvalues

    A = array([[1,2],[2,3]],'d')
    I = identity(2,'d')
    print "A = \n",A
    print "I = \n",I
    print "A*I = \n",matrixmultiply(A,I)
    print "eig(A) = ",Heigenvalues(A)
    return

def numeric_test():
    from Numeric import array,matrixmultiply,identity
    from LinearAlgebra import Heigenvalues

    A = array([[1,2],[2,3]],'d')
    I = identity(2,'d')
    print "A = \n",A
    print "I = \n",I
    print "A*I = \n",matrixmultiply(A,I)
    print "eig(A) = ",Heigenvalues(A)
    return

def test():
    try:
        numpy_test()
    except:
        print "Numpy test failed; probably import error"

    try:
        numeric_test()
    except:
        print "Numeric test failed; probably import error"

if __name__ == '__main__': test()

