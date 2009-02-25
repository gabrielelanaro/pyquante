#!/usr/bin/env python
from PyQuante.DFunctionals import cpbe,cpbe_huub,cpbe_rpm
from math import sqrt,exp,log,pi

# The cpbe_huub and cpbe_rpm are now in NewStuff/OtherPBE.py

def main():
    data = [(1.,1.,0.,0.,0.),
            (1.,0.1,0.,0.,0.),
            (1.,1.,0.5,0.5,0.5),
            (1.,0.1,0.5,0.5,0.5),            
            ]
    for rhoa,rhob,ga,gb,gab in data:
        ec,vca,vcb = cpbe_huub(rhoa,rhob,ga,gb,gab)
        ec2 = cpbe_rpm(rhoa,rhob,ga,gb,gab)
        ec3,vca3,vcb3 = cpbe(rhoa,rhob,ga,gb,gab)
        print ec,ec2,ec3
    return

if __name__ == '__main__': main()

