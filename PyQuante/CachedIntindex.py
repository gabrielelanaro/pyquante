# Experiment to see if caching the two electron integral indices
# would speed things up. It didn't speed things up appreciably,
# if at all, but I'm keeping the code around for reference.

class CachedIntindex:
    # Yuck! N4 storage
    def __init__(self):
        self.cache = {}

    def initialize(self,nbf):
        # This is even slower...
        #from numpy import zeros
        #self.cache = zeros((nbf,nbf,nbf,nbf),'i')
        for i in xrange(nbf):
            for j in xrange(i+1):
                ij = i*(i+1)/2+j
                for k in xrange(nbf):
                    for l in xrange(k+1):
                        kl = k*(k+1)/2+l
                        if ij < kl: continue
                        val = ij*(ij+1)/2+kl
                        self.cache[i,j,k,l] = val
                        self.cache[i,j,l,k] = val
                        self.cache[j,i,k,l] = val
                        self.cache[j,i,l,k] = val
                        self.cache[k,l,i,j] = val
                        self.cache[l,k,i,j] = val
                        self.cache[k,l,j,i] = val
                        self.cache[l,k,j,i] = val
        return

    def __call__(self,i,j,k,l): return self.cache[i,j,k,l]
                        
class CachedIndexList:
    def initialize(self,nbf):
        self.nbf = nbf
        self.intindex = CachedIntindex()
        self.intindex.initialize(nbf)
        self.jints = {}
        self.kints1 = {}
        self.kints2 = {}
        for i in xrange(self.nbf):
            for j in xrange(i+1):
                self.jints[i,j] = []
                self.kints1[i,j] = []
                self.kints2[i,j] = []
                for k in xrange(nbf):
                    for l in xrange(nbf):
                        self.jints[i,j].append(self.intindex(i,j,k,l))
                        self.kints1[i,j].append(self.intindex(i,k,j,l))
                        self.kints2[i,j].append(self.intindex(i,l,j,k))
        return
