# Experiment to see if caching the two electron integral indices
# would speed things up. Didn't speed things appreciably.

class CachedIntindex:
    # Yuck! N4 storage
    def __init__(self,nbf):
        self.cache = {}
        # This is even slower...
        #from numpy import zeros
        #self.cache = zeros((nbf,nbf,nbf,nbf),'i')
        for i in range(nbf):
            for j in range(i+1):
                ij = i*(i+1)/2+j
                for k in range(nbf):
                    for l in range(k+1):
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
                        
