# A formula that works of Boys Function
from scipy.integrate import quad
from scipy import exp
def Fm(m,x):
    return quad(lambda t: t**(2*m) * exp(-x*t**2), 0, 1)[0]
