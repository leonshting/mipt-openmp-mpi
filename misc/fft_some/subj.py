import numpy as np
import matplotlib.pyplot as plt
import random
import math
import cmath

def f(x):
    return int(np.abs(10-x))


pi2 = cmath.pi * 2.0
def DFT(fnList):
    N = len(fnList)
    FmList = []
    for m in range(N):
        Fm = 0.0
        for n in range(N):
            Fm += fnList[n] * cmath.exp(- 1j * pi2 * m * n / N)
        FmList.append(Fm / N)
    return FmList

def InverseDFT(FmList):
    N = len(FmList)
    fnList = []
    for n in range(N):
        fn = 0.0
        for m in range(N):
            fn += FmList[m] * cmath.exp(1j * pi2 * m * n / N)
        fnList.append(fn)
    return fnList


N = 100
a = 2.0
f = 2.0
p = float(random.randint(0, 360))
print "frequency = " + str(f)
print "amplitude = " + str(a)
print "phase ang = " + str(p)
T = 1
X = []
freqs = []
fnList = []
for n in range(N):
    t = float(n) / N * pi2
    fn = a * math.sin(f * t + p / 360 * pi2) + math.sin(2*f * t + p / 360 * pi2)
    fnList.append(fn)
    X.append(f * t + p / 360 * pi2)
    freqs.append(n/T)



plt.plot(X, [i.real for i in fnList])
plt.show()
FmList = DFT(fnList)
plt.plot(freqs,[i.real for i in FmList])
plt.plot(freqs,[i.imag for i in FmList])
plt.show()
plt.plot(X, [i.real for i in InverseDFT(FmList)])
plt.show()

