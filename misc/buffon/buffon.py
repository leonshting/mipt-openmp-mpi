__author__ = 'leonsht'

import random
import sys
def kx(x1,y1,x2,y2):
    return (y2-y1)/(x2-x1)

def x_closest(x_c, freq):
    a = x_c % freq
    if a > freq/2:
        return freq - a
    else:
        return a


needle_length = 0.001
num_lines = 100
counter = 0
if((1./num_lines) >= needle_length):
    for i in range(10000000):
        rand_cx = random.random()
        rand_cy = random.random()
        rand_somex = random.random()
        rand_somey = random.random()
        k = kx(rand_cx, rand_cy, rand_somex, rand_somey)
        x_dist = x_closest(rand_cx, 1./float(num_lines))
        if ((k*x_dist)**2 + x_dist**2) <= needle_length**2:
            counter +=1
        pi = 2 * float(2*needle_length*(i+1))/((1./num_lines) * (counter+1))
        sys.stdout.write("PI estimation is = %lf   \r" % (pi) )
        sys.stdout.flush()
