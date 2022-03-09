__author__ = 'bigbrother'
import random
def a(b):
    strret = ""
    for i in b:
        strret += str(b) + "\t"
    return strret+"\n"
file = open("input.txt", "w")
x = [int(input()) for i in range(4)]

file.write(a(x))
for j in range(x[0]):
    for k in range(x[1]):
        for l in range(x[2]):
            file.write(str(random.randint(0,1))+" "+str(j)+" " + str(k)+ " "+ str(l)+ "\n");