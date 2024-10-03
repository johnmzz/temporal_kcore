import random
import math

def parse(name):

    f_in = open("data/" + name + ".txt" , 'r')
    f_out = open("data/" + name + "_sorted.txt" , 'w')
    edges = []

    for line in f_in.readlines():
        lst = line.split()

        u = lst[0]
        v = lst[1]
        t = int(lst[2])

        edge = [t, u, v]
        edges.append(edge)

    edges.sort()
    
    for e in edges:
        f_out.write(e[1] + " " + e[2] + " " + str(e[0]) + "\n")

parse("bitcoin_alpha")