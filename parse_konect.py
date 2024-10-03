import random
import math

def parse_konect(name):

    f_in = open("data/" + name + "_org.txt" , 'r')
    f_out = open("data/" + name + ".txt", 'w')

    for line in f_in.readlines():
        lst = line.split()

        u = lst[0]
        v = lst[1]
        t = lst[3]

        f_out.write(u + " " + v +  " " + t + "\n")

parse_konect("wikiDyn_it")