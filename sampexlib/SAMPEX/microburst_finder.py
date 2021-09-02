
# reading in electron counts file: day 137 of year 1994

t = []
r1 = []
r2 = []
r3 = []
r4 = []
counter = 0

with open(r"C:\Users\jacob\Documents\CU_Boulder\Fall2020\Research\LASP Research\Lauren Blum\SAMPEX\hhrr1994137.txt", "r") as f:
    # with open(): creates a block that opens and closes a file automatically
    for line in f:
        if counter > 0:
            f_list = line.split()   # split() function takes a string and splits elements by existing white space
            t.append(float(f_list[0]))
            r1.append(int(f_list[1]))
            r2.append(int(f_list[2]))
            r3.append(int(f_list[3]))
            r4.append(int(f_list[4]))
        counter += 1
        
# detect when microbursts occur using the microburst algorithm

import math
import numpy as np

N_100 = np.add(np.add(r1,r2),np.add(r3,r4))  # total counts from all 4 SSDs every 100ms

A_500 = []  # vector of average counts every 500ms
i = 0

for i in range(0, len(N_100)-4, 5):
    avg = ( N_100[i] + N_100[i+1] + N_100[i+2] + N_100[i+3] + N_100[i+4] ) / 5
    A_500.append(avg)

t_microburst = []
i = 0
ii = 0

while i < len(N_100):
    algorithm = ( N_100[i] - A_500[ii] ) / math.sqrt(1 + A_500[ii])
    if  algorithm > 10:
        t_microburst.append(t[i])
        
    i += 1
    if (i % 5) == 0:   # the % symbol finds the remainder of the division (i / 5)
        ii += 1

print(t_microburst)

        