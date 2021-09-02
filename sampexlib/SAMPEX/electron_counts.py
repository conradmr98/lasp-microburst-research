# reading in electron counts file: day 137 of year 1994

t = []
r1 = []
r2 = []
r3 = []
r4 = []
counter = 0

with open(r"C:\Users\jacob\Documents\CU_Boulder\Fall2020\Research\LASP Research\Lauren Blum\SAMPEX\hhrr1994137.txt", "r") as f:
    for line in f:
        if counter > 0:
            f_list = line.split()   # split() function takes a string and splits elements by existing white space
            t.append(float(f_list[0]))
            r1.append(int(f_list[1]))
            r2.append(int(f_list[2]))
            r3.append(int(f_list[3]))
            r4.append(int(f_list[4]))
        counter += 1
        


import matplotlib.pyplot as plt
import numpy as np

# plotting electron counts on day 137 of year 1994 on a linear plot

plt.plot(t, r1, 'b-', t, r4, 'r-')
plt.xlabel('seconds')
plt.ylabel('count')
plt.grid(True)
plt.show()


# plotting electron counts on day 137 of year 1994 on a log plot plot

plt.semilogy(t, r1, 'b-', t, r4, 'r-')
plt.xlabel('seconds')
plt.ylabel('count')
plt.grid(True)
plt.show()

