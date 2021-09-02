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
        

# calulation of isotropy index of electron counts

import numpy as np

r1 = np.array(r1)
r4 = np.array(r4)
iso_indices = []

for i, ii in zip(r1, r4):
    if ii != 0:
        iso_index = np.abs(ii - i) / ii                    
        iso_indices = np.append(iso_indices, iso_index)   
    if ii == 0 and i != 0:
        iso_index = np.abs(i - ii) / i
        iso_indices = np.append(iso_indices, iso_index)

iso_index = np.average(iso_indices)
print(f" Isotropy index = {iso_index:.4f}")