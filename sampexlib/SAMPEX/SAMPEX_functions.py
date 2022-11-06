#!/usr/bin/env python
# coding: utf-8

# In[1]:


def read_counts(file):
    import numpy as np
    
    t = []; r1 = []; r2 = []; r3 = []; r4 = []
    counter = 0

    with open(r"C:\Users\jacob\Documents\CU_Boulder\Research\LASP\lasp-research\sampexlib\SAMPEX_Data\1993_data\{}".format(file), "r")     as f:
        # with open(): creates a block that opens and closes a file automatically
        # .format(a, b) replaces {} within preceding quote with variables (starting with a, then b, etc)
        for line in f:
            if counter > 0:
                f_list = line.split()       # split() function takes a string and splits elements by existing white space
                t.append(float(f_list[0]))
                r1.append(int(f_list[1]))
                r2.append(int(f_list[2]))
                r3.append(int(f_list[3]))
                r4.append(int(f_list[4]))
            counter += 1

    return np.array(t), np.array(r1), np.array(r2), np.array(r3), np.array(r4)


def read_days(year, start_day, num_days, month):
    import numpy as np
    from SAMPEX_functions import read_counts as read
    
    # reading in electron counts files:
    file = ['']*num_days

    for i in np.arange(len(file)):
        if len(str(int(start_day) + i)) == 1:
            file[i] = 'hhrr' + year + '00' + str(int(start_day)+i) + '.txt'
        elif len(str(int(start_day) + i)) == 2:
            file[i] = 'hhrr' + year + '0' + str(int(start_day)+i) + '.txt'
        elif len(str(int(start_day) + i)) == 3:
            file[i] = 'hhrr' + year + str(int(start_day)+i) + '.txt'

    t = np.array([])        
    r1 = np.array([])
    r2 = np.array([])
    r3 = np.array([])
    r4 = np.array([])

    for i in np.arange(len(file)):
        t_i, r1_i, r2_i, r3_i, r4_i = read(file[i])
    
        t = np.append(t, (t_i+86400*i))
        r1 = np.append(r1, r1_i)
        r2 = np.append(r2, r2_i)
        r3 = np.append(r3, r3_i)
        r4 = np.append(r4, r4_i)
    
    return np.array(t), np.array(r1), np.array(r2), np.array(r3), np.array(r4)


def mb_finder(t, r1, r2, r3, r4):
    import math
    import numpy as np

    ##### Uses microburst algorithm on a data set of electorn counts to find all microbursts within time interval #####
    # calculate N_100 and A_500 values for microburst algorithm
    N_100 = np.add(np.add(r1,r2),np.add(r3,r4))  # total counts from all 4 SSDs every 100ms
    A_500 = []  # array of average counts every 500ms
    i = 0

    for i in np.arange(len(N_100)):
        if i < 2:
            avg = ( N_100[0] + N_100[1] + N_100[2] + N_100[3] + N_100[4] ) / 5
            A_500.append(avg)
        
        elif i > len(N_100)-3:
            avg = ( N_100[len(N_100)-5] + N_100[len(N_100)-4] + N_100[len(N_100)-3] + N_100[len(N_100)-2] +                    N_100[len(N_100)-1] ) / 5
            A_500.append(avg)
        
        else:
            avg = ( N_100[i-2] + N_100[i-1] + N_100[i] + N_100[i+1] + N_100[i+2] ) / 5
            A_500.append(avg)
        
    # find microburst times and N_100, SSD1, SSD4 counts using algorithm
    t_microburst = []
    N_100_microburst = []
    r1_microburst = []
    r4_microburst = []
    mb_index = []
    MB_mask = []

    for i in np.arange(len(N_100)):
        algorithm = ( N_100[i] - A_500[i] ) / math.sqrt(1 + A_500[i])
        if algorithm > 10:
            t_microburst.append(t[i])
            N_100_microburst.append(N_100[i])
            r1_microburst.append(r1[i])
            r4_microburst.append(r4[i])
            mb_index.append(i)
            MB_mask.append(True)
        else:
            MB_mask.append(False)
        
    return np.array(t_microburst), np.array(N_100_microburst), np.array(r1_microburst), np.array(r4_microburst), np.array(mb_index), np.array(MB_mask), np.array(N_100), np.array(A_500)



def mb_magnitude(N_100, N_100_microburst, A_500):
    import numpy as np
    import math as m
    import pandas as pd
    
    # calculation of B_3 bin percentiles
    B_3 = list(pd.Series(N_100).rolling(30, center=True).quantile(0.1))

    # handles boundary conditions
    percentile_start = np.percentile(N_100[0:29], 10)
    percentile_end = np.percentile(N_100[len(N_100)-29:len(N_100)], 10)

    for i in np.arange(len(N_100)):
        if np.isnan(B_3[i]) and i < 15:
            B_3[i] = percentile_start

        elif np.isnan(B_3[i]) and i > len(N_100)-15:
            B_3[i] = percentile_end

    # calculation of microburst magnitudes using B_3 bin percentiles during microbursts
    B_3_microburst = []

    for i in np.arange(len(N_100)):
        algorithm = ( N_100[i] - A_500[i] ) / m.sqrt(1 + A_500[i])
        if algorithm > 10:
            B_3_microburst.append(B_3[i])
            
    B_3_microburst = np.array(B_3_microburst)
    B_3 = np.array(B_3)
    
    y_microburst = N_100_microburst - B_3_microburst  # magnitude of microburst
    return y_microburst, B_3_microburst, B_3



def iso_calculator(r1, r4):
    import numpy as np
    import math
    
    # calculation of the isotropy indices of electron counts
    iso_indices = []        # for data analysis purposes
    #iso_indices_plot = []   # for plotting purposes

    for i in np.arange(len(r1)):
        if r4[i] >= 20 and r1[i] >= 20:
            if  r4[i] != 0 and r4[i] > r1[i]:
                iso_index = r1[i]/r4[i]
                iso_indices.append(iso_index)
                #iso_indices_plot.append(iso_index)
            elif  r1[i] != 0 and r1[i] >= r4[i]:
                iso_index = r4[i]/r1[i]
                iso_indices.append(iso_index)
                #iso_indices_plot.append(iso_index)
        else:
            # when counts are too low, they do not contribute to the isotropy day-average isotropy index
            iso_index = np.nan
            iso_indices.append(iso_index)
        
    iso_indices = np.array(iso_indices)             # only the iso indices true values
    #iso_indices_plot = np.array(iso_indices_plot)   # returns 0 for indices that do not pass the 20 count test (for plotting)
    
    return iso_indices#, iso_indices_plot



def OrbAtt_augment(t_electrons, OrbAtt_file, num_days):
    import numpy as np
    
    # Augment OrbAtt data to fit counts data
    L_Shell = OrbAtt_file['L_Shell'].values; sec_of_day = OrbAtt_file['sec_of_day'].values; MLT = OrbAtt_file['MLT'].values;
    Pitch = OrbAtt_file['Pitch'].values; Lat = OrbAtt_file['GEO_Lat'].values; Long = OrbAtt_file['GEO_Long'].values;
    Radius = OrbAtt_file['GEO_Radius'].values

    # fix sec_of_day file to correct times
    count = 0
    start_i = 0

    for i in np.arange(len(sec_of_day)):
        if sec_of_day[i+1]-sec_of_day[i] < 0:
            sec_of_day[start_i:i+1] += 86400*count
        
            count += 1
            start_i = i+1
        if count >= (num_days-1):
            sec_of_day[start_i:] += 86400*count
            break

    # trim the start of OrbAtt data to match that of the counts data
    while sec_of_day[0] < t_electrons[0]:
        sec_of_day = sec_of_day[1:]; L_Shell = L_Shell[1:]; MLT = MLT[1:]; Pitch = Pitch[1:]; Lat = Lat[1:]; Long = Long[1:];
        Radius = Radius[1:]
    # trim the start of the counts data to match that of OrbAtt data
    while t_electrons[0] < sec_of_day[0]:
        t_electrons = t_electrons[1:]
        
    # create OrbAtt arrays with the same length as microburst arrays
    t_OrbAtt = [sec_of_day[0]]; LS_OrbAtt = [L_Shell[0]]; MLT_OrbAtt = [MLT[0]];
    P_OrbAtt = [Pitch[0]]; Lat_OrbAtt = [Lat[0]]; Long_OrbAtt = [Long[0]]; R_OrbAtt = [Radius[0]]

    a = 0  # OrbAtt index
    b = 0  # Counts index
    break_key = False

    while True:
        b += 1
        if b >= len(t_electrons) or a >= len(sec_of_day):
            break
        #counter
        if a+1 < len(sec_of_day):    # this handles end conditions
        
            if t_electrons[b] - t_electrons[b-1] > 1:       # this handles breaks in the electron data (catches up the OA data index)
                while sec_of_day[a] < int(t_electrons[b]):
                    a += 1
                    if a+1 >= len(sec_of_day):
                        while len(t_OrbAtt) < len(t_electrons):
                            t_OrbAtt.append(sec_of_day[a]); LS_OrbAtt.append(L_Shell[a]); MLT_OrbAtt.append(MLT[a]); 
                            P_OrbAtt.append(Pitch[a]); Lat_OrbAtt.append(Lat[a]); Long_OrbAtt.append(Long[a]); R_OrbAtt.append(Radius[a])
                        break_key = True
                        break
                if break_key == True:
                    break
        
            currentsec = sec_of_day[a]
            nextsec = sec_of_day[a+1]
    
            if abs(currentsec - t_electrons[b]) >= abs(nextsec - t_electrons[b]):  # this increments the OA data to match the counts data
                a += 1

        t_OrbAtt.append(sec_of_day[a]); LS_OrbAtt.append(L_Shell[a]); MLT_OrbAtt.append(MLT[a]); 
        P_OrbAtt.append(Pitch[a]); Lat_OrbAtt.append(Lat[a]); Long_OrbAtt.append(Long[a]); R_OrbAtt.append(Radius[a])

    t_OrbAtt = np.array(t_OrbAtt); LS_OrbAtt = np.array(LS_OrbAtt); MLT_OrbAtt = np.array(MLT_OrbAtt);
    P_OrbAtt = np.array(P_OrbAtt); Lat_OrbAtt = np.array(Lat_OrbAtt); Long_OrbAtt = np.array(Long_OrbAtt);
    R_OrbAtt = np.array(R_OrbAtt)
    
    return t_OrbAtt, LS_OrbAtt, MLT_OrbAtt, P_OrbAtt, Lat_OrbAtt, Long_OrbAtt, R_OrbAtt



def OrbAtt_augment_loop(t_electrons, r1, r2, r3, r4, OrbAtt_file, day_of_year):
    import numpy as np
    
    # Augment OrbAtt data to fit counts data
    L_Shell = OrbAtt_file['L_Shell'].values; MLT = OrbAtt_file['MLT'].values; Pitch = OrbAtt_file['Pitch'].values;
    Lat = OrbAtt_file['GEO_Lat'].values; Long = OrbAtt_file['GEO_Long'].values; Radius = OrbAtt_file['GEO_Radius'].values

    sec_of_day = (OrbAtt_file['hr'].values)*3600 + (OrbAtt_file['min'].values)*60 + OrbAtt_file['sec'].values
    sec_of_year = ((OrbAtt_file['day'].values)-1)*86400 + sec_of_day
    
    # extract the single day of OrbAtt data
    sec_of_day = sec_of_day[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    L_Shell = L_Shell[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    MLT = MLT[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    Pitch = Pitch[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    Lat = Lat[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    Long = Long[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]
    Radius = Radius[(sec_of_year > (day_of_year-1)*86400) & (sec_of_year < day_of_year*86400)]

    # trim the start of OrbAtt data to match that of the counts data
    while sec_of_day[0] < t_electrons[0]:
        sec_of_day = sec_of_day[1:]; L_Shell = L_Shell[1:]; MLT = MLT[1:]; Pitch = Pitch[1:]; Lat = Lat[1:]; Long = Long[1:];
        Radius = Radius[1:]
    # trim the start of the counts data to match that of OrbAtt data
    while t_electrons[0] < sec_of_day[0]:
        t_electrons = t_electrons[1:]; r1 = r1[1:]; r2 = r2[1:]; r3 = r3[1:]; r4 = r4[1:]
        
    # create OrbAtt arrays with the same length as microburst arrays
    t_OrbAtt = [sec_of_day[0]]; LS_OrbAtt = [L_Shell[0]]; MLT_OrbAtt = [MLT[0]];
    P_OrbAtt = [Pitch[0]]; Lat_OrbAtt = [Lat[0]]; Long_OrbAtt = [Long[0]]; R_OrbAtt = [Radius[0]]

    a = 0  # OrbAtt index
    b = 0  # Counts index
    break_key = False

    while True:
        b += 1
        if b >= len(t_electrons) or a >= len(sec_of_day):
            #print(1)
            break
        #counter
        if a+1 < len(sec_of_day):    # this handles end conditions
        
            if t_electrons[b] - t_electrons[b-1] > 1:       # this handles breaks in the electron data (catches up the OA data index)
                while sec_of_day[a] < int(t_electrons[b]):
                    a += 1
                    if a+1 >= len(sec_of_day):
                        while len(t_OrbAtt) < len(t_electrons):
                            t_OrbAtt.append(sec_of_day[a]); LS_OrbAtt.append(L_Shell[a]); MLT_OrbAtt.append(MLT[a]); 
                            P_OrbAtt.append(Pitch[a]); Lat_OrbAtt.append(Lat[a]); Long_OrbAtt.append(Long[a]); R_OrbAtt.append(Radius[a])
                        break_key = True
                        break
                if break_key == True:
                    break
                    #print(2)
        
            currentsec = sec_of_day[a]
            nextsec = sec_of_day[a+1]
    
            if abs(currentsec - t_electrons[b]) >= abs(nextsec - t_electrons[b]):  # this increments the OA data to match the counts data
                a += 1

        t_OrbAtt.append(sec_of_day[a]); LS_OrbAtt.append(L_Shell[a]); MLT_OrbAtt.append(MLT[a]); 
        P_OrbAtt.append(Pitch[a]); Lat_OrbAtt.append(Lat[a]); Long_OrbAtt.append(Long[a]); R_OrbAtt.append(Radius[a])

    t_OrbAtt = np.array(t_OrbAtt); LS_OrbAtt = np.array(LS_OrbAtt); MLT_OrbAtt = np.array(MLT_OrbAtt);
    P_OrbAtt = np.array(P_OrbAtt); Lat_OrbAtt = np.array(Lat_OrbAtt); Long_OrbAtt = np.array(Long_OrbAtt);
    R_OrbAtt = np.array(R_OrbAtt)
    
    return t_OrbAtt, LS_OrbAtt, MLT_OrbAtt, P_OrbAtt, Lat_OrbAtt, Long_OrbAtt, R_OrbAtt, t_electrons, r1, r2, r3, r4



# Code to make a dial plot with the sun facing up.
import numpy as np
import matplotlib.pyplot as plt

class Dial:
    def __init__(self, ax, angular_bins, radial_bins, H):
        """ 
        This class makes a dial (polar) plot were MLT is the azimuthal
        coordinate and L shell is the radial coordinate. 
        """
        self.ax = ax
        self.angular_bins = angular_bins
        self.radial_bins = radial_bins
        self.H = H

        if 'Polar' not in str(type(ax)):
            raise ValueError('Subplot is not polar. For example, '
                'create ax with \n ax[0] = plt.subplot(121, projection="polar")')
        return

    def draw_dial(self, cb_label, colorbar=True, L_labels=[2,4,6,8], cb_ticksize=15, cb_labelsize=15, mesh_kwargs={}, colorbar_kwargs={}):
        """
        Draws a dial plot on the self.ax subplot object (must have projection='polar' kwarg). 
        colorbar=True - Plot the colorbar or not.
        L_labels=[2,4,6,8] - What L labels to plot
        mesh_kwargs={} - The dictionary of kwargs passed to plt.pcolormesh() 
        colorbar_kwargs={} - The dictionary of kwargs passed into plt.colorbar()
        """
        self.L_labels = L_labels

        angular_grid, radial_grid = np.meshgrid(self.angular_bins, self.radial_bins)

        # Try-except block deals with the dimensions of the mesh and taransposes it
        # if necessary.
        try:
            p = self.ax.pcolormesh(angular_grid*np.pi/12, radial_grid, self.H.T, **mesh_kwargs)
        except TypeError as err:
            if 'Dimensions of C' in str(err):
                p = self.ax.pcolormesh(angular_grid*np.pi/12, radial_grid, self.H, **mesh_kwargs)
            else:
                raise

        self.draw_earth()
        self._plot_params()

        if colorbar:
            cb = plt.colorbar(p, ax=self.ax, **colorbar_kwargs)
            cb.set_label(label=cb_label, size=cb_labelsize)
            cb.ax.tick_params(labelsize=cb_ticksize)
        return

    def draw_earth(self, earth_resolution=50):
        """ Given a subplot object, draws the Earth with a shadow"""
        # Just x,y coords for a line (to map to polar coords)
        earth_circ = (np.linspace(0, 2*np.pi, earth_resolution), np.ones(earth_resolution)) 
        # x, y_lower, y_upper coords for Earth's shadow (maps to polar).
        earth_shadow = (
                        np.linspace(-np.pi/2, np.pi/2, earth_resolution), 
                        0, 
                        np.ones(earth_resolution)
                        )
        self.ax.plot(*earth_circ, c='k')
        self.ax.fill_between(*earth_shadow, color='k')
        return

    def _plot_params(self):
        # Draw L shell contours and get L and MLT labels 
        L_labels_names = self._draw_L_contours()
        mlt_labels = (self.ax.get_xticks()*12/np.pi).astype(int)
        # Sun facing up.
        self.ax.set_xlabel('MLT')
        self.ax.set_theta_zero_location("S") # Midnight at bottom
        self.ax.set_xticklabels(mlt_labels) # Transform back from 0->2pi to 0->24.
        self.ax.set_yticks(self.L_labels)
        self.ax.set_yticklabels(L_labels_names, fontdict={'horizontalalignment':'right'})
        return

    def _draw_L_contours(self, earth_resolution=50):
        """ Plots a subset of the L shell contours. """
        # Draw azimuthal lines for a subset of L shells.
        L_labels_names = [str(i) for i in self.L_labels[:-1]] + [f'L = {self.L_labels[-1]}']
        # L_labels_names = [str(i) for i in self.L_labels]
        for L in self.L_labels:
            self.ax.plot(np.linspace(0, 2*np.pi, earth_resolution), 
                        L*np.ones(earth_resolution), ls=':', c='k')
        return L_labels_names


# define pcolormesh function to handle nan values
def pcolormesh_nan(x: np.ndarray, y: np.ndarray, c: np.ndarray, 
                    ax, projection, cmap=None, norm=None, zorder=1, alpha=1):
    """handles NaN in x and y by smearing last valid value in column or row out,
    which doesn't affect plot because "c" will be masked too
    Stolen from:
    https://github.com/scivision/python-matlab-examples/blob/0dd8129bda8f0ec2c46dae734d8e43628346388c/PlotPcolor/pcolormesh_NaN.py
    
    Personal Notes: Added zorder as an arguement.
    """

    mask = np.isfinite(x) & np.isfinite(y)
    top = None
    bottom = None

    for i, m in enumerate(mask):
        good = m.nonzero()[0]

        if good.size == 0:
            continue
        elif top is None:
            top = i
        else:
            bottom = i

        x[i, good[-1] :] = x[i, good[-1]]
        y[i, good[-1] :] = y[i, good[-1]]

        x[i, : good[0]] = x[i, good[0]]
        y[i, : good[0]] = y[i, good[0]]

    x[:top, :] = np.nanmax(x[top, :])
    y[:top, :] = np.nanmax(y[top, :])

    x[bottom:, :] = np.nanmax(x[bottom, :])
    y[bottom:, :] = np.nanmax(y[bottom, :])

    ax.pcolormesh(x, y, np.ma.masked_where(~mask[:-1, :-1], c)[::-1, ::-1], 
                cmap=cmap, shading='flat', transform=ccrs.PlateCarree(), 
                norm=norm, zorder=zorder, alpha=alpha)
    return