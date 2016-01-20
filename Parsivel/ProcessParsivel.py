import numpy as np
import pdb
import datetime
import scipy.stats.mstats

def process_parsivel(raw_parsivel_object,time_interval=2,remove_bins=None):

    processed_object = ProcessParsivel(raw_parsivel_object)

    processed_object.apply_matrix(remove_bins=remove_bins)

    #processed_object.plot_diam_fspd()

    processed_object.time_averaging(time_interval=time_interval)

    pdb.set_trace()

    return processed_parsivel


class ProcessParsivel(object):

    '''
    ProcessParsivel takes a RawParsivel object and applies various processing
    to remove erroneous drops and time-average the data

    It returns a processed Parsivel object which includes the corrected matrix,
    error_roce, temperature, and weather code averaged to a new time interval

    It does NOT calculated the DSD. The ProcessParsivel instance needs to be fed
    into the ParsivelDSD class to get the DSD, reflectivity, LWC, etc. 

    Based on Parsivel wxcode, it applies either the standard V vs D conditional
    matrix, or a "frozen" matrix that uses a flat 6 m/s cutoff while still
    including the possibility of mixed rain. I use 6 m/s instead of 4 m/s because
    2 DVD observations indicate fall velocities > 4 m/s in periods of mixed precip
    '''

    def __init__(self,raw_parsivel):
        self.raw_parsivel = raw_parsivel #raw parsivel object (input)
        self.processed_matrix = []
        self.time = []
        self.error_code = []
        self.temperature = []
        self.wxcode = []
        self.matrix = []
        
        self.ndrops_10s = [] # processed numbers of drops (10s)

    def apply_matrix(self,remove_bins = None):
        '''
        Following IDL code written by Patrick Gatlin and Ali Tokay
        Assumes data has correct number of elements

        remove_bins
        2-element array of first and last bin to remove
        None: Remove the 2-smallest (0.064, 0.193 mm) bin per Parsivel default
        [0,2]: Removes the 3rd bin as well (useful if error codes are present)
        [3,5]: Remove bins 3-5 (as well as the 0,1 bin)
        Might be useful for looking at the contribution to rainrate from
        small, medium, or large drops
        '''      
        dmm  = self.drop_diameter # drop size bin (parsivel_diameter)
        delta = self.drop_spread # width of drop size bin
        velt = self.v_theoretical # theoretical terminal velocity
        velo =  self.v # measured velocity (corrected)
        vel1 =  self.vel1 # 50% > Velt
        vel2 = self.vel2 # 50% < Velt

        for i,wx in enumerate(self.raw_parsivel.wxcode):

            # apply conditional matrix to filter out questionable drops
            # rain if wxcode < 66, frozen if wxcode > 65
            if wx < 66:
                matrix_tmp = np.array(self.liquid_matrix)
            else:
                matrix_tmp = np.array(self.frozen_matrix)

            #get array removal elemants using the remove_bins parameter
            if remove_bins != None:
                try:
                    ri = np.array(remove_bins)*32
                    ri[1] +=1 
                    matrix_tmp[ri[0]:ri[1]] = 0 
                except InputError:
                    print('remove_bins must be a 2-element list')
        
            #apply conditional matrix with remove_bins
            matrix_tmp = self.raw_parsivel.matrix[i,:] * matrix_tmp
            if len(self.processed_matrix) == 0:
                self.processed_matrix.append(matrix_tmp)
            else:
                self.processed_matrix = np.vstack((self.processed_matrix,matrix_tmp))
            self.ndrops_10s.append((np.sum(self.raw_parsivel.matrix[i,:])))
   
    def time_averaging(self, time_interval=1, remove_missing=True):
        '''
        Takes processed Parsivel data in default 10s interval
        Time-averages data to the desired time_interval
        Parameters:
        time_interval = 1
        Float that represents the time (in minutes) to average the data
        Note that the data interval is skipped if there is any missing data within
        the time interval--for instance, 1-min averaging requires six 10-s telegrams
        Allowed values are 0.5 (30s) and any integer (1=1 min, 2= 2 min, etc)
        Recommended values are 1,2,5,10,20,30 in order to ensure even time bins

        remove_missing = True
        Flag to remove time periods with missing data
        False will leave time intervals with partial data and ignore empty intervals
        True will fill in missing periods with float('nan')
        '''

        pt = self.raw_parsivel.pytime
        if time_interval >= 1: #averaging to minutes
            for k,t in enumerate(pt):
                #seconds always 0, for minute, subtract off modulus of time interval
                #pdb.set_trace()
                pt[k] = t.replace(second=0)
                min_diff = datetime.timedelta(minutes=t.minute % time_interval)
                pt[k] -= min_diff
        if time_interval < 1: #averaging to 30 seconds
            for k,t in enumerate(pt):
                #subtract off modulus of seconds
                sec_diff = datetime.timedelta(seconds=t.second % (time_interval*60))
                pt[k] -= sec_diff
                        
        #iterate through averaged data, calculate new parameters
        nowtime = pt[0]
        dt = datetime.timedelta(minutes=time_interval)
        while nowtime <= pt[-1]:
            ind = np.where((pt == nowtime))[0]
            matrix_temp = np.zeros((1024)) 
            if len(ind) == time_interval*6: #need exact # of 10s data, otherwise assumptions will be wront
                self.time.append(nowtime)
                #error code is the maximum value in each data chunk
                self.error_code.append(max(self.raw_parsivel.error_code[ind]))
                self.temperature.append(np.mean(self.raw_parsivel.temperature[ind]))
                #wxcode is the mode in each data chunk
                wxcode_tmp = scipy.stats.mstats.mode(self.raw_parsivel.wxcode[ind])[0][0]
                self.wxcode.append(wxcode_tmp)
                for index in ind:
                    matrix_temp += self.raw_parsivel.matrix[index,:]
            else: #if missing data
                self.time.append(nowtime)
                self.error_code.append(float('nan'))
		self.wxcode.append(float('nan'))
		self.temperature.append(float('nan'))
                matrix_temp[:] = float('nan')
            if len(self.matrix) == 0:
                self.matrix.append(matrix_temp)
            else:
                self.matrix = np.vstack((self.matrix,matrix_temp))
          
                
            nowtime += dt  #jump to next time step   

    def plot_diam_fspd(self):
        # Make a 2D histogram of D vs V comparing raw and processed data
        # Created to test if the liquid vs. frozen filters are reasonable
        # and if the remove_bins parameter is working properly

        import matplotlib
        import matplotlib.pyplot as plt

        pm = self.processed_matrix
        rm = self.raw_parsivel.matrix
        pm_sum = np.zeros((32,32))
        rm_sum = np.zeros((32,32))

        # Iterate through samples and sum elements into 32x32 arrays
        for j,row in enumerate(pm[:,0]):
            pm_sum += np.reshape(pm[j,:],(32,32))
            rm_sum += np.reshape(rm[j,:],(32,32))
        
        #mask zeros for plot
        dsd_arr = np.ma.masked_where(pm_sum==0,pm_sum) 
        #dsd_arr = np.ma.masked_where(pm_sum==0,pm_sum) 

        #Atlas et al 1973 theoretical curve
        eqdiam_arr = np.linspace(0,15,1501)
        fspd_arr = 9.65-10.3*np.exp(-0.6*eqdiam_arr)
        fspd_arr[fspd_arr < 0] = float('nan')

        # Plot 2D histogram using pcolor
        fig, ax = plt.subplots()
        fig.suptitle('Heavy Snow/aggregates, processed data', fontsize=15)
        fig.set_size_inches(11,6) #x,y size, default is 8,6
        plt.pcolor(self.drop_diameter,self.v,dsd_arr,vmin = 0, vmax = 600)#vmax=np.max(dsd_arr))
        plt.plot(eqdiam_arr,fspd_arr,'-',lw=2.5,color="black")
             
        plt.xlabel('Diameter (mm)')
        plt.ylabel('Velocity (m s$^{-1}$)')
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('# of Particles')
        ax.set_ylim(0,20)
        ax.set_xlim(0,20)
    
        plt.show()
        pdb.set_trace()
        

    def info(self):
        print 'Raw Parsivel: '
        print 'Pytime length: '+str(len(self.pytime))
       

    drop_diameter = [
        0.064, 0.193, 0.321, 0.45, 0.579, 0.708, 0.836, 0.965, 1.094, 1.223, 1.416, 1.674,
        1.931, 2.189, 2.446, 2.832, 3.347, 3.862, 4.378, 4.892, 5.665,
        6.695, 7.725, 8.755, 9.785, 11.330, 13.390, 15.45, 17.51, 19.57, 22.145, 25.235] #also dmm

    drop_diameter_ott = [
        0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625,
        1.875, 2.125, 2.375, 2.750, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 
        13, 15, 17, 19, 21.5, 24.5] #diameters from the OTT Parsivel2 manual

    drop_spread = [
        0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.257,
        0.257, 0.257, 0.257, 0.257, 0.515, 0.515, 0.515, 0.515, 0.515, 1.030, 1.030,
        1.030, 1.030, 1.030, 2.060, 2.060, 2.060, 2.060, 2.060, 3.090, 3.090] #also delta

    v_theoretical = [
        0.089, 0.659, 1.239, 1.803, 2.353, 2.889, 3.404, 3.892,
        4.329, 4.705, 5.217, 5.833, 6.389, 6.886, 7.326, 7.878,
        8.424, 8.785, 9.002, 9.117, 9.173, 9.248, 9.323, 9.398,
        9.473, 9.586, 9.735, 9.885, 10.035, 10.185, 10.372, 10.597] #also velt

    v = [
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9,
        2.2, 2.6, 3, 3.4, 3.8, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
        17.6, 20.8] #also velo

    v_spread = [.1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .2, .2, .2, .2, .2, .4,.4, 
                .4, .4, .4, .8, .8, .8, .8, .8, 1.6, 1.6, 1.6, 1.6, 1.6, 3.2, 3.2] #gets vel1, vel2

    vel1 = [
      0.045, 0.329, 0.6200, 0.901, 1.177, 1.444, 1.702, 1.946, 2.165, 2.352,
      2.608, 2.916,  3.194, 3.443, 3.663, 3.939, 4.212,  4.392, 4.501, 4.559,
      4.587, 4.624,  4.6620, 4.699, 4.737, 4.793, 4.868,  4.943, 5.018,  5.093,
      5.186, 5.299]

    vel2 = [
      0.134, 0.988, 1.859, 2.704, 3.530, 4.333, 5.106, 5.837, 6.494, 7.057,
      7.825, 8.749, 9.583, 10.33, 10.989, 11.816, 12.635, 13.177, 13.503, 13.676,
      13.76, 13.873, 13.985, 14.097, 14.21, 14.378, 14.603, 14.828, 15.053, 15.278,
      15.559, 15.896]

    #matrices (to right = larger diameter, down = faster fall velocity)

    liquid_matrix = [
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    #frozen matrix: velocity cutoff: 6 m/s, size cutoff 12 mm (Ali Tokay, personal communication)
    frozen_matrix = [
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
