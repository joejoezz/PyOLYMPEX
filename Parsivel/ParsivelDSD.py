import numpy as np
import pdb
import datetime
import scipy.stats.mstats
import RawParsivel as rp
import ProcessParsivel as pp
import glob

def calc_dsd(apu,sitename,date,time_interval=1):
    '''
    Funtion that wraps all three Parsivel classes together to make the
    final DSD object

    Use to get the data to make plots
    '''


    indir = '/home/disk/funnel/olympex/archive2/'+apu+'/Parsivel/'+date[0:6]+'/'
    searchfor= indir+apu+'_'+date+'*'     
    infiles = glob.glob(searchfor)

    rpdata = rp.read_parsivel(infiles)
    ppdata = pp.process_parsivel(rpdata,time_interval=time_interval)
    dsd = ParsivelDSD(ppdata)
    dsd.get_precip_params()
    return dsd


class ParsivelDSD(object):

    '''
    ParsivelDSD takes a ProcessParsivel instance and computes the DSD using
    the methods from Tokay et al. (2014). This code is based on the
    get_rainparams.pro code written by Ali Tokay and converted to IDL by 
    Dave Wolff. 

    The reason for keeping the ParsivelDSD as a separate class is that the
    user may be interested in seeing characteristics of the DSD that do
    not appear in the final processed dataset...for instance, the contribution
    to LWC from each Parsivel bin. By using this class, one can branch off from
    various steps on the processing process and make plots or analyze the data.

    In some cases, the data is saved as a tuple, with the contribution to rainrate
    from each drop bin range included as the second element in the tuple. 

    The PlotParsivel class can be used to make plots of this data.

    Notes:
    Currently uses the corrected fall velocities from Ali Tokay...changing the 
    assumed fall velocities will affect the drop parameters. 

    Future additions:
    Standard deviation of drop size
    '''

    def __init__(self,processed_parsivel):
        #some data stored as tuples if applicable
        #first element = average or total across all 32 DSD bins
        #second element = value for each bin individually (plus time-dimension)
        self.proc_p2 = processed_parsivel
        self.time = self.proc_p2.time 
        timedim = len(self.proc_p2.matrix[:,0]) #time dimension
        #ndrops is the number of drops (non-normalzed) in each volume of air
        self.ndrops = (np.zeros(timedim),np.zeros((timedim,32)))
        #DSD is the number of drops per volume of air for a given drop size interval (i.e. for the 32 bins)
        self.dsd = np.zeros((timedim,32)) #dsd 
        #drop_conc is the number of drops per volume of air. It is NOT normalized by drop size interval like DSD
        self.drop_conc = (np.zeros(timedim),np.zeros((timedim,32))) #first tuple is the total drop conc
        self.lwc = (np.zeros(timedim),np.zeros((timedim,32))) #lwc from each drop size bin
        self.z = (np.zeros(timedim),np.zeros((timedim,32))) #reflectivity factor from each drop size bin
        self.dbz = np.zeros(timedim)
        self.rainrate = (np.zeros(timedim),np.zeros((timedim,32))) #rainrate from each drop size bin
        self.rainaccum = (np.zeros(timedim),np.zeros((timedim,32))) #total rain accumulation from each bin
        self.dmax = np.zeros(timedim) #max drop size
        self.dm = np.zeros(timedim) #mass-weighted mean diameter
        self.sigma_m = np.zeros(timedim) #variance of mass spectrum
        self.moments = np.zeros((timedim,8)) #drop moments
        #self.rainrate_old = (np.zeros(timedim),np.zeros((timedim,32))) #testing
    

    def get_precip_params(self):
        
        #Takes a processed_parsivel instance and calculates the precip params
         
        # parsivel_matrix: 1ength 1024 matrix of Parsivel diam/fspd
        # timerain: time of period in seconds #self.proc_p2.time_interval
        # wxcode: needed to differentiate between rain and snow
        # Note: if frozen precipitation detected, no rain rate or LWC is returned
        
        #add loop here to go through each time dimension
        #timerain has to be calculated individually for each record in case data is missing 
        #what usually happens is that maybe 1 min of data per day is missing in 10-30s intervals
        timerain = self.proc_p2.time_interval
        nrecords_exp = timerain*6
        nrecords_actual = np.array(self.proc_p2.num_records)
        nrecords_missing = (nrecords_exp - nrecords_actual).astype(float) #missing records
        for td,dmax in enumerate(self.dmax):
            #get correct time multiplier
            time_mult = 60 * timerain - nrecords_missing[td]*10 #units: seconds
            time_div = 60 / (timerain - (nrecords_missing[td]/6)) #units: s/min
            #reshape to 32x32
            matrix = np.reshape(self.proc_p2.matrix[td,:],(32,32))
            ndrops = np.sum(matrix) #total drop count for this time step **need to output this** as self.ndrops
            #moments: 0) concen, 1) mean diam, 2) surface area conc, 3)lwc, 6)z
            x4 = 0.0 #for moments
            x3 = 0.0 #for moments
            #Process each drop bin (diameter, velocity):
            for dind,dbin in enumerate(self.drop_diameter):
                for vind,vbin in enumerate(self.v_idlcode):
                    drops = matrix[vind,dind]
                    #Next step uses equation (6) from Tokay et al. (2014)
                    #parsivel laser area, units: mm^2 (subtracting off partial drops on edge)
                    p2_area2 = 180.*(30.-(dbin/2.))/100. #not sure why we divide by 100
                    p2_area = 180.*(30.-(dbin/2.))
                    #denominators
                    denom2 = time_mult * p2_area2 * vbin * self.drop_spread[dind] * 100 #per m^3*mmbin
                    denom2_beard = time_mult * p2_area2 * self.v_theoretical[vind] * self.drop_spread[dind] * 100 
                    self.dsd[td,dind] += (1.e6 * drops)/denom2 #10^6 converts to m^3*mm instead of mm^3*m
                    if self.proc_p2.wxcode[td] < 65: #use theoretical fall speed for rain
                        #units: s*mm^2*m/s*mm = mm^3*m
                        denominator = time_mult * p2_area * vbin #per m^3 (not per bin)
                    else: #no change for snow as of right now...could add in later
                        denominator = time_mult * p2_area * vbin 

                    vol = np.pi*dbin**3/6 #volume of 1 drop in this size bin
                    if drops > 0: self.dmax[td] = dbin #make dmax for this bin nonzero if >0 drops
                    self.drop_conc[1][td,dind] += drops*1.e6/denominator #direct evaulation of Tokay eq 6 for drop conc
                    self.lwc[1][td,dind] += drops*vol*1.e3/denominator #units: g/m^3 per mm bin (rho=1000 g/m^3)
                    self.z[1][td,dind] += drops * 1.e6 * dbin**6/denominator #reflectivity factor (6th power of diameter, normalized for area, time)
                    self.rainrate[1][td,dind] += drops * vol / p2_area *time_div #rainrate
                    #self.rainrate[1][td,dind] += drops * vol / p2_area * 60 / timerain #rainrate old
                    #4th and 3rd moments
                    if drops > 0:
                        x4 += 1.e6 * drops * dbin**4 / denominator
                        x3 += 1.e6 * drops * dbin**3 / denominator
                        for ind,moment in enumerate(self.moments[td,:]):
                            self.moments[td,ind] += 1.e6 * drops * dbin**ind / denominator
                    
            #compute total drop_conc, average lwc, etc (1st dimension of tuples)
            self.drop_conc[0][td] = np.sum(self.drop_conc[1][td,:]) #this is good
            self.lwc[0][td] = np.sum(self.lwc[1][td,:])
            self.z[0][td] = np.sum(self.z[1][td,:])
            self.rainrate[0][td] = np.sum(self.rainrate[1][td,:])
            #mass-weighted mean diameter = ratio of 4th to 3rd moments
            self.dm[td] = self.moments[td,4] / self.moments[td,3] 
            if self.z[0][td] > 0:
                self.dbz[td] = 10 * np.log10(self.z[0][td]) #z to dbz
            else:
                self.dbz[td] = float('nan')

            #eventually need to add sigma_m here


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

    #note: drop concentrations are affected by which velocity table is used
    v_theoretical = [
        0.089, 0.659, 1.239, 1.803, 2.353, 2.889, 3.404, 3.892,
        4.329, 4.705, 5.217, 5.833, 6.389, 6.886, 7.326, 7.878,
        8.424, 8.785, 9.002, 9.117, 9.173, 9.248, 9.323, 9.398,
        9.473, 9.586, 9.735, 9.885, 10.035, 10.185, 10.372, 10.597] #also velt

    v_parsivel = [
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9,
        2.2, 2.6, 3, 3.4, 3.8, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
        17.6, 20.8]

    v_idlcode = [
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.854, 0.962, 1.128, 1.354, 1.588, 1.828, 2.075,
        2.398, 2.782, 3.15, 3.502, 3.838, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
        17.6, 20.8] #also velo

