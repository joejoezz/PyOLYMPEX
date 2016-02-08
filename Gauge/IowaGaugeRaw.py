'''
Program to read and plot Iowa tipping bucket gauge data
Created by Joe Zagrodnik 5-Feb-2016
Read data from Iowa server
Write tipping bucket data to daily text files
Write metadata files as well
My Iowa Gauge class is designed to read Iowa packets. The main purpose
of this code is to get the data off of the Iowa server and put it into a format
that is easily readable for other purposes.
The write_text_files() function will write in the exact same format as Dave W's IDL code

This code also handles 2 issues with Iowa gauge packets:
1) Duplicate tips
2) Tips from the last 15 minutes of the day included in the first packet
   from the following day

To handle the first issue, I include an optional 'nextdate' parameter where the
user can input the next day (yyyymmdd string). The extra tips from days that are 
not 'date' are removed, as well as duplicates. 
'''
import os
import pdb
import numpy as np
import datetime

def save_iowa_gauge(date,gauge,outdir,outdir_meta,nextdate = None):
    
    tmpdir = '/home/disk/funnel/olympex/archive2/wxdata/Gauge/tmpdir/'

    rawgauge = IowaGaugeRaw(gauge,date)

    try:
        rawgauge.import_packets(tmpdir,gauge,nextdate=nextdate)
    except:
        print 'No packets found for '+rawgauge.gauge+' '+rawgauge.date
        return

    rawgauge.read_packets(tmpdir)

    #delete all packets
    rawgauge.delete_packets(tmpdir)

    #clean data
    rawgauge.clean_data()

    #write text files
    rawgauge.write_text_files(outdir,outdir_meta)


class IowaGaugeRaw(object):
    '''
    Raw Iowa Gauge data
    '''

    def __init__(self,gauge,date):
        self.gauge = gauge
        self.date = date
        self.lat = 0.0
        self.lon = 0.0
        self.time_a = []
        self.time_b = []
        self.metatime = []
        self.voltage = []
        self.temperature = []
        self.wetness = []
        self.solar = []
        self.rssi = []
        self.rain_a = []
        self.rain_b = []
        
    
    def import_packets(self,tmpdir,gauge,nextdate=None):
        #import packets to a temporary directory 
        from bs4 import BeautifulSoup
        import urllib2
        import urllib

        if nextdate == None:
            dates = [self.date]
        else:
            dates = [self.date,nextdate]

        for d in dates:
            path = 'http://s-iihr61.iihr.uiowa.edu/sensors/'+gauge+'/'+d[0:4]+'/'+d[4:6]+'/'+d[6:8]+'/'
            page1 = urllib2.urlopen(path)
            soup = BeautifulSoup(page1)

            for tag in soup.findAll('a'):
                filename = tag.get('href')
                if filename[0:4] == "NASA":
                    urllib.urlretrieve(path+filename, \
                       filename = tmpdir+filename+'.txt')

                    #print 'Got '+filename

                             
    def read_packets(self,tmpdir):
        import glob
        searchfor= tmpdir+'NASA*'
        #read each text file, append variables
        for infile in sorted(glob.glob(searchfor)):
            with open(infile, 'r') as f:
                for i,line in enumerate(f):
                    if i == 1: #only grab metatime for current day
                        year = int(line.split()[0].split('-')[0])
                        mon = int(line.split()[0].split('-')[1])
                        day = int(line.split()[0].split('-')[2])
                        hr = int(line.split()[1].split(':')[0])
                        minute = int(line.split()[1].split(':')[1])
                        sec = int(line.split()[1].split(':')[2])
                        if day == int(self.date[6:8]):
                            self.metatime.append(datetime.datetime(year,mon,day,hr,minute,sec))
                    elif i == 2:
                        self.lat = float(line.split(',')[0]) #should stay the same
                        self.lon = float(line.split(',')[1]) 
                    elif i == 3 and day == int(self.date[6:8]): #only grab for current day
                        metadata = line.split(',')
                        try:
                            self.voltage.append(float(metadata[0]))
                            self.temperature.append(float(metadata[1]))
                            self.wetness.append(float(metadata[2]))
                            self.solar.append(float(metadata[3]))
                            self.rssi.append(float(metadata[4])) #has to be converted?
                        except: #sometimes bad data...
                            self.voltage.append(-99)
                            self.temperature.append(-99)
                            self.wetness.append(-99)
                            self.solar.append(-99)
                            self.rssi.append(-99)
                    elif i > 4:
                        td = line.split(',')[0]
                        year = int(td.split()[0].split('-')[0])
                        mon = int(td.split()[0].split('-')[1])
                        day = int(td.split()[0].split('-')[2])
                        hr = int(td.split()[1].split(':')[0])
                        minute = int(td.split()[1].split(':')[1])
                        sec = int(td.split()[1].split(':')[2])
                        if int(line.split(',')[1]) == 81: #gauge a
                            self.time_a.append(datetime.datetime(year,mon,day,hr,minute,sec))
                            self.rain_a.append(0.254)
                        if int(line.split(',')[1]) == 82: #gauge b
                            self.time_b.append(datetime.datetime(year,mon,day,hr,minute,sec))
                            self.rain_b.append(0.254)
                    else:
                        continue

    def clean_data(self):
        #remove duplicates
        self.time_a = sorted(list(set(self.time_a)))
        self.rain_a = np.array(self.rain_a[0:len(self.time_a)])
        self.time_b = sorted(list(set(self.time_b)))
        self.rain_b = np.array(self.rain_b[0:len(self.time_b)])

        #extract the day in order to filter the arrays
        day_a = []
        day_b = []
        for tt in self.time_a:
            day_a.append(tt.day)
        for tt in self.time_b:
            day_b.append(tt.day)

        self.time_a = np.array(self.time_a)
        self.time_b = np.array(self.time_b)
        day_a = np.array(day_a)
        day_b = np.array(day_b)

        #remove data from wrong day
        good_a = np.where((day_a == int(self.date[6:8])))[0]
        good_b = np.where((day_b == int(self.date[6:8])))[0]

        self.time_a = self.time_a[good_a]
        self.time_b = self.time_b[good_b]
        self.rain_a = self.rain_a[good_a]
        self.rain_b = self.rain_b[good_b]

        #print 'cleaned data'

    def delete_packets(self,tmpdir):
        #empties the temporary directory
        os.system('rm -rf '+tmpdir+'NASA*')
        #print 'deleted all tmp packets'

    def write_text_files(self,outdir,outdir_meta):
        #writes meta and gauge text files
        #write in same format as NASA files
        time_a_out = []
        time_a_out.append(self.gauge+' '+str(self.lat)+' '+str(self.lon))
        time_b_out = []
        time_b_out.append(self.gauge+' '+str(self.lat)+' '+str(self.lon))
        for t in self.time_a:
            time_a_out.append(t.strftime("%Y-%m-%d %H:%M:%S"))
        for t in self.time_b:
            time_b_out.append(t.strftime("%Y-%m-%d %H:%M:%S"))

        time_a_out = np.array(map(str,time_a_out))
        time_b_out = np.array(map(str,time_b_out))
        title_a = self.gauge+'_A_'+self.date[0:4]+'-'+self.date[4:6]+'-'+self.date[6:8]+'.txt'
        title_b = self.gauge+'_B_'+self.date[0:4]+'-'+self.date[4:6]+'-'+self.date[6:8]+'.txt'

        #save gauge files
        os.chdir(outdir)
        os.system('mkdir -p ' + self.date[0:6])
        outdir2 = outdir+self.date[0:6]+'/'
        os.chdir(outdir2)
        os.system('mkdir -p ' + self.date)
        np.savetxt(outdir2+self.date+'/'+title_a,time_a_out,delimiter=' ',fmt="%s")
        np.savetxt(outdir2+self.date+'/'+title_b,time_b_out,delimiter=' ',fmt="%s")
        print 'saved '+title_a
        print 'saved '+title_b

        #save metadata
        metatime_out = np.array(map(str,self.metatime))
        title_meta = self.gauge+'_'+self.date[0:4]+'_'+self.date[4:8]+'_meta.txt'
        metastr = []
        for i,m in enumerate(metatime_out):
            metastr.append(self.gauge+' '+m+' '+str(self.lat)+' '+str(self.lon)+' '+str(self.voltage[i])+','+
                           str(self.temperature[i])+','+str(self.wetness[i])+','+str(self.solar[i])+','+str(self.rssi[i]))
        metastr_out = np.array(map(str,metastr))

        os.chdir(outdir_meta)
        os.system('mkdir -p ' + self.date[0:6])
        outdir_meta2 = outdir_meta+self.date[0:6]+'/'
        os.chdir(outdir_meta2)
        os.system('mkdir -p ' + self.date)     
        np.savetxt(outdir_meta2+self.date+'/'+title_meta,metastr_out,delimiter=' ',fmt="%s")
        print 'saved '+title_meta


