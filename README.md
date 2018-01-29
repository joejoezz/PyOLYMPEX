PyOLYMPEX is a Python module to handle ground data from the OLYMPEX field campaign.

Required modules:
numpy, scipy

Features
Parsivel
- Processes Raw Parsivel APU (.dat) files from the OLYMPEX field campaign. 
- Allows user to specify time interval (dt) for computing N(D)
- Allows user to specify maximum allowable drop size
- Identifies and removes periods of snow contamination
- Identifies and removes periods with instrument error codes 
- Easily allows for slicing data by drop size bin
Iowa Rain Gauges
- Processes Iowa Gauge packets from the Iowa Gauge server
- Writes daily text and metadata files
- Computes rain rates at user-specified time intervals (e.g 30 min, 1 hour, 1 day). 

Instructions:

Parsivel:

```
import ParsivelDSD as pdsd
```

Download raw Parsivel data at: ftp://trmm-fc.gsfc.nasa.gov/distro/apu/

```
apu = 'apu06'
sitename = 'Prairie Creek'
date = '20160127'
time_interval = 30' 
```

Compute DSD and derived parameters
```dsd = pdsd.calc_dsd(indir,apu,sitename,date,time_interval=time_interval)```

Iowa Gauges:

import methods
```
import IowaGaugeRaw as igr
```
Specify list of gauge numbers, date, and directories to save in
```
gaugelist = ['NASA0043','NASA0028']
yyyymmdd = '20151101'
yyyymmdd2 = '20151102' #specify following date unless using for real-time
outdir = 'directory for gauge files'
outdir_meta = 'directory for metadata files'
for gauge in gaugelist:
    igr.save_iowa_gauge(yyyymmdd,gauge,outdir,outdir_meta,nextdate=yyyymmdd2)
```





