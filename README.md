###PyOLYMPEX is a Python module to handle data from the OLYMPEX field campaign

####Required modules:
numpy, scipy

####Current features: 
- Processes Raw Parsivel APU (.dat) files from the OLYMPEX field campaign. 
- Allows user to specify averating time (usually 1 min but can be any multiple of 10 seconds)
- Different filtering for rain vs. snow
- Easily allows for slicing data by drop size bin

####Future Features:
#####Priority #1: Parsivel
- Write data to text files
- Plot timeseries, drop size distributions, and fall velocities

#####Other planned additions
- 2DVD Processing (to compare with Parsivel)
- Integrate PyART to plot NPOL, DOW, and EC X-band radars
- Plot GPM overpasses
- Plot MRR
- Plot PIP disdrometer data
- Process and plot rain gauge data (both Iowa and standalone gauges)

####Instructions for use:

######import methods
```
import ParsivelDSD as pdsd
```

######Download raw Parsivel data at: ftp://trmm-fc.gsfc.nasa.gov/distro/apu/

######Specify directory and site name
```
apu = 'apu06'
sitename = 'Prairie Creek'
date = '20160127'
time_interval = 30' 
indir = [user should specify indir]
```

######Compute DSD and derived parameters
```dsd = pdsd.calc_dsd(indir,apu,sitename,date,time_interval=time_interval)```



