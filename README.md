PyOLYMPEX is a Python module to handle data from the OLYMPEX field campaign

Required modules:
numpy, scipy

Current features: 
-Processes Raw Parsivel APU (.dat) files from the OLYMPEX field campaign. 
-Allows user to specify averating time (usually 1 min but can be any multiple of 10 seconds)
-Different filtering for rain vs. snow
-Removal of any range of drop size bins

Future Features:
Priority #1: Parsivel
 -ParsivelDSD class that uses Tokay et al. (2014) method to compute DSD, LWC, Z, etc. (under construction)
 -Write data to text files
 -Plot timeseries, drop size distributions, and fall velocities

Other planned additions
2DVD Processing (to compare with Parsivel)
Integrate PyART to plot NPOL, DOW, and EC X-band radars
Plot GPM overpasses
Plot MRR
Plot PIP disdrometer data
Process and plot rain gauge data (both Iowa and standalone gauges)

Instructions for use:

import methods
<import RawParsivel as rp>
<import ProcessParsivel as pp>

Download raw Parsivel data at: ftp://trmm-fc.gsfc.nasa.gov/distro/apu/
<infile = '/dir/apu01_2015110100.dat'>

Create Raw Parsivel instance (by reading data):
<rpdata = rp.read_parsivel(infile)>

Create Processed Parsivel instance
<ppdata = pp.process_parsivel(rpdata,time_interval=2)>

Make a test plot of diameter vs. fall speed
<ppdata.plot_diam_fspd()>


