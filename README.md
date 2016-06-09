# PEP
Matlab and ncl codes associated with the paper 
"Long-lead predictions of Eastern US hot days from Pacific sea surface temperatures"
by Karen A. McKinnon, Andrew Rhines, Martin P. Tingley, and Peter Huybers

For inquiries about code, please contact Karen McKinnon at mckinnon AT ucar DOT edu

Most codes are .m files, to be run in Matlab
There are two .ncl files, which run in NCL

While we have made our best attempts at commenting code and providing all dependencies,
it is unlikely that we have been 100% successful. If you are trying to use the code to
reproduce the results and have identified missing files or dependencies, please contact 
Karen McKinnon.

Some additional scripts not produced by the authors, but freely available online, are 
used but not provided in this repository. These include:
m_map
cdfdate2num
runmean
lbmap
rgb

Note that certain information (primarily directory paths) must be input into main.m and
the .ncl scripts in order for them to run.

In addition, the script reads from .mat files that contain GHCND data, rather than the 
.dly files which are provided at http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ as a 
compressed .tar file (ghcnd_all.tar.gz). The .dly files are plain text. We have provided
two .dly files and our correspondingly processed .mat files so that the user can transform
the .dly files to .mat files in similar forms if desired. If you would like further 
information about how to process the .dly files in Matlab, please contact Karen McKinnon.

We have also provided the time series of T95 (PEP-T95TimeSeries.txt), which is produced 
from the daily station data. The .txt file can be found in the supplementary information,
as well as in this github repo. Note that, in the codes, the .txt file is never called, 
but you can use it in order to skip over some of the data processing steps. Also note that
T95 is *technically* defined for all days of the year (and so the values are provided at
daily resolution); however, all of our analysis only relies on the peak summer values.

The code 'main.m' serves as a wrapper for the remainder of the subroutines.
