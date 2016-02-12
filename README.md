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

Note that certain information (primarily directory paths) must be input into main.m and
the .ncl scripts in order for them to run.

In addition, the script reads from .mat files that contain GHCND data (rather than the 
.txt files which are provided at http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/.
If you would like access to these .mat files, or information about how to produce them,
please contact Karen McKinnon.

We have also provided the time series of T95, which is produced from the daily station 
data, as part of the supplementary information.

The code 'main.m' serves as a wrapper for the remainder of the subroutines.
