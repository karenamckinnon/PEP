load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;*******************************************************
;*** Calculate pressure anomalies from seasonal      ***
;*** cycle and running mean                          ***
;*******************************************************

;-------------------------------------------------------------
; User specifications
;-------------------------------------------------------------

ymdStrt = 19790101                         ; start yyyymmdd
ymdLast = 20150731                         ; last
yrStrt  = ymdStrt/10000
yrLast  = ymdLast/10000

nhar    = 3                                ; number of fourier comp

var     = "uwnd"                            ; name of file

dir    = "DATADIR"							; ** put your data directory here **
fnameAnom = var+".10m.anomalies."+yrStrt+"."+yrLast+".nc"

;***********************************************************
; Read user specified time and create required yyyyddd
;***********************************************************
files = systemfunc("ls "+dir+var+".10m.gauss.*.nc")
f       = addfiles(files, "r")

;***********************************************************
; Read user specified time and create required yyyyddd
;***********************************************************
time    = f[:]->time          ; time:units = "hours since"
TIME    = cd_calendar(time, 0)             ; type float
year    = floattointeger( TIME(:,0) )
month   = floattointeger( TIME(:,1) )
day     = floattointeger( TIME(:,2) )
ddd     = day_of_year(year, month, day)
yyyyddd = year*1000 + ddd                  ; needed for input

;***********************************************************
; Read data: short2flt
;***********************************************************
x       =  short2flt( f[:]->$var$ )    ; convert to float
xnew    = rm_single_dims(x)
delete(x)
x = xnew
delete(xnew)
printVarSummary( x )

;**********************************************************
; Compute daily climatology: raw and then 'smoothed'
;***********************************************************
xLevel = x
xClmDay = clmDayTLL(xLevel, yyyyddd)     ; daily climatology at each grid point

;***********************************************************
; Compute smoothed daily climatology using 'nhar' harmonics
;***********************************************************
xClmDay_sm = smthClmDayTLL(xClmDay, nhar)

;***********************************************************
; Compute daily anomalies using smoothed climatologies
;***********************************************************
xAnom  = calcDayAnomTLL (xLevel, yyyyddd, xClmDay_sm)
xAnom@long_name = "Anomalies from Smooth Daily Climatology"

copy_VarMeta(x, xAnom)
delete(x)

;***********************************************************
; Create netCDF
;***********************************************************

system("/bin/rm -f "+dir+fnameAnom)      ; rm any pre-exist file, if any
ncdf    = addfile (dir+fnameAnom, "c")
; make time and UNLIMITED dimension      ; recommended  for most applications
filedimdef(ncdf,"time",-1,True)
ncdf->uwnd = xAnom

end