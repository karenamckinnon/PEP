load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin 

;*******************************************************
;*** Using the weighted composites of uwnd and vwnd  ***
;*** calculate the composite Plumb fluxes using the  ***
;*** equations from Takaya and Nakamura (1997)       ***
;*******************************************************

; At command line, provide lagMin, lagMax, and hotdays filename

; Read in days that are to be composited on
cacheDir = "" ; ** put your cache directory here! **
dataDir = "" ; ** put your data directory here! **
hotDaysfname = cacheDir+hotDaysFname
hotDays = asciiread(hotDaysfname,-1,"double")

; Load U and V
uwnd_files = systemfunc("ls " + dataDir + "/NCEP_R2/uwnd/uwnd.1979.2015.nc")
vwnd_files = systemfunc("ls " + dataDir + "/NCEP_R2/vwnd/vwnd.1979.2015.nc")
hgt_files = systemfunc("ls " + dataDir + "/NCEP_R2/hgt/hgt.*.nc")
temp_files = systemfunc("ls " + dataDir + "/NCEP_R2/air/air.*.nc")
fuwnd = addfiles(uwnd_files(0),"r")
fvwnd = addfiles(vwnd_files(0),"r")
fhgt = addfiles(hgt_files,"r")
ftemp = addfiles(temp_files,"r")

; Get dimensions (same for uwnd and vwnd)
lat = fuwnd[0]->lat
lon = fuwnd[0]->lon
time = fuwnd[:]->time
level = fuwnd[:]->level

lat = lat(::-1) ; switch to S-->N
uwnd = short2flt(fuwnd[0]->uwnd(:,:,::-1,:))
vwnd = short2flt(fvwnd[0]->vwnd(:,:,::-1,:))
hgt = short2flt(fhgt[:]->hgt(:,:,::-1,:))
temp = short2flt(ftemp[:]->air(:,:,::-1,:))

print("Got data")

; Gas constant
gc = 290
; Gravitational acceleration
ga = 9.80665
; Radius of the earth
re = 6378388
; scale height
sclhgt=8000.
; pi
pi = atan(1.0)*4.
; Coriolis parameter
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0 = "lat"
f&lat = lat
f@_FillValue = default_fillvalue(typeof(f))

; cosine
coslat = cos(lat(:)*pi/180.)

; GET JJA FOR MEAN STATE
tvec = cd_calendar(time,0)
mmDay = tvec(:,1)
ixDayJJA = ind(((mmDay.eq.6).or.(mmDay.eq.7).or.(mmDay.eq.8)))

U = dim_avg_n_Wrap(uwnd(ixDayJJA,:,:,:),0)
V = dim_avg_n_Wrap(vwnd(ixDayJJA,:,:,:),0)
T = dim_avg_n_Wrap(temp(ixDayJJA,:,:,:),0)

print("Got base state")

; Calculate anomalies from climatology
nlevels = dimsizes(fuwnd[0]->level)
uwndAnom = new(dimsizes(uwnd),typeof(uwnd),uwnd@_FillValue)
vwndAnom = new(dimsizes(vwnd),typeof(vwnd),vwnd@_FillValue)

TIME    = cd_calendar(time, 0)
year    = floattointeger( TIME(:,0) )
month   = floattointeger( TIME(:,1) )
day     = floattointeger( TIME(:,2) )
ddd     = day_of_year(year, month, day)
yyyyddd = year*1000 + ddd
nhar    = 3

do i = 0,(nlevels-1),1
   print(sprintf("%2g",i))
;***********************************************************
; Compute daily climatology: raw and then 'smoothed'
;***********************************************************
   uLevel = uwnd(:,i,:,:)
   vLevel = vwnd(:,i,:,:)
   uClmDay = clmDayTLL(uLevel, yyyyddd)     ; daily climatology at each grid point
   vClmDay = clmDayTLL(vLevel, yyyyddd)

;***********************************************************
; Compute smoothed daily climatology using 'nhar' harmonics
;***********************************************************
   uClmDay_sm = smthClmDayTLL(uClmDay, nhar)
   vClmDay_sm = smthClmDayTLL(vClmDay, nhar)

;***********************************************************
; Compute daily anomalies using smoothed climatologies
;***********************************************************
	uwndAnom(:,i,:,:)   = calcDayAnomTLL (uLevel, yyyyddd, uClmDay_sm)
	vwndAnom(:,i,:,:)   = calcDayAnomTLL (vLevel, yyyyddd, vClmDay_sm)

	delete(uLevel)
	delete(vLevel)
	delete(uClmDay)
	delete(vClmDay)
	delete(uClmDay_sm)
	delete(vClmDay_sm)

end do
copy_VarMeta(uwnd, uwndAnom)
copy_VarMeta(vwnd, vwndAnom)

; Remove trend
uwndAnom2 = dtrend_n(uwndAnom, True, 0)
vwndAnom2 = dtrend_n(vwndAnom, True, 0)

delete(uwndAnom)
delete(vwndAnom)

dsizes = dimsizes(U)
nlevel = dsizes(0)
nlat = dsizes(1)
nlon = dsizes(2)

levelMat = conform_dims(dimsizes(T),level,0)
coslatMat = conform_dims(dimsizes(T),coslat,1)
fMat = conform_dims(dimsizes(T),f,1)

; buoyancy frequency, calculated using climatology
theta = T*(1000./levelMat)^0.286
; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(theta, -sclhgt*log(level/1000), False, 0, 0)

NN = (gc*(levelMat/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = default_fillvalue(typeof(NN))

; Compute |U|
modU = U
modU = (/modU@_FillValue/)
modU = (/sqrt(U^2 + V^2)/)
modU@long_name = "Magnitude of velocity field"
modU@units = "m/s"

; find time indices associated with hot days
idxHot = get1Dindex(time, hotDays)
idxHotShort = idxHot(ind(.not.ismissing(idxHot)))
delete(idxHot)
idxHot = idxHotShort

lags = ispan(lagMin, lagMax, 1)
lags = lags(::-1)
nlags = dimsizes(lags)

Wx = new((/nlags,nlevel,nlat,nlon/),float,T@_FillValue)
Wx = (/Wx@_FillValue/)
Wx@long_name = "Zonal wave activity flux"
Wx@units = "m2/s2"

Wy = new((/nlags,nlevel,nlat,nlon/),float,T@_FillValue)
Wy = (/Wy@_FillValue/)
Wy@long_name = "Meridional wave activity flux"
Wy@units = "m2/s2"

Wz = new((/nlags,nlevel,nlat,nlon/),float,T@_FillValue)
Wz = (/Wz@_FillValue/)
Wz@long_name = "Vertical wave activity flux"
Wz@units = "m2/s2"

; Want to compute WAF for different lead times
lag = lagMax ; start with longest lag
indSave = 0

do while(lag.ge.lagMin)
	print(lag)
	idxUse = idxHot - lag
	; remove negative indices
	idxUseNew = idxUse(ind(idxUse > 0))
	delete(idxUse)
	idxUse = idxUseNew
	delete(idxUseNew)
	uComp = dim_avg_n_Wrap(uwndAnom2(idxUse, :, :, :), 0)
	vComp = dim_avg_n_Wrap(vwndAnom2(idxUse, :, :, :), 0)
	uComp@_FillValue = default_fillvalue(typeof(uComp))
	vComp@_FillValue = default_fillvalue(typeof(vComp))

	sf = uComp
	vp = uComp
	sf = (/sf@_FillValue/)
	vp = (/vp@_FillValue/)
	uv2sfvpf(uComp,vComp,sf,vp)
	delete(vp)
	delete(uComp)
	delete(vComp)

	sf@long_name = "streamfunction"
	sf@units = "m^2/s"

	; Compute all first and second (plus cross) derivatives of the streamfunction
	sf_x = sf
	sf_y = sf
	sf_xx = sf
	sf_yy = sf
	sf_xy = sf
	sf_yx = sf
	sf_x = (/sf_x@_FillValue/)
	sf_y = (/sf_y@_FillValue/)
	sf_xx = (/sf_xx@_FillValue/)
	sf_yy = (/sf_yy@_FillValue/)
	sf_xy = (/sf_xy@_FillValue/)
	sf_yx = (/sf_yx@_FillValue/)
	gradsf(sf,sf_x,sf_y)
	gradsf(sf_x,sf_xx,sf_xy)
	gradsf(sf_y,sf_yx,sf_yy)
	delete(sf_yx) ; Redundant

	; do vertical derivatives
	sf_z = center_finite_diff_n(sf, -sclhgt*log(level/1000),False,0,0)
	sf_xz = center_finite_diff_n(sf_x, -sclhgt*log(level/1000),False,0,0)
	sf_yz = center_finite_diff_n(sf_y, -sclhgt*log(level/1000),False,0,0)

	Wx(indSave, :, :, :) = (/ ( U * (sf_x^2 - sf * sf_xx) + V * (sf_x * sf_y - sf * sf_xy ) ) * (levelMat/1000) / (2 * modU) /)
	Wy(indSave, :, :, :) = (/ ( U * ( sf_x * sf_y - sf * sf_xy ) + V * ( sf_y^2 - sf * sf_yy ) ) * (levelMat/1000) / (2 * modU) /)
	Wz(indSave, :, :, :) = (/ ( U * ( sf_x * sf_z - sf * sf_xz ) + V * ( sf_y * sf_z - sf * sf_yz ) ) * (fMat * fMat / NN) * ((levelMat/1000) / (2 * modU)) /)

	lag = lag - 1
	indSave = indSave + 1

	delete(sf_x)
	delete(sf_y)
	delete(sf_xx)
	delete(sf_xy)
	delete(sf_yy)
	delete(sf_z)
	delete(sf_xz)
	delete(sf_yz)

end do

; Clean up
delete(U)
delete(V)
delete(modU)

; for output
Wx!0 = "time"
Wx&time = lags
Wx!1 = "level"
Wx&level = level
Wx!2 = "lat"
Wx&lat = lat
Wx!3 = "lon"
Wx&lon = lon

Wy!0 = "time"
Wy&time = lags
Wy!1 = "level"
Wy&level = level
Wy!2 = "lat"
Wy&lat = lat
Wy!3 = "lon"
Wy&lon = lon

Wz!0 = "time"
Wz&time = lags
Wz!1 = "level"
Wz&level = level
Wz!2 = "lat"
Wz&lat = lat
Wz!3 = "lon"
Wz&lon = lon

; Save output
 system("/bin/rm -f "+ncName)
 ncdf     = addfile(ncName,"c")  ; open output netCDF file
 ; make time and UNLIMITED dimension      ; recommended  for most applications
 filedimdef(ncdf,"time",-1,True)
 ncdf->time = lags
 ncdf->lat = lat
 ncdf->lon = lon
 ncdf->Wx = Wx
 ncdf->Wy = Wy
 ncdf->Wz = Wz
end