function time=cdfdate2num(units,calendar,time)
% Convert netcdf times to datenums
%
% USAGE: time = cdfdate2num(units,calendar,time)
%
% Warning: strange calendars like 360_day are simply stretched by 365.24237/360
% This means datevec will no longer return the right hour and minutes etc. 
%
% Example:
%  time= cdfdate2num('hours since 1856-01-03 -07:00','360_day',5)
%  datestr(time)
%
%
%
% Aslak Grinsted 2013


%inspired by reference: http://netcdf4-python.googlecode.com/svn/trunk/docs/netcdftime.netcdftime-pysrc.html#utime.num2date
% units=n.vars.time.atts.units.value;
% calendar=n.vars.time.atts.calendar.value;

pat='(?<resolution>\w+)\s+since\s+(?<offset>.*)';

units=regexpi(units,pat,'names','once');

if isempty(units)
    error('Unknown time units.')
end

switch lower(units.resolution)
    case 'days'
    case 'hours'
        time=time/24;
    case 'minutes'
        time=time/24/60;
    case 'seconds'
        time=time/24/60/60;
    otherwise
        error('unsupported units');
end

offset=datenum(units.offset);

yearlen=365.24237;
switch lower(calendar)
    case {'julian','standard','gregorian','proleptic_gregorian'}
        time=offset+time;
    case {'noleap','365_day'}
        time=offset+time*yearlen/365;
    case {'all_leap','366_day'}
        time=offset+time*yearlen/366;
    case '360_day'
        time=offset+time*yearlen/360;
    otherwise
        error('unsupported calendar');
end
