function [T,d,flag] = get_data(t,varname,stationRecord)
% [T,dates,flag] = get_data(yr_begin,yr_end,varname,stationRecord)
% Extracts data for a given variable from the stationRecord structure,
% subsetting for the desired time period.
%
% In:
%    yr_begin -
%    yr_end   -
%    varname  -
%    stationRecord -
% Out:
%    T - the data for varname
%    dates - the dates for those data, in matlab datenum format
%    flag  - returns zero if the record doesn't contain data
%            matching the input criteria.

flag = 1; % set to zero if record doesn't fit date criteria
d = stationRecord.dates;
T = 0;
% consider stations that have data over entire set of years we consider
if (min(d) <+ t(1) && max(d) >= t(end))

    try % Some weird files have numel(yr_range) = 0
        t1 = find(d == t(1));
        t2 = find(d == t(end));

        fn = char(['stationRecord.' varname '.data']);
        T = eval(fn);
        T = double(T(t1:t2));
        d = d(t1:t2);

    catch err
        flag = 0;
    end

else
    flag = 0;
end

return