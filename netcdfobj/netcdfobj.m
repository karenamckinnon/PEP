classdef netcdfobj < handle
%% Easy interface for reading netcdf files
%
%
% Example:
% o=netcdfobj('e:\meanstationpressure.cdf');
% o.vars
% pressure=o.vars.pressure.value;
%
%
% Aslak Grinsted 2009
    properties(SetAccess=protected) %setting not supported yet.
        fname= '';
        dims = [];
        vars = [];
        ncid = [];
        mode = '';
        convertfun=@(x)double(x);
    end
    properties(Hidden,SetAccess = private,GetAccess = private)
        ndims =[];
        nvars =[];
        natts =[];
        unlimdimID= [];
    end
    
    methods
        function obj = netcdfobj(fname,mode)
            if nargin<2
                mode='NC_NOWRITE';
            end
            obj.mode = upper(mode);
            obj.fname = fname;
            obj.ncid = netcdf.open(obj.fname,obj.mode);
            [obj.ndims,obj.nvars,obj.natts,obj.unlimdimID]= netcdf.inq(obj.ncid);
        end
        
        function value=get.ncid(obj)
            value=obj.ncid;
        end
        
        function value=get.dims(obj)
            if isempty(obj.dims)
                obj.dims=netcdfnamedlist(obj.ndims);
                for ii=1:obj.ndims
                    obj.dims.listdata{ii}=netcdfdim(obj,ii-1);
                end
            end
            value=obj.dims;
        end
        
        function value=get.vars(obj)
            if isempty(obj.vars)
                obj.vars=netcdfnamedlist(obj.nvars);
                for ii=1:obj.nvars
                    obj.vars.listdata{ii}=netcdfvar(obj,ii-1);
                end
            end
            value=obj.vars;
        end

%         function value=subsref(obj,s)
%             if s(1).type='.'
%             value=subsref(obj.vars,s);
%         end
              
        function close(obj)
            if ~isempty(obj.ncid)
                try
                    netcdf.close(obj.ncid);
                catch
                end
                obj.ncid=[];
            end
        end
        
        function delete(obj)
            obj.close();
        end
    end
    

end
