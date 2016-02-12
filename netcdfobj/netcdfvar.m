classdef netcdfvar < handle
    % Helper netcdfobj.m
    % 
    %
    % Aslak Grinsted 2009
    properties(SetAccess=private) %setting not supported yet.
        name= '';
        id = [];
        xtype = [];
        dimids = [];
        atts = [];
        parent = [];
        value;
        rawvalue;
    end
    

    methods 
        function obj = netcdfvar(parent,id) 
            obj.parent=parent;
            obj.id=id;
            [name, xtype, dimids, numatts] = netcdf.inqVar(obj.parent.ncid,id); %get all meta data
            obj.xtype=xtype;
            obj.dimids=dimids;
            %obj.name.struct('id',ii-1,'xtype',xtype,'dimids',dimids,'getValue',@()netcdf.getVar(ncid,ii-1));
            obj.name=name;
            obj.atts=netcdfnamedlist(numatts);
            for jj=1:numatts
               obj.atts.listdata{jj}=netcdfatt(obj,jj-1);
            end
        end
        
        function value=get.rawvalue(obj)
            value=netcdf.getVar(obj.parent.ncid,obj.id); 
            %TODO: implement mechanism to only load a subset of data.
        end
        
        function value=get.value(obj)
            value=netcdf.getVar(obj.parent.ncid,obj.id,'double'); %TODO:allow a setting on netcdfobj which holds datatype. single as default?
            %INSERT nans: --------------- (todo: move to filler to private function?
            if ~isempty(obj.atts.missing_value); %TODO: allow a little slack due to roundoff errors
                misval = obj.atts.missing_value.value;
                if ischar(misval), misval=str2double(misval);end
                value(value == misval) = NaN;
            end
            if ~isempty(obj.atts.('_FillValue'))
                misval =obj.atts.('_FillValue').value;
                if ischar(misval), misval=str2double(misval);end
                value(value == misval) = NaN;
            end
            if ~isempty(obj.atts.('valid_range'))
                valid_range = obj.atts.('valid_range').value;
                value(value<min(valid_range))=nan;
                value(value>max(valid_range))=nan;
            end
            if ~isempty(obj.atts.('valid_min'))
                valid_min = obj.atts.('valid_min').value;
                value(value<valid_min)=nan;
            end
            if ~isempty(obj.atts.('valid_max'))
                valid_max = obj.atts.('valid_max').value;
                value(value>valid_max)=nan;
            end
            %---------------scale and offset
            if ~isempty(obj.atts.('scale_factor'))
                scale = obj.atts.('scale_factor').value;
                value=value*scale;
            end
            if ~isempty(obj.atts.('add_offset'))
                offset = obj.atts.('add_offset').value;
                value=value+offset;
            end
            %--------------
            if ~isempty(obj.atts.('calendar'))&&~isempty(obj.atts.('units'))
                value=cdfdate2num(obj.atts.units.value,obj.atts.calendar.value,value);
            end

            
            %TODO: make mechanism to select preferred data type 
            
            %TODO: implement mechanism to only load a subset of data.
        end

        
        
%         function set.value(obj,val) %NOT IMPLEMENTED YET
%             obj.val=val;
%             netcdf.putVar(val);
%         end
        
        
        function prettydisp(obj)
             s='';
             for ii=1:length(obj.dimids)
                 s=[s ' ' obj.parent.dims{obj.dimids(ii)+1}.name];
                 s=[s '(' num2str(obj.parent.dims{obj.dimids(ii)+1}.len) ')'];
             end
            disp(sprintf('%-15s \t[%s]',obj.name,deblank(s)))
        end
        
    end
    

end
