classdef netcdfatt < handle
    % Helper netcdfobj.m
    %
    %
    % Aslak Grinsted 2009
    properties(SetAccess=private) %setting not supported yet.
        name= '';
        id = [];
        value = [];
        parent = [];
    end
    
    methods
        function obj = netcdfatt(parent,id)
            obj.parent=parent;
            obj.id=id;
            ncid=parent.parent.ncid;
            obj.name=netcdf.inqAttName(ncid,parent.id,id);
            obj.value=netcdf.getAtt(ncid,parent.id,obj.name);
        end
        
        function prettydisp(obj)
            if isnumeric(obj.value)&&(numel(obj.value)<10)
                disp([obj.name ': ' mat2str(obj.value)]);
            elseif ischar(obj.value)&&(size(obj.value,1)==1)
                disp([obj.name ': ' obj.value(1,1:min(end,60))]);
            else
                disp([obj.name ': [' class(obj.value) ']'] );
            end
        end
        
    end
    
    
    
end
