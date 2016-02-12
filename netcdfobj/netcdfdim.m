classdef netcdfdim < handle
    % Helper netcdfobj.m
    % 
    %
    % Aslak Grinsted 2009
    properties(SetAccess=private) %setting not supported yet.
        name= '';
        id = [];
        len = [];
        parent = [];
    end
    

    methods
        function obj = netcdfdim(parent,id) 
            obj.parent=parent;
            obj.id=id;
            [name, len] = netcdf.inqDim(obj.parent.ncid,id); %get all meta data
            obj.len=len;
            obj.name=name;
        end
        
        function prettydisp(obj)
            disp([obj.name ': ' num2str(obj.len)])
        end
        
    end
    

end
