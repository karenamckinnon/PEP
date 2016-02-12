classdef netcdfnamedlist < handle
    % Helper netcdfobj.m
    %
    %
    % Aslak Grinsted 2009
    properties(Hidden)
        listdata = [];
    end
    
    methods
        function obj = netcdfnamedlist(noitems)
            obj.listdata=cell(noitems,1);
        end
        
        function value=subsref(obj,s)
            value=[];
            switch s(1).type
                case {'()','{}'}
                    value = obj.listdata{cell2mat(s(1).subs)};
                case '.'
                    fieldidx=obj.fieldindex(s(1).subs);
                    if ~isempty(fieldidx)
                        value=obj.listdata{fieldidx};
                    end
%                     for ii=1:length(obj.listdata)
%                         if strcmp(s(1).subs,obj.listdata{ii}.name)
%                             value=obj.listdata{ii};
%                             break
%                         end
%                     end
            end
            if isempty(value)
                return
                %error(sprintf('''%s'' not found.',strcat(s.subs)));
            end
            for jj=2:length(s)
                value=subsref(value,s(jj));
            end
        end
        
        function value=fieldindex(obj,s)
            value=[];
            for ii=1:length(obj.listdata)
                if strcmp(s,obj.listdata{ii}.name)
                    value=ii;
                    return
                end
            end
        end
        
        %
        %         function value=subsindex(obj,name)
        %             for ii=1:length(listdata)
        %                 if strcmp(name,listdata{ii}.name)
        %                     value=ii;
        %                     return
        %                 end
        %             end
        %             error(sprintf('''%s'' not found.',name));
        %         end
        
        
        function display(obj)
            for ii=1:length(obj.listdata)
                prettydisp(obj.listdata{ii});
            end
        end
        
    end
    
    
end
