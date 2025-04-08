classdef Tophat

    properties
        SE;
        result;
    end
    
    methods
        
        function obj = Tophat
            obj.SE = strel( 'disk', 3 );
        end
        
        % inImg must range from 0 to 1
        function obj = process(obj, inImg)
            obj.result = imtophat( inImg, obj.SE );
        end
        
    end


end