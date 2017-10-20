classdef dh
    %DataHelper useful functions 
    
    properties
    end
    
    methods (Static)
        
        function raster = rasterize(data, dur, onsets)
            %Rasterize creates a raster from data
            %   raster is n x length(onsets) x dur
            raster = zeros(size(data,1), dur, length(onsets));
            for i = 1 : length(onsets)
                raster(:,:,i) = data(:, onsets(i):onsets(i)+dur-1);
            end
        end
        
        function new_dat = split_concat(dat)
             % opposite of rasterize
             n = size(dat,3);
             new_dat = zeros(size(dat,1),n*size(dat,2));
             for i = 1 : n
                 on = (i-1)*size(dat,2)+1;
                 off = on+size(dat,2)-1;
                 new_dat(:,on:off) = dat(:,:,i);
             end
        end
        
        function resps = max_in_windows( responses, duration, onsets )
            %max_in_windows finds the max responses in each window
            %   window = onset : onset + duration - 1
            resps = zeros(size(responses,1),length(onsets));
            for i = 1 : size(responses,1)
                for j = 1 : length(onsets)
                    resps(i,j) = max(...
                        responses(i,onsets(j) : onsets(j)+duration-1));
                end
            end
        end
    end
    
end

