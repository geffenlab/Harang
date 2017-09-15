function [ raster ] = rasterize( data, dur, onsets )
%Rasterize creates a raster from data
%   raster is n x length(onsets) x dur

raster = zeros(size(data,1), dur, length(onsets));
for i = 1 : length(onsets)
    raster(:,:,i) = data(:, onsets(i):onsets(i)+dur-1);
end

end
