function [ raster ] = rasterize_by_epoch( data, dur, onsets )
%Rasterize creates a raster from data
%   raster is n x length(onsets) x dur x ...

raster = zeros(size(data,1), size(data,2), dur, length(onsets));
for i = 1 : length(onsets)
    raster(:,:,:,i) = data(:, :, onsets(i):onsets(i)+dur-1);
end

end