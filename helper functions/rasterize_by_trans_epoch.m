function [ raster ] = rasterize_by_trans_epoch( data, ...
    trans, uT, dur, onsets )
%Rasterize creates a raster from data
%   raster is n x length(onsets) x dur x ...

raster = zeros(size(data,1), size(data,2),...
    length(uT), length(onsets));
for e = 1 : length(onsets)
    stimuli = onsets(e) : onsets(e)+dur-1;
    raster(:,:,:,e) = ...
        rasterize_by_trans(data(:,:,stimuli), trans(stimuli,:), uT);
end

end