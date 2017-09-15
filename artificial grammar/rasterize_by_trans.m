function [ raster, uT, nT ] = rasterize_by_trans( data, trans )
%Rasterize creates a raster from data
%   raster is n x 
%   uT - unique transitions, sorted
%   nT - num transitions

% unique transitions
uT = unique(trans,'rows');
uT = sortrows(uT,[2,1]);
raster = zeros(size(data,1), size(data,2), length(uT));
nT = zeros(1,length(uT));
for t = 1:length(uT)
    rows = find(trans(:,1)==uT(t,1) & trans(:,2)==uT(t,2));
    nT(t) = length(rows);
    raster(:,:,t) = mean(data(:,:,rows),3);
end

end