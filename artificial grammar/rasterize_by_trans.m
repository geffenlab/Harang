function [ raster, num_trans ] = rasterize_by_trans(data, trans, uniq_trans)
%Rasterize_by_trans sorts responses by unique transitions and averages the
%responses across frames
%   raster is #neurons x stimulus duration (ie size(data,2)) x #unique
%   transitions
%   num_trans - num transitions

% unique transitions
raster = zeros(size(data,1), size(data,2), length(uniq_trans));
num_trans = zeros(1,length(uniq_trans));
for t = 1:length(uniq_trans)
    rows = find(trans(:,1)==uniq_trans(t,1) & trans(:,2)==uniq_trans(t,2));
    num_trans(t) = length(rows);
    raster(:,:,t) = mean(data(:,:,rows),3);
end

end