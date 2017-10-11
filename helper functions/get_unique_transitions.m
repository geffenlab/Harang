function [ uT ] = get_unique_transitions( trans )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

uT = unique(trans,'rows');
uT = sortrows(uT,[2,1]);

end
