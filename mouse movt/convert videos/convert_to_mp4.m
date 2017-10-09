function [  ] = convert_to_mp4( reader, filepath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

v = VideoWriter(filepath, 'MPEG-4');

open(v)
while hasFrame(reader)
    writeVideo(v, readFrame(reader));
end
close(v)

end

