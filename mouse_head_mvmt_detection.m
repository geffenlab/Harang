%% Optical flow of mouse head videos
% Harang Ju
% 09/05/2017
% Goal: find times of large head movement of a mouse

%% Initialization
directory = 'F:/videoData/';
filename = '20170905K074_spont10mins';
filetype = '.avi';

filepath = [directory filename filetype];
addpath(directory);
vidReader = VideoReader(filepath);

%% Optical Flow
% optFlow = opticalFlowFarneback;
% scale = 2;
optFlow = opticalFlowLKDoG('NumFrames',3);
scale = 25;
% optFlow = opticalFlowHS;

while hasFrame(vidReader)
    frame = readFrame(vidReader);
    flow = estimateFlow(optFlow, frame);
    imshow(frame);
    hold on
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',scale)
    hold off
end
