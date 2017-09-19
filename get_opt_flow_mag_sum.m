function flowMags = get_opt_flow_mag_sum(filepath, starttime, ...
    numFrames, framesPerFlow)
%% Optical flow of mouse head videos
% Harang Ju
% 09/05/2017
% Goal: find times of large head movement of a mouse

%% Initialization
vidReader = VideoReader(filepath, 'CurrentTime', starttime);

%% Optical Flow
optFlow = opticalFlowFarneback('NumPyramidLevels', 3,...
    'PyramidScale', 0.5, 'NumIterations', 3,...
    'NeighborhoodSize', 7, 'FilterSize', 15);
% Farneback for dense feature set
scale = 5;
% optFlow = opticalFlowLKDoG('NumFrames',3);
% Lukas-Kanade for sparse feature set
% assumes displacement of image contents between two nearby frames is small
% scale = 25;
% optFlow = opticalFlowHS;
% Horn-Schunck assumes smoothness in flow over whole image

flowMags = zeros(1,numFrames);

for i = 1 : numFrames
    for j = 1 : framesPerFlow - 1
        % skip frames
        readFrame(vidReader);
    end
    frame = readFrame(vidReader);
    flow = estimateFlow(optFlow, frame);
    imshow(frame);
    hold on
    plot(flow,'DecimationFactor',[80 80],'ScaleFactor',scale)
    hold off
    flowMags(i) = sum(sum(flow.Magnitude));
end

