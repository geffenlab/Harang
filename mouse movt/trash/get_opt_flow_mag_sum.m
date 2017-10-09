function flow_mag = get_opt_flow_mag_sum(frames, frames_per_flow)
%% Optical flow of mouse head videos
% Harang Ju
% 09/05/2017
% Goal: find times of large head movement of a mouse

%% Optical Flow
opt_flow = opticalFlowFarneback('NumPyramidLevels', 3,...
    'PyramidScale', 0.5, 'NumIterations', 3,...
    'NeighborhoodSize', 7, 'FilterSize', 15);
% Farneback for dense feature set
% Lukas-Kanade for sparse feature set
% assumes displacement of image contents between two nearby frames is small
% Horn-Schunck assumes smoothness in flow over whole image

num_frames = size(frames,3);
flow_mag = zeros(1,num_frames);
for i = 1 : frames_per_flow : num_frames
    disp(['frame: ' num2str(i) '/' num2str(num_frames)])
    flow = estimateFlow(opt_flow, frames(:,:,i));
    flow_mag(i) = sum(sum(flow.Magnitude));
end

