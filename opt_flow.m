
directory = 'F:/videoData/';
filename = '20170905K074_FRA_log_986_02';
filetype = '.avi';
filepath = [directory filename filetype];
starttime = 20*60 + 12; % seconds into video
numFrames = 2 * 30;
framesPerFlow = 2;

addpath(directory);
flow = get_opt_flow_mag_sum(filepath, starttime, numFrames,...
    framesPerFlow);

