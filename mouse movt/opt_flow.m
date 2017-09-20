
directory = 'F:/videoData/';
addpath(directory);
file_name = '20170905K074_FRA_log_986_02';
file_type = '.avi';
file_path = [directory file_name file_type];
% start_time = 20*60 + 12;
start_time = 20*60 + 10;
% start_time = 1;
% start_time = 1;
% num_frames = 2 * 30;
num_frames = 30 * 20;
frames_per_flow = 1;

%%
frames = get_frames(file_path, start_time, num_frames);

%%
% opt_flow = get_opt_flow_mag_sum(frames, frames_per_flow);
diff = get_frame_diffs(frames, frames_per_flow);
motion = diff;

%%
plot(motion./mean(motion(1:100)),'.-')

%% find frames to remove
threshold = 2.5;
frames_to_remove = find(motion./mean(motion(1:100)) > threshold);
frames_filtered = frames;
frames_filtered(:,:,frames_to_remove) = [];
plot(motion./mean(motion(1:100)),'.-'); hold on
plot(frames_to_remove, motion(frames_to_remove)./mean(motion(1:100)),'o')
hold off

%% show skipped frames
j = 1;
for i = 1 : size(frames,3)
    subplot(1,2,1); imshow(frames(:,:,i)); title('unfiltered');
    subplot(1,2,2);
    if i == frames_to_remove(j)
        j = j + 1;
        title(['SKIP ' num2str(j)]); pause(0.2);
    else
        imshow(frames(:,:,i));
    end
    pause(0.01)
end


