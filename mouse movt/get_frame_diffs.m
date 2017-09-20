function [ flow ] = get_frame_diffs( frames, frames_per_flow )
%get_frame_diffs calculates the sum of
%   a subtraction of the previous subtraction from all frames

num_frames = size(frames,3);
flow = zeros(1,num_frames);
prev_frame = frames(:,:,1);
for i = 1 : frames_per_flow : num_frames
    disp(['frame: ' num2str(i) '/' num2str(num_frames)])
    flow(i) = sum(sum(frames(:,:,i) - prev_frame));
    prev_frame = frames(:,:,i);
end

end

