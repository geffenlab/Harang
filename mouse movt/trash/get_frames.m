function [ frames ] = get_frames( file_path, start_time, num_frames )
%get_frames gets frames

vid_reader = VideoReader(file_path);%,'CurrentTime',start_time);
f = readFrame(vid_reader);
vid_reader = VideoReader(file_path);
frames = uint8(zeros(size(f,1),size(f,2),num_frames));
frames(:,:,1) = f;

for i = 2 : num_frames
    if ~hasFrame(vid_reader)
        disp(['No more frames after reading ' num2str(i) ' frames.'])
        continue
    end
    disp(['Getting frame ' num2str(i) '/' num2str(num_frames)]);
    frames(:,:,i) = readFrame(vid_reader);
end

end

