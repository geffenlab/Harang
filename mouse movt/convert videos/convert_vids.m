function convert_vids

directory = 'F:/videoData/';
new_dir = 'F:/mp4/';
addpath(directory);

dirs = dir(directory);
for i = 1 : length(dirs)
    tic
    d = dirs(i).name;
    if d(1) == '.'
        continue
    end
    disp(['Converting file ' d '.'])
    v = VideoReader([directory d]);
    disp([new_dir d(1:end-4)])
    convert_to_mp4(v, [new_dir d(1:end-4)]);
    toc
end
