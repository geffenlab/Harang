function convert_vids_from_dir(from_dir, to_dir)
% convert_vids converts videos to mp4 in from_dir and saves to to_dir

addpath(from_dir);

dirs = dir(from_dir);
for i = 1 : length(dirs)
    tic
    d = dirs(i).name;
    if d(1) == '.'
        continue
    end
    dest_d = [to_dir d(1:end-4) '.mp4'];
    if exist(dest_d,'file')
        disp(['File ' dest_d ' already exists. Skipping video conversion.'])
        continue
    end
    disp(['Converting file ' d '.'])
    v = VideoReader([from_dir d]);
    disp([to_dir d(1:end-4)])
    convert_to_mp4(v, [to_dir d(1:end-4)]);
    toc
end
