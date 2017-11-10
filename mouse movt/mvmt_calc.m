
aiv_dir = 'F:/videoData/';
mp4_dir = 'F:/mp4/';
mot_dir = 'F:/motionData/';
addpath(aiv_dir)
addpath(mp4_dir)
addpath(mp4_dir)

convert_vids_from_dir(aiv_dir, mp4_dir);

calc_mvmt_from_dir(mp4_dir, mot_dir);
