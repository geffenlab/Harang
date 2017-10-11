
aiv_dir = 'F:/videoData/';
mp3_dir = 'F:/mp4/';
mot_dir = 'F:/motionData/';
addpath(aiv_dir)
addpath(mp3_dir)
addpath(mp3_dir)

convert_vids_from_dir(aiv_dir, mp3_dir);

calc_mvmt_from_dir(mp3_dir, mot_dir);
