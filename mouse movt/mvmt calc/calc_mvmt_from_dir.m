function calc_mvmt_from_dir( vid_dir, mot_dir )
%calc_mvmt_from_dir 

dirs = dir(vid_dir);
for i = 1 : length(dirs)
    tic
    d = dirs(i).name;
    if d(1) == '.'
        continue
    end
    mot_file = [mot_dir d(1:end-4) '.mat'];
    disp(mot_file)
    if exist(mot_file,'file')
        disp(['File ' mot_file ' already exists. Skipping motion analysis.'])
        continue
    end
    disp(['Analyzing file ' d '...'])
    fm = calc_mvmt([vid_dir d]);
    save(mot_file, 'fm')
    toc
end

end

