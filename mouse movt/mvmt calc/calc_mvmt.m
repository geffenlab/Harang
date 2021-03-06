function fm = calc_mvmt(file)
% get_mvmt finds magnitude of movement in video file

v = VideoReader(file);
nf = v.NumberOfFrames;
v = VideoReader(file);
of = opticalFlowFarneback('NumPyramidLevels', 3,...
    'PyramidScale', 0.5, 'NumIterations', 3,...
    'NeighborhoodSize', 7, 'FilterSize', 15);
fm = zeros(1,nf);
for i = 1 : nf
    disp(['frame: ' num2str(i) '/' num2str(nf)])
    fr = readFrame(v);
    f = estimateFlow(of, fr(:,:,1));
    fm(i) = sum(sum(f.Magnitude));
end