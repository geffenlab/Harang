load('F:\HarangData\K074_20170905_2P_FRA_01.mat')
load('F:\HarangData\20170905K074_FRA_log_986_02_movt.mat')
stimInfo = stimInfo.stimInfo;
%%
raster = spikes.raster;
subplot(2,1,1); plot(mean(raster)); axis([0 size(raster,2) 0 inf])
subplot(2,1,2); plot(fm); axis([0 length(fm) 0 inf])
%%
fr_sp = exptInfo.fr; % spike data frame rate
dur_sp = size(raster,2) / fr_sp; % spikes duration in seconds
fr_mot = 20;
dur_mot = size(fm,2) / fr_mot;
extra_frames_sp = (dur_sp-dur_mot)*fr_sp;
%% interpolate motion data to fit spike data 
x = 1/fr_mot:1/fr_mot:dur_mot;
xq = 1/fr_mot:1/fr_sp:dur_mot;
mot_interp = [0 interp1(x,fm,xq)];
% %% normalize motion data
% mot_interp = zscore(mot_interp);
%%
subplot(2,1,1); plot(mean(raster)); axis([0 size(raster,2) 0 inf])
subplot(2,1,2); plot(mot_interp); axis([0 length(mot_interp) -inf inf])
%%
rho = corr(raster(:,1:end-ceil(extra_frames_sp))',mot_interp');
clf
bar(rho)
mean(rho)
%%
subplot(2,1,1); plot(raster(1,:)); axis([0 size(raster,2) 0 inf])
subplot(2,1,2); plot(mot_interp); axis([0 length(mot_interp) 0 inf])
%%
subplot(2,1,1); plot(exptInfo.badFrames(1,:)); axis([0 size(raster,2) 0 inf])
subplot(2,1,2); plot(mot_interp); axis([0 length(mot_interp) 0 inf])

