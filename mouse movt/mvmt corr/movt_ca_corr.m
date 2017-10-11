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
%% see location of neurons correlated with motion
rho_corr = find(rho>0);
rho_corr_anti = find(rho<0);
loc_corr = spatialInfo.centroid(rho_corr,:);
loc_corr_anti = spatialInfo.centroid(rho_corr_anti,:);
scatter(loc_corr(:,1),loc_corr(:,2),50,rho(rho_corr),'+','LineWidth',2)
hold on
scatter(loc_corr_anti(:,1),loc_corr_anti(:,2),30,-1*rho(rho_corr_anti),...
    'o','LineWidth',2)
hold off; plotprefs; h = colorbar; axis square
ylabel(h, 'pos/neg correlation'); colormap jet;
title('location of pos/neg corr neurons on frame')
legend({'Ca-response pos corr w/ motion',...
    'Ca-response neg corr w/ motion'})
%% normalize motion
mot_interp_norm = mot_interp/std(mot_interp);
%% match bad frames
eps = 0.001;
thresh = 100; 
num_bad_f = sum(exptInfo.badFrames);
min_error = 1;
error = min_error + 1;
while abs(error) >= min_error
    num_thresh_mot = sum(mot_interp_norm > thresh);
    error = num_thresh_mot - num_bad_f;
    thresh = thresh + eps * error;
end
%% plot thresholded bad frames and bad frames
subplot(2,1,1); plot(exptInfo.badFrames); axis([0 size(raster,2) 0 inf])
index = find(mot_interp_norm>thresh);
subplot(2,1,2); plot(mot_interp_norm);
hold on; plot(index,mot_interp_norm(index),'*')
axis([0 length(mot_interp) 0 inf])
hold off
%% imagesc same thing
% imagesc([exptInfo.badFrames(1,1:end-ceil(extra_frames_sp));...
%     mot_interp_norm > thresh])
clf
imagesc([exptInfo.badFrames(1:end-ceil(extra_frames_sp));...
    mot_interp_norm > thresh])

