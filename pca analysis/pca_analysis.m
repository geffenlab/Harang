load('F:\HarangData\K074_20170905_2P_FRA_01.mat')
load('F:\HarangData\20170905K074_FRA_log_986_02_movt.mat')
%%
stimInfo = stimInfo.stimInfo;
raster = spikes.raster;
raster = calcium.npilSubTraces;
%% normalize spikes
raster_norm = log(raster);
%% pca
[coeff, score, ~] = pca(raster_norm');
imagesc(coeff(:,1:100)); title('principal components'); axis square
ylabel('coefficients'); xlabel('components');
colormap gray; colorbar
plotprefs
%% show projections
scatter3(score(:,1),score(:,2),score(:,3));%,32,1:size(score,1))
colormap jet; colorbar
box on; axis([-1 1 -1 1 -1 1]); axis vis3d
%% color by freq of stim
stim_atten = stimInfo.index(stimInfo.order,1);
stim_dur = ceil((stimInfo.tDur/1e3+.9) * exptInfo.fr);
events.eventsOn;
c = zeros(1,47000);
c_index = repelem(events.eventsOn',stim_dur) + ...
    repmat(0:stim_dur-1,[1 length(events.eventsOn)]);
c(c_index) = repelem(stim_atten, stim_dur);
s_index = c_index(1:10:length(c_index));
scatter3(score(s_index,1),score(s_index,2),score(s_index,3),40,...
    c(s_index),'o');
colormap jet; colorbar; title('colored by stimulus freqency')
box on; axis vis3d
%% color by attenuation of stim
stim_atten = stimInfo.index(stimInfo.order,2);
stim_dur = ceil((stimInfo.tDur/1e3+.9) * exptInfo.fr);
events.eventsOn;
c = zeros(1,47000);
c_index = repelem(events.eventsOn',stim_dur) + ...
    repmat(0:stim_dur-1,[1 length(events.eventsOn)]);
c(c_index) = repelem(stim_atten, stim_dur);
s_index = c_index(1:20:length(c_index));
scatter3(score(s_index,1),score(s_index,2),score(s_index,3),40,...
    c(s_index),'o');
colormap jet; colorbar; title('colored by stimulus attenuation')
box on; axis([-1 1 -1 1 -1 1]); axis vis3d
%% color by motion
% use movt_ca_corr code
c = log(mot_interp);
s_index = 1:100:size(score,1)-1000;
scatter3(score(s_index,1),score(s_index,2),score(s_index,3),40,...
    c(s_index),'o');
colormap jet; colorbar; title('colored by motion')
box on; axis vis3d
