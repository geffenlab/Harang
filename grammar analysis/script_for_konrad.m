
load('K070_20170903_artificialGrammar_02.mat')
stim_dur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);
prev_frames = 4;
% sort spikes by events (ie stimulus presentations)
% (#neurons X stim_dur+prev_frames X #events)
raster = rasterize(spikes.raster, stim_dur + prev_frames, ...
    events.eventsOn - prev_frames);
% get average spike value for each event
raster_stim = squeeze(mean(raster(:,prev_frames+1:end,:), 2));
% average spike for first half minus second half of stimulus & ISI
raster_stim_half = squeeze(mean(raster(:,5:10,:), 2)) - ...
    squeeze(mean(raster(:,11:end,:), 2));
% show mean responses to stimulus 1
imagesc(raster_stim(:,stimInfo.order==1));
xlabel('stimulus 1 "events"'); ylabel('neurons');
%% view pca with time and stimulus detrend
rs=raster_stim_half-mean(raster_stim_half,2)*ones(1,2000);
smooth_stim=rs-conv2(rs,ones(1,100)/100,'same');
[~, score, ~] = pca(smooth_stim');
% pca by time
figure
subplot(1,2,1); title('colored by time')
scatter3(score(:,1),score(:,2),score(:,3),32,1:size(score,1),'.')
colormap jet; colorbar
box on; axis([-1 1 -1 1 -1 1]); axis vis3d
% pca by stimulus
subplot(1,2,2); title('colored by stimulus');
for i = 1 : 3
    s = find(stimInfo.order==i);
    scatter3(score(s,1),score(s,2),score(s,3),32,'.');
    hold on
end
box on; axis([-1 1 -1 1 -1 1]); axis vis3d
legend({'stim 1','stim 2','stim 3'})
