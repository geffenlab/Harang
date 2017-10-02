
load('K070_20170903_artificialGrammar_02.mat')
stim_dur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);
prev_frames = 4;
% sort spikes by events (ie stimulus presentations)
% (#neurons X stim_dur+prev_frames X #events)
raster = rasterize(spikes.raster, stim_dur + prev_frames, ...
    events.eventsOn - prev_frames);
% get average spike value for each event
raster_stim = squeeze(mean(raster(:,prev_frames+1:end,:), 2));
% show mean responses to stimulus 1
imagesc(raster_stim(:,stimInfo.order==1));
xlabel('stimulus 1 "events"'); ylabel('neurons');
% sort data by transitions & average across instances of unique trans
% (#neurons X stim_dur+prev_frames X #unique transitions)
trans = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
uniq_trans = unique(trans,'rows');
uniq_trans = sortrows(uniq_trans,[2,1]);
[raster_trans, num_trans] = rasterize_by_trans(raster, trans, uniq_trans);
