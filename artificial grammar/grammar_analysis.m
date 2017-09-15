load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

uniqStim = unique(stimInfo.order); % how many words
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);

%%
raster = rasterize(spikes.raster, stimDur, events.eventsOn);
% epoch
epoch_size = 100;
epoch_onsets = epoch_size*(0:19)+1;
raster_epoch = rasterize_by_epoch(raster, epoch_size, epoch_onsets);
avg_resp_epoch = mean(mean(raster_epoch,4),3);
% transitions
trans = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
[raster_trans, uT, nT] = rasterize_by_trans(raster, trans);

%%
raster_trans_epoch = zeros(size(raster,1),size(raster,2),...
    );
for e = 1 : length(epoch_onsets)
    
end

%% Plot each transition raster order by most responsive neurons
figure
grammar_flat = reshape(stimInfo.grammar', [1 9]);
for t = 1:length(uT)
    subplot(3,3,t)
    [~,index] = sort(max(raster_trans(:,:,t),[],2),'descend');
    imagesc(raster_trans(index,:,t), [0 0.5])
    title([ num2str(uT(t, 1)) '->' num2str(uT(t,2))])
    xlabel('frames'); ylabel('neurons')
    title(['P(' num2str(uT(t,1)) '->' ...
        num2str(uT(t,2)) ')=' num2str(grammar_flat(t))]);
    colorbar
end

%% Work out mean response to each transition across neurons
mtr = squeeze(mean(mean(raster_trans,1),2));
mtr = reshape(mtr,[3,3])';
c = corr(stimInfo.grammar(:),mtr(:));
disp(['corr b/t grammar & mean response = ' num2str(c)])

