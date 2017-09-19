load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

uniqStim = unique(stimInfo.order); % how many words
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);

%%
raster = rasterize(spikes.raster, stimDur, events.eventsOn);
% epoch
epoch_size_f = 100;
epoch_onsets_f = epoch_size_f*(0:19)+1;
raster_epoch = rasterize_by_epoch(raster, epoch_size_f, epoch_onsets_f);
avg_resp_epoch = mean(mean(raster_epoch,4),3);
% transitions
trans = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
uT = get_unique_transitions(trans);
[raster_trans, nT] = rasterize_by_trans(raster, trans, uT);

%%
epoch_size_s = 200;
epoch_onsets_s = epoch_size_s*(0:9)+1; 
epoch_onsets_s(end) = epoch_onsets_s(end)-1;
raster_trans_epoch = rasterize_by_trans_epoch(raster,trans,uT,...
    epoch_size_s,epoch_onsets_s);
%%
trans_inst = zeros(length(uT),length(epoch_onsets_s));
for t = 1 : length(uT)
    for e = 1:length(epoch_onsets_s)
        stims = epoch_onsets_s(e):epoch_onsets_s(e)+epoch_size_s-1;
        trans_inst(t,e) = sum(trans(stims,1)==uT(t,1) & ...
            trans(stims,2)==uT(t,2));
    end
end
%%
subplot(1,2,1)
raster_trans_epoch_mean = squeeze(mean(mean(raster_trans_epoch,1),2));
imagesc([raster_trans_epoch_mean; ...
    mean(raster_trans_epoch_mean, 1)])
caxis([0 0.03])
colormap('gray'); xlabel('epoch'); colorbar; plotprefs;
title('mean spikes for 316 neurons')
set(gca,'yticklabel',{...
    'P(1|1)=0.1','P(2|1)=0.7','P(3|1)=0.3',...
    'P(1|2)=0.2','P(2|2)=0.2','P(3|2)=0.6',...
    'P(1|3)=0.7','P(2|3)=0.1','P(3|3)=0.1', 'average'})
subplot(1,2,2)
imagesc([trans_inst; mean(trans_inst, 1)]); colorbar; plotprefs;
colormap('gray'); xlabel('(each 200 stimuli presentations)');
title('number of stimuli transitions'); set(gca,'yticklabel',{})
%%
[~,index] = sort(max(raster_trans(:,:,2),[],2),'descend');
raster(index(51:end),:,:) = [];
%%


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

