load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

uniqStim = unique(stimInfo.order); % how many words
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);

%%
raster = rasterize(spikes.raster, stimDur+4, events.eventsOn-4);
% epoch
epoch_size_f = 100;
epoch_onsets_f = epoch_size_f*(0:19)+1;
raster_epoch = rasterize_by_epoch(raster, epoch_size_f, epoch_onsets_f);
avg_resp_epoch = mean(mean(raster_epoch,4),3);
% transitions
trans = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
uT = get_unique_transitions(trans);
[raster_trans, nT] = rasterize_by_trans(raster, trans, uT);
% stimulus
raster_stim = squeeze(mean(raster, 2));

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

%% plot individual neuron responses
clf
n = 5;
o = stimInfo.order;
up = 0.06;
for j = 1 : floor(316/n)
figure
up = max(mean(mean(raster(:,:,:),1),2));
for i = 1 : n
    neuron = n*(j-1)+i;
    r = squeeze(raster(neuron,:,:));
    % plot mean responses to stimuli
    subplot(n,4,4*(i-1)+1); hold on
    axis([0 inf 0 inf]); plotprefs;
    ylabel(['neuron ' num2str(neuron)]);
    plot(mean(r(:,o==1),2),'.-')
    plot(mean(r(:,o==2),2),'.-')
    plot(mean(r(:,o==3),2),'.-')
    if i == 1
        title('mean spikes'); legend({'word 1','word 2','word 3'})
    end
    % plot mean responses to stimuli 1
    subplot(n,4,4*(i-1)+2); hold on
    axis([0 inf 0 inf]); plotprefs
    plot(mean(r(:,(o(1:end-1)==1 & o(2:end)==1)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==2 & o(2:end)==1)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==3 & o(2:end)==1)==1),2),'.-')
    if i == 1
        title('spikes to stim 1 from...');
    end
    subplot(n,4,4*(i-1)+3); hold on
    axis([0 inf 0 inf]); plotprefs
    plot(mean(r(:,(o(1:end-1)==1 & o(2:end)==2)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==2 & o(2:end)==2)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==3 & o(2:end)==2)==1),2),'.-')
    if i == 1
        title('spikes to stim 2 from...');
    end
    subplot(n,4,4*(i-1)+4); hold on
    axis([0 inf 0 inf]); plotprefs
    plot(mean(r(:,(o(1:end-1)==1 & o(2:end)==3)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==2 & o(2:end)==3)==1),2),'.-')
    plot(mean(r(:,(o(1:end-1)==3 & o(2:end)==3)==1),2),'.-')
    if i == 1
        title('spikes to stim 3 from...');
    end
end
saveas(gcf,['spikes by transition from n=' num2str(neuron-n+1) ' to ' ...
        num2str(neuron) '.png'])
end

%% pca
p = pca(raster_stim');
p_orthnorm = inv(diag(std(raster_stim')))* p;
% rsp = p' * raster_stim;
rsp = raster_stim * p;
o = stimInfo.order;
%%
for i = 1 : 3
    s = find(stimInfo.order==i);
    scatter3(rsp(1,s),rsp(2,s),rsp(3,s),'.');
    hold on
end
hold off
legend({'ones','twos','threes'})
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
plotprefs
%%
tr = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2; 1 3; 2 3; 3 3];
for t = 1 : 9
    s = (o(1:end-1)==tr(t,1) & o(2:end)==tr(t,2)) == 1;
    scatter3(rsp(1,s),rsp(2,s),rsp(3,s),'.');
    hold on
end
hold off
legend({...
    '1->1 (P=0.1)','2->1 (P=0.7)','3->1 (P=0.3)',...
    '1->2 (P=0.2)','2->2 (P=0.2)','3->2 (P=0.6)',...
    '1->3 (P=0.7)','2->3 (P=0.1)','3->3 (P=0.1)'
    })
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
plotprefs

%% 
[~, score, ~] = pca(raster_stim');
scatter3(score(:,1),score(:,2),score(:,3),32,1:size(score,1))
colormap jet; colorbar
box on; axis([-1 1 -1 1 -1 1]); axis vis3d
plotprefs; hold on
for i = 1 : 3
    s = find(stimInfo.order==i);
    scatter3(score(s,1),score(s,2),score(s,3),'.');
end
hold off
legend({'time-sorted','ones','twos','threes'})
%% sort by time and cluster
range = 1:2000;
[coeff, score, latent] = pca(raster_stim(:,range)');
for i = 1 : 3
    s = find(stimInfo.order(range)==i);
    scatter3(score(s,1),score(s,2),score(s,3),'.');
    hold on
end
hold off; legend({'ones','twos','three'})
plotprefs; box on;
axis([-1 1 -1 1 -1 1]); axis vis3d
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')


