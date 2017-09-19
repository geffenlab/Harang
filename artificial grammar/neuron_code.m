load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

%%
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);
raster = rasterize(spikes.raster, stimDur, events.eventsOn);
raster_stim = squeeze(mean(raster, 2));
input_raster = [stimInfo.order==1;stimInfo.order==2;stimInfo.order==3];
%% standardize raster
raster_stim_z = raster_stim - repmat(mean(raster_stim,2),[1 2000]);
raster_stim_z = raster_stim_z ./ repmat(std(raster_stim,0,2),[1 2000]);
%% get most responsive neurons
[~,index] = sort(mean(raster_stim,2),'descend');
top = 316;
raster_stim_top = raster_stim(index(1:top),:);
figure(1);
bar(mean(raster_stim(index,:),2))
%%
c = corr(input_raster', raster_stim_z');
threshold = 0.05;
figure(2)
subplot(2,2,1);
imagesc(c); title('correlations b/t stimuli & neurons');
plotprefs; colorbar; colormap('gray'); set(gca,'ytick',1:3);
xlabel('neurons (avg spikes)'); ylabel('stimuli');
subplot(2,2,2);
imagesc(c>threshold); title(['correlations > ' num2str(threshold)]);
plotprefs; colorbar; colormap('gray'); set(gca,'ytick',1:3);
subplot(2,2,4);
imagesc((c<(-1*threshold))); title(['correlations < -' num2str(threshold)]);
plotprefs; colorbar; colormap('gray'); set(gca,'ytick',1:3);
subplot(2,2,3)
c_t = (c>threshold) - (c<(-1*threshold));
imagesc(c_t); plotprefs; colorbar; colormap('gray'); set(gca,'ytick',1:3);
title(['correlations > ' num2str(threshold) ' & '...
    'correlations < -' num2str(threshold)]);
%%
% scatter(stimInfo.order,raster_stim(133,:),'*');
figure(3);clf
stim = 1;
stim_resp = find(c_t(stim,:));
n_resp = length(stim_resp);
n_row = 3;
for r = 1 : n_resp
    subplot(n_row,ceil(n_resp/n_row),r)
    boxplot(raster_stim_top(stim_resp(r),:)',stimInfo.order');
    plotprefs; set(findobj(gca,'type','line'),'linew',2);
end
text(0.5,0.5,'spikes','Rotation',90,'FontSize',16)
text(0.5,0.5,'stimuli','FontSize',16)
text(0.5,0.5,'spikes from neurons most correlated to stimuli 1','FontSize',16)
%% population coding
figure(4);clf
top_corr = 10;
[~,top_corr_index] = sort(max(abs(c),[],1), 'descend');
[idx, C, sumd] = kmeans(raster_stim(top_corr_index(1:top_corr),:)',3);
plot(stimInfo.order(1:100),'o');
hold on; plot(idx(1:100),'.');
hold off
plotprefs; axis([0 100 0 4])

