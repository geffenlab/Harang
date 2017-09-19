load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

%%
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI) * exptInfo.fr);
raster = rasterize(spikes.raster, stimDur, events.eventsOn);
raster_stim = squeeze(mean(raster, 2));
input_raster = [stimInfo.order==1;stimInfo.order==2;stimInfo.order==3];
%%
[~,index] = sort(mean(raster_stim,2),'descend');
top = 100;
raster_stim_top = raster_stim(index(1:top),:);
figure(1);
bar(mean(raster_stim(index,:),2))
%%
c = corr(input_raster', raster_stim_top');
threshold = 0.07;
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

