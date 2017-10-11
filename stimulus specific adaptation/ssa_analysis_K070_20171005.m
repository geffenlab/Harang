
%% 
load('F:\HarangData\K070_20171005_SSA_01.mat');
%%
stim_dur = ceil((stimInfo.tDur+stimInfo.ICI)/1e3 * exptInfo.fr);
raster = rasterize(spikes.raster, stim_dur, events.eventsOn);
%%
for i = 1 : 20
    subplot(2,1,1)
    imagesc(squeeze(mean(raster(1:10,:,8*(i-1)+1:8*i),2)))
    title(['AAAAAAAA Trial ' num2str(i)])
    subplot(2,1,2)
    imagesc(squeeze(mean(raster(1:10,:,8*(i-1+20)+1:8*(i+20)),2)))
    title(['BBBBBBBB Trial ' num2str(i+20)])
    pause
end
