
%% 
load('F:\HarangData\K070_20171005_SSA_01.mat');
%%
stim_dur = ceil((stimInfo.tDur+stimInfo.ICI)/1e3 * exptInfo.fr);
raster = rasterize(spikes.raster, stim_dur, events.eventsOn);
%%
