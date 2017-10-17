
%%
load('F:/HarangData/K070_20171013_SSA_03.mat')
%%
stim_dur = stimInfo.tDur*1e-3*exptInfo.fr;
ici_dur = stimInfo.ICI*1e-3*exptInfo.fr;
dur = ceil(stim_dur+ici_dur);
%%
stim = zeros(1,size(spikes.raster,2));
stim_ind = repelem(events.eventsOn',ceil(stim_dur));
stim_ind = stim_ind + repmat(1:ceil(stim_dur),[1 length(events.eventsOn)]);
stim(stim_ind) = 1;
%%
predur = 10;
r_ca = dh.rasterize(calcium.npilSubTraces,dur+predur,...
    events.eventsOn-predur);
r_sp = dh.rasterize(spikes.raster,dur+predur,events.eventsOn-predur);
r_in = dh.rasterize(stim,dur+predur,events.eventsOn-predur);
seq_r_ca = dh.rasterize(calcium.npilSubTraces,dur*8+predur,...
    events.eventsOn(1:8:end)-predur);
seq_r_sp = dh.rasterize(spikes.raster,dur*8+predur,...
    events.eventsOn(1:8:end)-predur);
seq_r_in = dh.rasterize(stim,dur*8+predur,...
    events.eventsOn(1:8:end)-predur);
%% plot neurons for seq
seq = 1;
for n = 1 : 40
    clf; hold on
    ph.pltsqz(seq_r_sp(n,:,seq),'LineWidth',2)
    in = seq_r_in(1,:,seq);
    ph.pltsqz(find(in),0*in(in>0),'r*','LineWidth',6);
    title(['neuron ' num2str(n) ' seq ' num2str(seq)])
    hold off; ph.prefs; pause
end
%% plot mean seq
for n = 1 : 40
    clf; hold on
    for s = 1 : 4
        ph.pltsqz(smooth(mean(...
            seq_r_sp(n,:,stimInfo.order==s)...
            ,3)),'LineWidth',2)
    end
    in = seq_r_in(1,:,seq);
    ph.pltsqz(find(in),0*in(in>0),'r*','LineWidth',6);
    legend({'AA','BB','AB','BA','stimulus'})
    hold off; ph.prefs; pause
end
%% plot var, mean of seqs
subplot(1,5,1)
imagesc(var(calcium.npilSubTraces,0,2))
subplot(1,5,2)
ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order==1),3))
subplot(1,5,3)
ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order==2),3))
subplot(1,5,4)
ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order==3),3))
subplot(1,5,5)
ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order==4),3))
%% find neurons with most variable activity
var_n = var(calcium.npilSubTraces,0,2);
n_hi_var = find(var_n > .9e6)';
%% plot those neurons
subplot(1,4,1)
ph.imgsqz(mean(seq_r_ca(n_hi_var,:,stimInfo.order==1),3))
subplot(1,4,2)
ph.imgsqz(mean(seq_r_ca(n_hi_var,:,stimInfo.order==2),3))
subplot(1,4,3)
ph.imgsqz(mean(seq_r_ca(n_hi_var,:,stimInfo.order==3),3))
subplot(1,4,4)
ph.imgsqz(mean(seq_r_ca(n_hi_var,:,stimInfo.order==4),3))
%% plot those neurons individually
for n = 1 : length(n_hi_var)
    clf; hold on;
    ph.pltsqz(mean(seq_r_ca(n_hi_var(n),:,stimInfo.order==1),3),'LineWidth',2)
    ph.pltsqz(mean(seq_r_ca(n_hi_var(n),:,stimInfo.order==2),3),'LineWidth',2)
    ph.pltsqz(mean(seq_r_ca(n_hi_var(n),:,stimInfo.order==3),3),'LineWidth',2)
    ph.pltsqz(mean(seq_r_ca(n_hi_var(n),:,stimInfo.order==4),3),'LineWidth',2)
    in = seq_r_in(1,:,seq);
    ph.pltsqz(find(in),0*in(in>0),'r*','LineWidth',6);
    title(['neuron ' num2str(n_hi_var(n)) ' var ' num2str(var_n(n_hi_var(n)),'%1.2g')])
    legend({'AA','BB','AB','BA','stimulus'})
    hold off; ph.prefs; pause
end
%% variance at baseline
var_n_base = var(calcium.npilSubTraces(:,1:events.eventsOn(1)),0,2);
find(var_n_base > .9e6)';
%% plot var & mean at baseline
subplot(1,2,2)
imagesc(calcium.npilSubTraces(:,1:events.eventsOn(1)))
subplot(1,2,1)
imagesc(var(calcium.npilSubTraces(:,1:events.eventsOn(1)),0,2))
%% show mean response to XXXXXXX vs X/Y
titles = {'AAAAAAAA' 'BBBBBBBB' 'AAAAAAAB' 'BBBBBBBA'};
for s = 1 : 4
    subplot(1,4,s)
    imagesc([mean(mean(seq_r_ca(n_hi_var,1:220,stimInfo.order==s),3),2)...
        mean(mean(seq_r_ca(n_hi_var,221:end,stimInfo.order==s),3),2)])
    title(titles{s})
end
%% 

