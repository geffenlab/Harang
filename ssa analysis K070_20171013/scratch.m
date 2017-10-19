
%%
load('F:/HarangData/K070_20171013_SSA_03.mat')
%%
stim_dur = stimInfo.tDur*1e-3*exptInfo.fr;
ici_dur = stimInfo.ICI*1e-3*exptInfo.fr;
dur = ceil(stim_dur+ici_dur);
%%
stim = zeros(1,size(spikes.raster,2));
stim_ind = repelem(events.eventsOn',ceil(stim_dur));
stim_ind = stim_ind + repmat(0:ceil(stim_dur)-1,[1 length(events.eventsOn)]);
stim(stim_ind) = 1;
%%
predur = 50;
post_stim_dur = 50;
post_seq_dur = 200;
r_ca = dh.rasterize(calcium.npilSubTraces,dur+predur+post_stim_dur,...
    events.eventsOn-predur);
r_sp = dh.rasterize(spikes.raster,dur+predur+post_stim_dur,...
    events.eventsOn-predur);
r_in = dh.rasterize(stim,dur+predur+post_stim_dur,...
    events.eventsOn-predur);
seq_r_ca = dh.rasterize(calcium.npilSubTraces,dur*8+predur+post_seq_dur,...
    events.eventsOn(1:8:end)-predur);
seq_r_sp = dh.rasterize(spikes.raster,dur*8+predur+post_seq_dur,...
    events.eventsOn(1:8:end)-predur);
seq_r_in = dh.rasterize(stim,dur*8+predur+post_seq_dur,...
    events.eventsOn(1:8:end)-predur);
%% find neurons with most variable activity
var_n = var(calcium.npilSubTraces,0,2);
% seq_r_ca(n_hi_var(n),:,stimInfo.order==1);
n_hi_var = find(var_n > .5e6)';
%%
seq_r_in_3_4 = seq_r_in;
seq_r_in_3_4(:,:,stimInfo.order==1 | stimInfo.order==2) = 0;
seq_r_in_4 = seq_r_in;
seq_r_in_4(:,:,stimInfo.order==1 | stimInfo.order==2 | ...
    stimInfo.order==3) = 0;
seq_r_in_1_2 = seq_r_in;
seq_r_in_1_2(:,:,stimInfo.order==3 | stimInfo.order==4) = 0;
seq_r_in_2 = seq_r_in;
seq_r_in_2(:,:,stimInfo.order==1 | stimInfo.order==3 | ...
    stimInfo.order==4) = 0;


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
titles = {'AAAAAAAA' 'AAAAAAAB' 'BBBBBBBB' 'BBBBBBBA'};
subplot(1,5,1)
imagesc(var(calcium.npilSubTraces,0,2))
title('var in Ca2+')
for s = 1 : 4
    subplot(1,5,s+1)
    ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order==s),3))
    title(titles{s})
end
%% plot those neurons
for s = 1 : 4
    subplot(1,4,s)
    ph.imgsqz(mean(seq_r_ca(n_hi_var,:,stimInfo.order==s),3))
    title(titles{s})
end
%% plot those neurons individually
for n = 1 : length(n_hi_var)
    clf; hold on;
    hs = zeros(1,5);
    for s = 1 : 4
        data = seq_r_ca(n_hi_var(n),:,stimInfo.order==s);
%         data = seq_r_ca(n,:,stimInfo.order==s);
        y = smooth(mean(data,3),5)';
        std_y = std(data,0,3)/sqrt(length(y));
        x = 1:length(y);
        hs(s) = ph.pltsqz(y,'LineWidth',2);
        color = get(hs(s),'Color');
        patch([x fliplr(x)], [y-std_y, fliplr(y)+std_y], color,...
            'facealpha',0.3,'edgecolor','none')
    end
    in = seq_r_in(1,:,1);
    fig = gca;
    hs(5) = ph.pltsqz(find(in),0*in(in>0),'b*','LineWidth',6);
    title(['neuron ' num2str(n_hi_var(n)) ' var ' num2str(var_n(n_hi_var(n)),'%1.2g')])
    onsets = strfind(in,01);
    for i = 1 : length(onsets)  
        plot([onsets(i)-0.001 onsets(i)+0.001],[0 fig.YLim(2)],'k--')
    end
    plot([])
%     title(['neuron ' num2str(n)])
    legend(hs,{'AA','AB','BB','BA','stimulus'})
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
for s = 1 : 4
    subplot(1,4,s)
    imagesc([mean(mean(seq_r_ca(n_hi_var,1:220,stimInfo.order==s),3),2)...
        mean(mean(seq_r_ca(n_hi_var,221:end,stimInfo.order==s),3),2)])
    title(titles{s})
end
%% plot var, mean of seqs -- 1st half
titles = {'AAAAAAAA' 'BBBBBBBB' 'AAAAAAAB' 'BBBBBBBA'};
subplot(1,5,1)
imagesc(var(calcium.npilSubTraces,0,2))
title('var in Ca2+')
half1 = 1:10;
for s = 1 : 4
    subplot(1,5,s+1)
    ph.imgsqz(mean(seq_r_ca(:,:,stimInfo.order(half1)==s),3))
    title(titles{s})
end
%%
clf; hold on
ph.prefs
ph.pltsqz(dh.split_concat(seq_r_sp(62,:,:)),'LineWidth',2);
ph.pltsqz(dh.split_concat(seq_r_in_3_4));
ph.pltsqz(dh.split_concat(seq_r_in_4));
legend({'spikes (neuron 62)','B...B','B...A'})
%%
clf; hold on
ph.prefs
c = .7;
ph.pltsqz(dh.split_concat(seq_r_sp(37,:,:)),'LineWidth',2);
ph.pltsqz(dh.split_concat(c*seq_r_in_1_2));
ph.pltsqz(dh.split_concat(c*seq_r_in_2));
% ph.pltsqz(dh.split_concat(c*seq_r_in_3_4),'LineWidth',2);
% ph.pltsqz(dh.split_concat(c*seq_r_in_4),'LineWidth',2);
legend({'spikes (neuron 37)','A...A','A...B'})
%%
clf; hold on
ph.prefs
ph.pltsqz(dh.split_concat(seq_r_sp(42,:,:)),'LineWidth',2);
ph.pltsqz(dh.split_concat(seq_r_in_3_4));
ph.pltsqz(dh.split_concat(seq_r_in_4));
legend({'spikes (neuron 42)','B...B','B...A'})
%%
% clf; hold on; ph.prefs
% ph.pltsqz(smooth(mean(r_sp(62,:,8*find(stimInfo.order==1)),3),10),'LineWidth',2);
% ph.pltsqz(smooth(mean(r_sp(62,:,8*find(stimInfo.order==1)+7),3),10),'LineWidth',2);
%%
clf; hold on; ph.prefs;
ph.pltsqz(smooth(mean(seq_r_sp(62,:,stimInfo.order==1),3),1),'LineWidth',2)
ph.pltsqz(.1*seq_r_in(:,:,1),'LineWidth',2)

