
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
%% calculate SSA for neuron 62, sequence 1
n = 62;
seq_avg_62 = mean(seq_r_sp(n,:,stimInfo.order==1),3);
seq_max_avg_62 = dh.max_in_windows(seq_avg_62, ceil(stim_dur+ici_dur),...
    strfind(seq_r_in(:,:,1),[1 1 1 1 1]));
ssa_avg_62 = seq_max_avg_62(1)/seq_max_avg_62(end);
%% calculate all SSA
% to stim 1
seq = 1;
seq_avg_stim1 = mean(seq_r_sp(:,:,stimInfo.order==seq),3);
seq_max_avg_stim1 = dh.max_in_windows(seq_avg_stim1, ceil(stim_dur+ici_dur),...
    strfind(seq_r_in(:,:,seq),[1 1 1 1 1]));
ssa_stim1 = seq_max_avg_stim1(:,1) ./ seq_max_avg_stim1(:,end);
% to stim 3
seq_avg_stim3 = mean(seq_r_sp(:,:,stimInfo.order==seq),3);
seq_max_avg_stim3 = dh.max_in_windows(seq_avg_stim3, ceil(stim_dur+ici_dur),...
    strfind(seq_r_in(:,:,seq),[1 1 1 1 1]));
ssa_stim3 = seq_max_avg_stim3(:,1) ./ seq_max_avg_stim3(:,end);




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
        data = seq_r_sp(n_hi_var(n),:,stimInfo.order==s);
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
ph.pltsqz(dh.split_concat(seq_r_sp(62,:,:)));
ph.pltsqz(dh.split_concat(seq_r_in_3_4));
ph.pltsqz(dh.split_concat(seq_r_in_4));
ph.pltsqz(dh.split_concat(seq_r_in_1_2));
legend({'spikes (neuron 62)','B...B','B...A'})
ph.prefs
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
ph.pltsqz(mean(seq_r_sp(62,:,stimInfo.order==1),3),'LineWidth',2);
ph.pltsqz(.1*seq_r_in(:,:,1),'LineWidth',2);




%%
n = 62; seq = 1;
clf; hold on;
pad = 30;
xs = strfind(seq_r_in(:,:,seq),[1 1 1 1 1]);
ph.pltsqz(smooth(mean(...
    seq_r_sp(n,xs(1)-pad:xs(1)+pad+ceil(stim_dur+ici_dur),...
    stimInfo.order==1),3),10));
ph.pltsqz(smooth(mean(...
    seq_r_sp(n,xs(end)-pad:xs(end)+pad+ceil(stim_dur+ici_dur),...
    stimInfo.order==1),3),10));
fig = gca;
plot([pad+1 pad+1.0001],[0 fig.YLim(2)],'k--')
legend({'first response','last response','stimulus onset'})
title('SSA to Sequence AAAAAAAA (neuron 62)')
annotation('textbox',[0.15 0.4 .3 .5],'String',...
    ['SSA = peak_{first}/peak_{last} = ' num2str(ssa_avg_62)],...
    'FontSize',16,'LineStyle','none','FitBoxToText','on');
fig = gca;
plot([pad+1 pad+1.0001],[0 fig.YLim(2)],'k--')
plot([pad+ceil(stim_dur)+1 pad+ceil(stim_dur)+1.0001],[0 fig.YLim(2)],'k--')
plot([pad+ceil(stim_dur)+dur+1 pad+ceil(stim_dur)+dur+1.0001],[0 fig.YLim(2)],'k--')
ph.prefs;
%% just one instance neuron 62 - one-by-one
n = 62;
for i = 1 : 80
    clf; hold on
    ph.pltsqz(smooth(seq_r_sp(n,:,i),5));
    xs = strfind(seq_r_in(:,:,i),[1 1 1 1 1]);
    fig = gca;
    for j = 1 : 8
        plot([xs(j) xs(j)+0.0001],[0 fig.YLim(2)],'k--')
    end
    title(['seq ' num2str(i)])
    ph.prefs % 24 31 44 72 79
    pause
end
%% neuron 62 plot at seq 31
n = 62; seq = 31;
subplot(2,1,1)
ph.pltsqz(seq_r_sp(n,:,seq));
hold on
xs = strfind(seq_r_in(:,:,seq),[1 1 1 1 1]);
fig = gca;
for j = 1 : 8
    plot([xs(j) xs(j)+0.0001],[0 fig.YLim(2)],'k--')
end
xlabel('t'); ylabel('spike')
title(['Neuron ' num2str(n) ' spikes at ' num2str(seq)...
    '^{st} sequence ("AAAAAAAA")'])
ph.prefs; hold off
%% mean spikes to sequence AAAAAAAA for neuron 62
n = 62; stim = 1;
subplot(2,1,2)
ph.pltsqz(smooth(mean(seq_r_sp(n,:,stimInfo.order==stim),3),10),'LineWidth',2);
hold on
ph.pltsqz(find(seq_r_in(:,:,stim)),0*seq_r_in(:,seq_r_in(:,:,stim)>0,stim),...
    'r*','LineWidth',2);
xs = strfind(seq_r_in(:,:,stim),[1 1 1 1 1]);
fig = gca;
for i = 1 : 8
    plot([xs(i) xs(i)+0.0001],[0 fig.YLim(2)],'k--')
end
title('averaged for all instances')
ph.prefs; hold off
%% neuron 62 ssa plot - seq 31 + raster of all other sequences
% clf; hold on
subplot(2,1,1);
n = 62; seq = 31;
pad = 20;
xs = strfind(seq_r_in(:,:,seq),[1 1 1 1 1]);
pts1 = xs(1)-pad:xs(1)+dur+pad;
pts8 = xs(end)-pad:xs(end)+dur+pad;
ph.pltsqz(seq_r_sp(n,pts1,seq)); hold on
ph.pltsqz(seq_r_sp(n,pts8,seq));
plot([pad+1 pad+1.0001],[0 1],'k--')
plot([pad+ceil(stim_dur)+1 pad+ceil(stim_dur)+1.0001],[0 1],'k--')
plot([pad+ceil(stim_dur)+dur+1 pad+ceil(stim_dur)+dur+1.0001],[0 1],'k--')
legend({'first response','last (8^{th}) response','stimulus on/offset'})
title(['Neuron ' num2str(n) ' SSA at ' num2str(seq) ...
    '^{st} sequence ("AAAAAAAA")'])
% axis([0 length(pts1) -1 1])
ph.prefs; hold off

subplot(2,1,2);
ph.pltsqz(smooth(mean(...
    seq_r_sp(n,xs(1)-pad:xs(1)+pad+ceil(stim_dur+ici_dur),...
    stimInfo.order==1),3),10));
hold on;
ph.pltsqz(smooth(mean(...
    seq_r_sp(n,xs(end)-pad:xs(end)+pad+ceil(stim_dur+ici_dur),...
    stimInfo.order==1),3),10));
fig = gca;
plot([pad+1 pad+1.0001],[0 fig.YLim(2)],'k--');
% legend({'first response','last response','stimulus onset'})
% title('SSA to Sequence AAAAAAAA (neuron 62)')
title('averaged & smoothed over all instances')
annotation('textbox',[0.15 .4 .3 .5],'String',...
    ['SSA = peak_{first}/peak_{last} = ' num2str(ssa_avg_62)],...
    'FontSize',16,'LineStyle','none','FitBoxToText','on');
fig = gca;
plot([pad+1 pad+1.0001],[0 fig.YLim(2)],'k--')
plot([pad+ceil(stim_dur)+1 pad+ceil(stim_dur)+1.0001],[0 fig.YLim(2)],'k--')
plot([pad+ceil(stim_dur)+dur+1 pad+ceil(stim_dur)+dur+1.0001],[0 fig.YLim(2)],'k--')
ph.prefs; hold off

% all seq
% fr_1 = squeeze(seq_r_sp(1:50,pts1,seq));
% fr_2 = squeeze(seq_r_sp(1:50,pts8,seq));
% bound = -.5; step = bound/50;%size(r_ca,1);
% plot([0 length(pts1)], [bound bound],'k--')
% imagesc(1:length(pts1),step:step:bound,fr_1,'AlphaData',~~fr_1);
% imagesc(1:length(pts1),step+bound:step:bound*2,fr_2,'AlphaData',~~fr_2);
% ph.prefs; colorbar; colormap(flipud(bone));
%% ^^^^^^^^^^ MAKE THIS A FUNCTION


%% NSF figures
%% neuron 62 plot, seq AAAAAAAA
n = 62; seq = 31; clf
pts = 45:350;
ph.pltsqz(smooth(mean(seq_r_sp(n,pts,stimInfo.order==stim),3),10),'b',...
    'LineWidth',2);
hold on
ph.pltsqz(find(seq_r_in(:,pts,stim)),...
    0*seq_r_in(:,seq_r_in(:,pts,stim)>0,stim),...
    'r*','LineWidth',2);
xs = strfind(seq_r_in(:,pts,stim),[1 1 1 1 1]);
fig = gca;
for i = 1 : 8
    plot([xs(i) xs(i)+0.0001],[0 fig.YLim(2)],'k--')
end
xlabel('time'); ylabel('average spike')
ax = gca;
ph.prefs; hold off
%% 
clf
activity = mean(mean(seq_r_sp(:,:,stimInfo.order==stim),3),2);
activity = activity / max(activity);
ph.plot_ca_img(spatialInfo,activity)
n62_patch = spatialInfo.ROIs{62};
ax = gca;
ax.XTick = []; ax.XTickLabel = {};
ax.YTick = []; ax.YTickLabel = {};
patch(n62_patch(:,1),n62_patch(:,2),'b','LineStyle','none')

%%
pts = 1:size(seq_r_in,2);
for i = 1 : length(pts)
    activity = seq_r_sp(:,i,31);
    activity = smooth(activity,10);
    activity = activity/max(activity);
    ph.plot_ca_img(spatialInfo,activity)
    title(['t=' num2str(i)])
    pause(0.01)
end

%% plot
subplot(2,1,1);
hist(ssa_stim1, 30); ph.prefs;
title('SSA for neurons firing to sequence AAAAAAAA')
ylabel('# neurons');
subplot(2,1,2);
title('SSA for neurons firing to sequence BBBBBBBB')
hist(ssa_stim3, 30); ph.prefs;
xlabel('SSA'); ylabel('# neurons');


