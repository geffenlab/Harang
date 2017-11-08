
%%
load('E:/HarangData/K070_20171027_artificialGrammar_03.mat')

%% plot stimulus
clf; hold on
plot(1,stimInfo.chordTones(:,1)/1e3,'ok')
plot(2,stimInfo.chordTones(:,2)/1e3,'ok')
plot(3,stimInfo.chordTones(:,3)/1e3,'ok')
plot(4,stimInfo.chordTones(:,4)/1e3,'ok')
ph.prefs; hold off
axis([0 5 0 60])
xlabel('chord index'); ylabel('frequency (kHz)')


%%
f = calcium.npilSubTraces;
% add minimum to f
f_min = min(f,[],2);
if max(f_min) > 0
    error('There is a neuron with f_min > 0');
end
f = f - repmat(f_min,[1 size(f,2)]);
len_trial = ceil(stimInfo.t_dur * exptInfo.fr);
onsets = events.eventsOn;
n_neur = size(f,1);
dur = ceil((stimInfo.ISI + stimInfo.t_dur) * exptInfo.fr);








%% 
% delta F / F0 per trial
len_pre = floor(stimInfo.ISI * exptInfo.fr/2);
len_post = len_pre;
window_samp = dh.window_pretrial(onsets, len_pre);
f0_vals = dh.f0_min_movmean(f, window_samp);
window_f0 = dh.window_paratrial(onsets, len_pre, len_trial, len_post);
f0 = dh.f0_apply(f, f0_vals, window_f0);
dff0 = (f - f0) ./ f0;
dff0_sm = dh.denoise(dff0);
%% 
% delta F / F0 for all trials
window_samp = dh.window_pretrial(onsets(1), onsets(1)-1);
f0_vals = dh.f0_min_movmean(f, window_samp);
window_f0 = dh.window_paratrial(onsets, len_pre, len_trial, len_post);
f0 = dh.f0_apply(f, f0_vals, window_f0);
dff0 = (f - f0) ./ f0;
dff0_sm = dh.denoise(dff0);
%%
% delta F / F0 mean of past S seconds
samp_dur = 19.9; % sec
len_samp = floor(samp_dur * exptInfo.fr);
window_samp = dh.window_pretrial(onsets, len_samp);
% f0_vals = dh.f0_min_movmean(f, window_samp);
f0_vals = dh.f0_mean(f, window_samp);
len_post = floor((stimInfo.t_dur + stimInfo.ISI) * exptInfo.fr);
window_f0 = dh.window_paratrial(onsets, 0, len_trial, len_post);
f0 = dh.f0_apply(f, f0_vals, window_f0);
dff0 = (f - f0) ./ f0;
dff0_sm = dh.denoise(dff0);
clear samp_dur len_samp window_samp f0_vals len_post window_f0

%% plot sample F
clf
n = 2;
plot(f(n,:));
hold on; plot([0 6.5e4], [0 0], '--r')
ph.prefs; axis([0 6.5e4 -inf inf]); ylabel('F'); xlabel('t'); title('neuron 2')
%% plot sample F0
plot(f0(n,:),'r')
ph.prefs
%% plot sample DF/F0
clf
plot(dff0(n,:))
ph.prefs; xlabel('t'); ylabel('\DeltaF/F0'); title(['neuron ' num2str(n)])
%% plot sample smoothed DF/F0
hold on;
plot(dff0_sm(n,:),'r')
ph.prefs; hold off




%% transitions
trans = dh.find_transitions(stimInfo.order);
trans_uniq = dh.get_unique_transitions(trans);
%% transition index
trans_idx = cell(1,size(trans_uniq,1));
for i = 1 : length(trans_idx)
    trans_idx{i} = strfind(stimInfo.order, trans_uniq(i,:));
end
% empirical P(transition)
trans_prob = cellfun(@length, trans_idx) / size(trans,1);
clear i



%% rasterize df/f0
f_ras = dh.rasterize(dff0_sm,dur,onsets);
%% average in each trial
f_ras_avg = squeeze(mean(f_ras,2));



%% correlation with chords
subplot(2,2,1);
hist(corr(stimInfo.order'==1,f_ras_avg'))
title('chord 1'); ph.prefs
subplot(2,2,2);
hist(corr(stimInfo.order'==2,f_ras_avg'))
title('chord 2'); ph.prefs
subplot(2,2,3);
hist(corr(stimInfo.order'==3,f_ras_avg'))
title('chord 3'); ph.prefs
xlabel('correlation'); ylabel('# neurons')
subplot(2,2,4);
hist(corr(stimInfo.order'==4,f_ras_avg'))
title('chord 4'); ph.prefs;

%% which neurons most correlated?
thresh = 0.07;
r1 = corr(stimInfo.order'==1,f_ras_avg');
r2 = corr(stimInfo.order'==2,f_ras_avg');
r3 = corr(stimInfo.order'==3,f_ras_avg');
r4 = corr(stimInfo.order'==4,f_ras_avg');
n_corr = {find(r1>thresh), find(r2>thresh), find(r3>thresh), find(r4>thresh)};
%% plot average response to stimuli
dur_s = dur / exptInfo.fr;
for s = 1 : length(stimInfo.words)
    ni = 1:n_neur;%n_corr{s};
    si = stimInfo.order==s;
    subplot(2,2,s)
%     x = linspace(0,dur_s,dur);
%     y = mean(f_ras(ni, :, si), 3);
%     se = std(f_ras(ni, :, si), 0, 3) / length(si);
%     ph.error_shade(x,y,se,'b')
    ph.pltsqz(linspace(0,dur_s,dur), mean(f_ras(ni, :, si), 3)');
    title(['word ' num2str(s) ' (n=' num2str(length(ni)) ')'])
    ph.prefs
    axis([0 dur_s 0 0.2])
end
subplot(2,2,3); xlabel('sec'); ylabel('avg \DeltaF/F0')
clear s ni si
%% plot average response with individual response for all neurons
for n = 1 : n_neur
clf3


for i = 1 : length(stimInfo.words)
    subplot(2,2,i); hold on
    si = find(stimInfo.order==i);
    x = linspace(0,dur_s,dur);
    for j = 1 : length(si)
        plot(x,squeeze(f_ras(n,:,si(j))));
    end
    plot(x,squeeze(mean(f_ras(n,:,si),3)),'k',...
       'LineWidth',2)
    title(['word ' num2str(i)])
    ph.prefs; hold off
    axis([0 dur_s -1 2])
end
subplot(2,2,3); xlabel('sec'); ylabel('avg \DeltaF/F0')
disp(['n ' num2str(n)])
pause(0.1)
end
clear n i j x

%% get average response for transitions
for t = 1 : size(trans_uniq,1)
    ni = union(n_corr{trans_uniq(t,1)},n_corr{trans_uniq(t,2)});
    si = trans_idx{t};
    x = [mean(f_ras(ni,:,si),3) mean(f_ras(ni,:,si+1),3)];
    subplot(3,3,t)
    ph.pltsqz(linspace(0,2*dur_s,2*dur),x');
    title([num2str(trans_uniq(t,1)) '->' num2str(trans_uniq(t,2))...
        ' P=' num2str(trans_prob(t),'%3.2f')])
    axis([0 2*dur_s -0.1 0.4])
    ph.prefs
end
subplot(3,3,7); xlabel('sec'); ylabel('\DeltaF/F0')
clear t ni si x
%% get average difference is response for transitions
dt_avg = zeros(1,size(trans_uniq,1));
dt_se = zeros(size(dt_avg));
for t = 1 : size(trans_uniq,1)
    si = trans_idx{t};
    d = (f_ras_avg(:,si+1) - f_ras_avg(:,si)) / f_ras_avg(:,si) * 100;
    dt_avg(t) = mean(d(:));
    dt_se(t) = std(d(:)) / sqrt(length(d(:)));
end
clear t si d
%% plot
clf; hold on
bar(dt_avg,'w')
for i = 1 : size(trans_uniq,1)
    plot([i i+0.001], [dt_avg(i)-dt_se(i) dt_avg(i)+dt_se(i)], 'k')
end
xlabel('transitions'); ylabel('percent change')
axis([0 10 -1 1])
ph.prefs
clear i

%% 
y = tsne(f_ras_avg');
%% plot tsne
clf
plot(y(:,1),y(:,2),'.')
ph.prefs
%% plot tsne by stimulus
clf; hold on
s = 1;
plot(y(stimInfo.order==s,1),y(stimInfo.order==s,2),'.')
s = 2;
plot(y(stimInfo.order==s,1),y(stimInfo.order==s,2),'.')
s = 3;
plot(y(stimInfo.order==s,1),y(stimInfo.order==s,2),'.')
s = 4;
plot(y(stimInfo.order==s,1),y(stimInfo.order==s,2),'.')
legend({'1','2','3','4'})
hold off; ph.prefs; colorbar
clear s
%% plot tsne by time
clf
t = 1 : size(y,1);
scatter(y(t,1),y(t,2),32,t,'.')
colorbar
ph.prefs
%% plot tsne by transitions
clf; hold on
for i = 1 : length(trans_idx)
    idx = trans_idx{i};
    scatter(y(idx,1),y(idx,2),40,i*ones(1,length(idx)),'.');
end
colormap jet
tx = cell(size(trans_uniq,1),1);
for i = 1 : size(trans_uniq,1)
    s1 = trans_uniq(i,1);
    s2 = trans_uniq(i,2);
    tx{i} = [num2str(s1) '->' num2str(s2) '; P='...
        num2str(stimInfo.grammar(s2,s1))];
end
ph.prefs; hold off; colorbar
legend(tx)
clear i idx tx
%% find peaks of df/f0





%% plot some sample neurons

s = [3 4 3 4 3 4 3 4 3 4];
sx = 1;
si = strfind(stimInfo.order, s);
ni = n_corr{3};%union(n_corr{s(1)},n_corr{s(2)});
clf
plot(f_ras_avg(ni,si(sx):si(sx)+length(s)-1)')
ylabel('\DeltaF/F0'); xlabel('trials'); title(['seq: ' num2str(s)])
ph.prefs; axis([0 length(s)+1 -inf inf])
clear si ni

%% plot neurons
s = 2;
os = find(stimInfo.order==s);
for n = 1 : size(f,1)
    clf; hold on
    plot(f_ras_avg(n,:))
    ph.prefs; title(['n = ' num2str(n)])
    ax = gca;
    for i = 1 : length(os)
        plot([os(i) os(i)+1e-3],[ax.YLim(1) ax.YLim(2)],'--r')
    end
    hold off
    waitforbuttonpress
end
clear n s os ax





%% plot neurons adapting
n = 51:60;
s = [1 3];
si = strfind(stimInfo.order,s);
si = si(:);
clf
plot(si,f_ras_avg(n,si),'-*')
ph.prefs; axis([si(1) si(end) -1 2])
clear n s si st
%% compare responses with different preceding stimulus
% trans_resp_avg = zeros(size(f,1), size(trans_uniq,1));
trans_resp_avg = zeros(1, size(trans_uniq,1));
trans_resp_se = zeros(size(trans_resp_avg));
for t = 1 : size(trans_resp_avg,2)
    idx = trans_idx{t} + 1;
    x = f_ras_avg(:,idx);
    trans_resp_avg(1,t) = mean(x(:));
    trans_resp_se(1,t) = std(x(:)) / sqrt(length(idx)*size(f_ras_avg,1));
end
clear t i x
%% plot different responses
errorbar(1:length(trans_resp_avg),trans_resp_avg,trans_resp_se,'*',...
    'LineWidth',2)
tx = cell(size(trans_uniq,1),1);
for i = 1 : size(trans_uniq,1)
    tx{i} = num2str(trans_prob(i),'%3.2f');
end
ax = gca; ax.XTick = 1:9; ax.XTickLabels = tx;
xlabel('transitions')
ph.prefs
%%
errorbar(trans_prob,trans_resp_avg,trans_resp_se,'*',...
    'LineWidth',2)
xlabel('transition probabilities')
ph.prefs







%%
r_sp = dh.rasterize(spikes.raster, dur, events.eventsOn);
% r_sp_avg = squeeze(mean(r_sp,2));
% r_sp_max = squeeze(max(r_sp,2));
%% spikes averaged over each trial
% max spike, sum spike
trans_sp_avg = cell(1, size(trans_uniq,1));
trans_sp_max = cell(1, size(trans_uniq,1));
trans_sp_sum = cell(1, size(trans_uniq,1));
for i = 1 : length(trans_sp_avg)
    idx = trans_idx{i};
    trans_sp_avg{i} = squeeze(mean(r_sp(:,:,idx), 2));
    trans_sp_max{i} = squeeze(max(r_sp(:,:,idx), [], 2));
    trans_sp_sum{i} = squeeze(sum(r_sp(:,:,idx), 2));
end
clear i idx
%% find some transitions
% t = [3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4];
% t = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
t = [3 1];
n = 1;
seq = 4;
t_idx = strfind(stimInfo.order, t);
t_sp = dh.split_concat(f_ras_avg(:,:,t_idx(seq):t_idx(seq)+length(t)-1));
% t_sp = dh.split_concat(r_sp(:,:,t_idx(seq):t_idx(seq)+length(t)-1));
clf; hold on
plot(t_sp(n,:)','k')
ax = gca;
for i = 1 : length(t)
    on = dur * (i-1) + 1;
    color = '--r';
    if mod(i,2) == 0
        color = '--b';
    end
    plot([on on+1e-4], [0 ax.YLim(2)], color)
end
ph.prefs; hold off
clear t n seq t_idx t_sp ax i on color
%%
ph.imgsqz(mean(dh.rasterize(t_sp,dur,1:dur:size(t_sp,2)),2))













%%
clf; hold on
sp_tr = squeeze(max(r_sp_tr,[],2));
% sp_tr = squeeze(sum(r_sp_tr,2));
% sp_tr = squeeze(mean(r_sp_tr,2));
bar(mean(sp_tr,1),'w','LineWidth',2)
title('max spike burst averaged over trials')
for i = 1 : size(trans_uniq,1)
    sp = sp_tr(:,i);
    sp_m = mean(sp);
    sp_se = std(sp) / sqrt(length(sp));
    plot([i i+0.001], [sp_m-sp_se sp_m+sp_se], 'k')
end
tx = cell(size(trans_uniq,1),1);
for i = 1 : size(trans_uniq,1)
    s1 = trans_uniq(i,1);
    s2 = trans_uniq(i,2);
    tx{i} = [num2str(s1) '->' num2str(s2) '; P='...
        num2str(stimInfo.grammar(s2,s1))];
end
ax = gca; ax.XTick = 1:9; ax.XTickLabels = tx;
xlabel('transitions')
ph.prefs; hold off
clear tx i s1 s2 ax sp sp_m sp_se

%%
n_trans_u = size(trans_uniq,1);
n_trans = size(trans,1);
prob_trans = zeros(1,n_trans_u); % 1 X n_trans_u
for i = 1 : n_trans_u
    prob_trans(i) = length(find(trans(:,1)==trans_uniq(i,1) & ...
        trans(:,2)==trans_uniq(i,2))) / n_trans;
end
sp = squeeze(max(r_sp_tr,[],2)); % n_neurons X n_trans_u
sp_m = mean(sp,1);
sp_se = std(sp,0,1) / sqrt(size(sp,1));
% sp_m = mean(squeeze(mean(r_sp_tr,2)),1);
% scatter(prob_trans,sp_m,500,'.')
errorbar(prob_trans, sp_m, sp_se,'*','LineWidth',2)
axis([0 0.4 0 0.05])
xlabel('P(transition)'); ylabel('max spike rate (avg over neurons)')
n_neurons = size(f,1);
sp_ = sp';
[rho, pval] = corr(repmat(prob_trans',[n_neurons 1]),1./sp_(:));
title([])
ph.prefs; hold off
clear n_trans_u n_trans i
%% plot probability versus max spike burst
clf; hold on
for i = 1 : size(stimInfo.grammar,1)
    s2_i = find(trans_uniq(:,1) == i);
    s2 = trans_uniq(s2_i,2);
    prob = stimInfo.grammar(s2,i);
    scatter(prob,sp_m(s2_i),500,'.')
end
axis([0 1.1 0 0.05])
ph.prefs; hold off
xlabel('P(transition)'); ylabel('max spike rate')
clear sp_m i s2_i s2 prob
%% t-test


%% raw spikes


