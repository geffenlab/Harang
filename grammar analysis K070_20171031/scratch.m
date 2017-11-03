
%%
load('F:HarangData/K070_20171027_artificialGrammar_03.mat')


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

%%
clf; hold on
n = 6;
s = 2;
plot(dff0_sm(n,:));
% plot(f0(n,:));
% plot(f(n,:));
% ph.prefs;
si = events.eventsOn(stimInfo.order==s);
yl = ylim;
% for i = 1 : length(si)
%     plot([si(i) si(i)+1e-4], [yl(1) yl(2)], 'k--')
% end
axis([0 inf -1 2])
hold off
for i = 1 : size(f,1)
    plot(dff0_sm(i,:))
    pause
end
clear n s si yl i

%%
% rasterize df/f0
f_ras = dh.rasterize(dff0_sm,dur,onsets);
%%
s = 3;
r = 1:100;
for n = 1 : size(f,1)
    clf; hold on
    ph.pltsqz(dh.split_concat(f_ras(n,:,r)));
    si = events.eventsOn(stimInfo.order(r)==s);
    yl = ylim;
    for i = 1 : length(si)
        plot([si(i) si(i)+1e-4], [yl(1) yl(2)], 'k--')
    end
    ph.prefs
    pause
end
hold off
clear s r n si yl
%% average in each trial
f_ras_avg = squeeze(mean(f_ras,2));
%%
for n = 1 : size(f,1)
    clf; hold on
    plot(f_ras_avg(n,:))
    si = find(stimInfo.order==1);
    scatter(si,zeros(1,length(si)),'r.')
    si = find(stimInfo.order==2);
    scatter(si,.1*ones(1,length(si)),'c.')
    si = find(stimInfo.order==3);
    scatter(si,.2*ones(1,length(si)),'y.')
    si = find(stimInfo.order==4);
    scatter(si,.3*ones(1,length(si)),'g.')
    hold off
    pause
end
clear n












%%
r_sp = dh.rasterize(spikes.raster, dur, events.eventsOn);
% r_sp_avg = squeeze(mean(r_sp,2));
% r_sp_max = squeeze(max(r_sp,2));

%% transitions
trans = dh.find_transitions(stimInfo.order);
trans_uniq = dh.get_unique_transitions(trans);
% r_sp_tr = dh.rasterize_by_trans(r_sp, trans, trans_uniq);
%% transition index
trans_idx = cell(1,size(trans_uniq,1));
for i = 1 : length(trans_idx)
    trans_idx{i} = strfind(stimInfo.order, trans_uniq(i,:));
end
% empirical P(transition)
trans_prob = cellfun(@length, trans_idx) / size(trans,1);
clear i
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
t_sp = dh.split_concat(r_sp(:,:,t_idx(seq):t_idx(seq)+length(t)-1));
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
clear i
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


