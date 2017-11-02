
%%
load('F:HarangData/K070_20171027_artificialGrammar_03.mat')

%% 
% delta F / F0 per trial
f = calcium.npilSubTraces;
len_trial = ceil(stimInfo.t_dur * exptInfo.fr);
len_pretrial = floor(stimInfo.ISI * exptInfo.fr/2);
len_posttrial = len_pretrial;
window_sampling = dh.window_pretrial(events.eventsOn, len_pretrial);
window_f0 = dh.window_paratrial(events.eventsOn, len_pretrial, ...
    len_trial, len_posttrial);
f0 = dh.f0_min_movmean(f, window_sampling, window_f0);
dff0 = (f - f0) ./ f0;
dff0_sm = dh.denoise(dff0);

%% 
% delta F / F0 for all trials


%%
clf; hold on
n = 100;
s = 2;
plot(dff0_sm(n,:)); ph.prefs;
si = events.eventsOn(stimInfo.order==s);
yl = ylim;
for i = 1 : length(si)
    plot([si(i) si(i)+1e-4], [yl(1) yl(2)], 'k--')
end
hold off
clear n s si yl i

%%
dur = ceil((stimInfo.ISI + stimInfo.t_dur) * exptInfo.fr);
r_sp = dh.rasterize(spikes.raster, dur, events.eventsOn);
r_sp_avg = squeeze(mean(r_sp,2));
r_sp_max = squeeze(max(r_sp,2));

%% transitions
trans = dh.find_transitions(stimInfo.order);
trans_uniq = dh.get_unique_transitions(trans);
r_sp_tr = dh.rasterize_by_trans(r_sp, trans, trans_uniq);
%%
bar(squeeze(mean(max(r_sp_tr,[],2),1)),'k')
% bar(squeeze(mean(mean(r_sp_tr,2),1)))
ticks = cell(size(trans_uniq,1),1);
for i = 1 : size(trans_uniq,1)
    s1 = trans_uniq(i,1);
    s2 = trans_uniq(i,2);
    ticks{i} = [num2str(s1) '->' num2str(s2) '; P='...
        num2str(stimInfo.grammar(s2,s1))];
end
ax = gca; ax.XTickLabels = ticks;
xlabel('transitions')
ph.prefs
clear ticks i s1 s2 ax
%% t-test


%% raw spikes

