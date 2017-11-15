%%
load('E:/HarangData/K070_20171013_SSA_03.mat')
%% plot stimulus
clf; hold on
plot(1,stimInfo.chordFreqs(:,1)/1e3,'ok')
plot(2,stimInfo.chordFreqs(:,2)/1e3,'ok')
plot(3,stimInfo.chordFreqs(:,3)/1e3,'ok')
plot(4,stimInfo.chordFreqs(:,4)/1e3,'ok')
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
len_trial = ceil(stimInfo.tDur/1e3 * exptInfo.fr);
onsets = events.eventsOn;
n_neur = size(f,1);
dur = ceil((stimInfo.ICI + stimInfo.tDur)/1e3 * exptInfo.fr);
%%
% delta F / F0 mean of past S seconds
samp_dur = 19; % sec
len_samp = floor(samp_dur * dur);
window_samp = dh.make_windows(onsets-len_samp, len_samp);
f0_vals = dh.f0_mean(f, window_samp);
len_post = floor(stimInfo.ICI/1e3 * exptInfo.fr);
len_pre = len_post;
window_f0 = dh.make_windows(onsets-len_pre,len_pre+len_trial+len_post);
f0 = dh.f0_apply(f, f0_vals, window_f0);
dff0 = (f - f0) ./ f0;
dff0_sm = dh.denoise(dff0);
clear samp_dur len_samp window_samp f0_vals len_post window_f0
%% rasterize df/f0
f_ras = dh.rasterize(dff0_sm,dur,onsets);
%% average in each trial
f_ras_avg = squeeze(mean(f_ras,2));
%%
stim_per_seq = 8;
dur_seq = stim_per_seq * dur;
order = repelem(stimInfo.order',1,stim_per_seq);
dur_seq_pre = dur;
dur_seq_post = dur;
dur_tot = dur_seq + dur_seq_pre + dur_seq_post;
f_ras_seq = dh.rasterize(dff0_sm, dur_tot,...
    onsets(1:stim_per_seq:length(order)) - dur_seq_pre);
% f_ras_seq = dh.rasterize(dff0_sm, dur_seq, onsets(1:stim_per_seq:length(order)));
%%
n = 62;
for i = 1 : length(stimInfo.order)
    clf; hold on
    ph.pltsqz(...
        linspace(-1*dur_seq_pre/exptInfo.fr,...
        dur_tot/exptInfo.fr, dur_tot),...
        f_ras_seq(n,:,i));
    %     axis([-1*dur_seq_pre/exptInfo.fr ...
    %         (dur_seq+dur_seq_post)/exptInfo.fr -2 4])
    axis([-inf inf -2 4])
    %     ph.pltsqz(linspace(0,dur_seq/exptInfo.fr, dur_seq), f_ras_seq(n,:,i));
    %     axis([-1 dur_seq/exptInfo.fr -2 4])
    ax = gca;
    for j = 1 : stim_per_seq
        idx = j-1;
        plot([idx idx+1e-3], [ax.YLim(1) ax.YLim(2)], '--k')
    end
    title(['@' num2str(i) ' stim ' num2str(stimInfo.order(i))])
    ph.prefs; hold off
    waitforbuttonpress
end
clear n i ax j idx
%%
for n = 1 : size(f,1)
    clf; hold on
    for s = 1 : length(stimInfo.index)
        x = linspace(-1*dur_seq_pre/exptInfo.fr, dur_tot/exptInfo.fr,...
            dur_tot);
        y_ = f_ras_seq(n,:,stimInfo.order==s);
        y = squeeze(mean(y_,3));
        e = squeeze(std(y_,0,3)/sqrt(size(y_,3)));
%         plot(x,y,'LineWidth', 2);
        ph.error_shade(x,y,e,'r');
    end
    axis([-inf inf -2 2])
    ax = gca;
    for j = 1 : stim_per_seq
        idx = j-1;
        plot([idx idx+1e-3], [ax.YLim(1) ax.YLim(2)], '--k')
    end
    legend({'1','2','3','4'})
    title(['n=' num2str(n)])
    ph.prefs; hold off
    waitforbuttonpress
end
clear n s ax idx

%% saved

%%







