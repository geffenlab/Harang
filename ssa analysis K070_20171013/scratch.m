%%
load('F:/HarangData/K070_20171013_SSA_03.mat')
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
ns = [29 32 37 40 46 51 52 54 62 64 66 77 82 84 94 105 114 116 123 125 ...
    128 134 140 143 149 152 161 162 168 178 186 190 207 225 242 256 260 ...
    266 274];
for n = ns %1 : size(f,1)
    clf; hold on
    for s = 1 : length(stimInfo.index)
        x = linspace(-1*dur_seq_pre/exptInfo.fr, dur_tot/exptInfo.fr,...
            dur_tot);
        y_ = f_ras_seq(n,:,stimInfo.order==s);
        y = squeeze(mean(y_,3));
        e = squeeze(std(y_,0,3)/sqrt(size(y_,3)));
        switch s
            case 1; c = 'r';
            case 2; c = 'b';
            case 3; c = 'g';
            case 4; c = 'c';
            otherwise; c = 'k';
        end
        ph.error_shade(x,y,e,c,'LineWidth',2);
    end
    xlabel('time (s)'); ylabel('\DeltaF/F0')
    axis([-inf inf -1 1.5])
    ax = gca;
    for j = 1 : stim_per_seq
        idx = j-1;
        plot([idx idx+1e-3], [ax.YLim(1) ax.YLim(2)], '--k')
    end
%     legend({'1','2','3','4'})
    title(['n=' num2str(n)])
    ph.prefs; hold off
%     filename = ['dff0 n' num2str(n)];
%     savefig(filename);
%     print(filename, '-dpng');
    waitforbuttonpress
end
clear ns n s x y_ y e c ax j idx

%% neuron 37 response
n = 37;
clf; hold on
for s = 1 : length(stimInfo.index)
    x = linspace(-1*dur_seq_pre/exptInfo.fr, dur_tot/exptInfo.fr,...
        dur_tot);
    y_ = f_ras_seq(n,:,stimInfo.order==s);
    y = squeeze(mean(y_,3));
    e = squeeze(std(y_,0,3)/sqrt(size(y_,3)));
    switch s
        case 1; c = 'r';
        case 2; c = 'b';
        case 3; c = 'g';
        case 4; c = 'c';
        otherwise; c = 'k';
    end
    ph.error_shade(x,y,e,c,'LineWidth',2);
end
xlabel('time (s)'); ylabel('\DeltaF/F0')
axis([-inf inf -1 1.5])
ax = gca;
for j = 1 : stim_per_seq
    idx = j-1;
    plot([idx idx+1e-3], [ax.YLim(1) ax.YLim(2)], '--k')
end
ph.prefs
clear n x y_ y e s c ax j idx

%% glm
% logistic
n = 37;
s = 1;
% generate data
x = linspace(-1*dur_seq_pre/exptInfo.fr, dur_tot/exptInfo.fr,...
    dur_tot);
y_ = f_ras_seq(n,:,stimInfo.order==s);
y = squeeze(mean(y_,3));
% plot
clf; hold on
idx = dur_seq_pre + 1 : dur_seq_pre + dur;
plot(x(idx), y(idx), '.')
ph.prefs
% glm fit
b = glmfit(x(idx), [y(idx)], 'binomial', 'logit');
v = glmval(b, x(idx), 'logit');
% plot fit
plot(x, v, 'r-')

clear n s x y_ y idx b v

%% saved

%%







