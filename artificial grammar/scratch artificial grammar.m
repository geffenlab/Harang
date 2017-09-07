
%% outputs / spikes
% normalize spikes
% spikes_norm = log(spikes.raster) + 3;
spikes_norm = spikes.raster;

%% inputs / stimuli
% event on
stimuli = zeros(3, length(spikes_norm));
for i = 1 : length(events.eventsOn)
    stimuli(stimInfo.order(i), events.eventsOn(i)) = 1;
end

%% find correlation bt neurons
corr_neurons = spikes_norm * spikes_norm';

%% find correlation bt stimuli and neurons
corr_stim_neur = stimuli * spikes_norm';

%% find mean Ca response
frames_per_sample = 11;
frames_ahead = 1;
ca_word_stim = zeros(316, 2000, frames_per_sample + frames_ahead);
for n = 1 : 316
    for s = 1 : 2000
        onset = events.eventsOn(s)-1;
        ca_word_stim(n, s, :) = calcium.npilSubTraces(n, ...
            onset - frames_ahead : onset + frames_per_sample - 1);
    end
end
ca_word = zeros(316, 3, frames_per_sample + frames_ahead);
for n = 1 : 316
    for w = 1 : 3
        ca_word(n, w, :) = mean(...
            ca_word_stim(n, events.eventsOn==w, :), 2);
    end
end
