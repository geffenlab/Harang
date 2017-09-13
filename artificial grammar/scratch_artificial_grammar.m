
%% load data
load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')
num_neur = size(spikes.raster, 1)
num_stim = length(stimInfo.order)
num_words = length(stimInfo.words)
frames_stim = ceil(exptInfo.fr * (stimInfo.t_dur+stimInfo.ISI))
frames_ahead = ceil(exptInfo.fr * 0.1) % ~0.1s
frames = frames_stim + frames_ahead;

%% separate spikes by stimulus
spikes_stim = zeros(num_neur, num_stim, frames);
for s = 1 : num_stim
    onset = events.eventsOn(s) - frames_ahead;
    offset = onset + frames_stim + frames_ahead - 1;
    spikes_stim(:, s, :) = spikes.raster(:, onset : offset);
end

%% average spikes for each word
spikes_word = zeros(num_neur, num_words, frames);
for w = 1 : num_words
    spikes_word(:, w, :) = mean(...
        spikes_stim(:, find(stimInfo.order==w), :), ...
        2);
end

%% plot spikes for each word
figure
max_spikes_word1 = max(spikes_word(:, 1, :), [], 3);
[~, max_spike_indices] = sort(max_spikes_word1, 'descend');
for w = 1 : num_words
    subplot(1, num_words, w)
    imagesc(squeeze(spikes_word(max_spike_indices, w, :)))
    title(['word ' num2str(w)]);
    xlabel('frame'); ylabel('neuron')
end

%% separate spikes by transition
spikes_trans = zeros(num_neur, num_words * num_words, frames);
trans_all = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
trans_unique = unique(trans_all, 'rows');
num_trans = zeros(1, length(trans_unique));
for t = 1 : num_words * num_words
    trans = trans_unique(t,:);
    indices_trans = find(trans_all(:,1)==trans(1) & ...
        trans_all(:,2)==trans(2));
    num_trans(t) = sum(indices_trans>0);
%     spikes_trans(:,t,:) = mean(spikes_stim(:,indices_trans(1:50),:), 2);
    spikes_trans(:,t,:) = mean(spikes_stim(:,indices_trans,:), 2);
end

%% plot spikes by transition
for t = 1 : num_words * num_words
    subplot(num_words, num_words, t)
    imagesc(squeeze(spikes_trans(max_spike_indices(2:100), t, :)))
    caxis([0 0.3]); xlabel('frames'); ylabel('neurons');
    title([num2str(trans_unique(t,1)) '->' ...
        num2str(trans_unique(t,2))]);
    colorbar
end