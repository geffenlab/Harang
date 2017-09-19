
%% load data
load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')
num_neur = size(spikes.raster, 1);
num_stim = length(stimInfo.order);
num_words = length(stimInfo.words);
% frames_stim = ceil(exptInfo.fr * (stimInfo.t_dur+stimInfo.ISI));
% frames_ahead = ceil(exptInfo.fr * 0.1); % ~0.1s
frames_stim = ceil(exptInfo.fr * (stimInfo.t_dur+stimInfo.ISI));
frames_ahead = 0; % ~0.1s
frames = frames_stim + frames_ahead;

%% separate spikes by stimulus
spikes_stim = zeros(num_neur, num_stim, frames);
for s = 1 : num_stim
    onset = events.eventsOn(s) - frames_ahead;
    offset = onset + frames - 1;
    spikes_stim(:, s, :) = spikes.raster(:, onset : offset);
end

%% average spikes for each word
spikes_word = zeros(num_neur, num_words, frames);
for w = 1 : num_words
    spikes_word(:, w, :) = mean(...
        spikes_stim(:, stimInfo.order==w, :), ...
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
spikes_stim = spikes_stim(:,2:end,:);
spikes_trans = zeros(num_neur, num_words * num_words, frames);
trans_all = [stimInfo.order(1:end-1)' stimInfo.order(2:end)'];
trans_uniq = unique(trans_all, 'rows');
trans_uniq = sortrows(trans_uniq,[2,1]);
num_trans = zeros(1, length(trans_uniq));
for ii = 1:length(trans_uniq)
    rows = find(trans_all(:,1)==trans_uniq(ii,1) &...
        trans_all(:,2)==trans_uniq(ii,2));
    num_trans(ii) = length(rows);
    spikes_trans(:,ii,:) = mean(spikes_stim(:,rows(20:40),:),2);
end

%% plot spikes by transition
figure
grammar_flat = reshape(stimInfo.grammar', [1 9]);
for t = 1 : length(trans_uniq)
    subplot(num_words, num_words, t)
    [~, max_spike_indices] = sort(max(spikes_trans(:,t,:), [], 3), 'descend');
    imagesc(squeeze(spikes_trans(max_spike_indices, t, :)), [0 0.5])
    xlabel('frames'); ylabel('neurons');
    title(['P(' num2str(trans_uniq(t,1)) '->' ...
        num2str(trans_uniq(t,2)) ')=' num2str(grammar_flat(t))]);
    colorbar
end
disp('neurons sorted by max spike')

%% average spikes for all neurons for each transition
spikes_trans_avg = mean(mean(spikes_trans, 3), 1);
spikes_trans_avg = reshape(spikes_trans_avg, [3 3])';
c = corr(stimInfo.grammar(:), spikes_trans_avg(:));
disp(['corr b/t grammar & mean response = ' num2str(c)])

% plot
figure
imagesc(spikes_trans_avg)
colormap('gray'); colorbar



