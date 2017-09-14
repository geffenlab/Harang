
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
uT = unique(trans_all, 'rows');
uT = sortrows(uT,[2,1]);
num_trans = zeros(1, length(uT));

for ii = 1:length(uT)
    rows = find(trans_all(:,1)==uT(ii,1) & trans_all(:,2)==uT(ii,2));
    num_trans(ii) = length(rows);
    spikes_trans(:,ii,:) = mean(spikes_stim(:,rows,:),2);
end
% for t = 1 : length(trans_unique)
%     trans = trans_unique(t,:);
%     indices_trans = find(trans_all(:,1)==trans(1) & ...
%         trans_all(:,2)==trans(2));
%     num_trans(t) = length(indices_trans);
% %     spikes_trans(:,t,:) = mean(spikes_stim(:,indices_trans(1:50),:), 2);
%     spikes_trans(:,t,:) = mean(spikes_stim(:,indices_trans,:), 2);
% end
% disp('each frame average of 50 stimuli')

%% plot spikes by transition
figure
grammar_flat = reshape(stimInfo.grammar', [1 9]);
for t = 1 : length(uT)
    subplot(num_words, num_words, t)
    [~, max_spike_indices] = sort(max(spikes_trans(:,t,:), [], 3), 'descend');
    imagesc(squeeze(spikes_trans(max_spike_indices, t, :)))
    caxis([0 0.5]); xlabel('frames'); ylabel('neurons');
    title(['P(' num2str(uT(t,1)) '->' ...
        num2str(uT(t,2)) ')=' num2str(grammar_flat(t))]);
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



