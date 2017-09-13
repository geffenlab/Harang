
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
for t = 1 : num_words * num_words
    trans = trans_unique(t, :);
    indices_trans = find(trans_all(:,1)==trans(1) & ...
        trans_all(:,2)==trans(2));
    spikes_trans(:, t, :) = mean(spikes_stim(:, indices_trans, :), 2);
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

%% separate neuron response by prior word
% A -> A, B -> A, C -> A
% A -> B, B -> B, C -> B
% A -> C, B -> C, C -> C
ca_word_stim_prev = {...
    zeros(num_neur, num_words, sum(stimInfo.order==1), 14),...
    zeros(num_neur, num_words, sum(stimInfo.order==2), 14),...
    zeros(num_neur, num_words, sum(stimInfo.order==num_words), 14)};
% ca_word_stim_prev = zeros(num_neur, num_words, num_words, 14);
for w = 1 : num_words
    for n = 1 : num_neur
        stim_order = find(stimInfo.order == w);
        for s = stim_order
            if s < 2
                continue
            end
            word_prev = stimInfo.order(s-1);
            ca_word_stim_prev{w}(n, word_prev, find(stim_order==s), :) = ...
                spikes_stim(n, s, :);
        end
    end
end

%% get averages
ca_word_stim_prev_mean = zeros(num_neur, num_words, num_words, 14);
for n = 1 : num_neur
    for w = 1 : num_words
        for w_prev = 1 : num_words
            ca_word_stim_prev_mean(n, w, w_prev, :) = ...
                mean(ca_word_stim_prev{w}(n, w_prev, :, :), num_words);
        end
    end
end

%% plot neuron response to A by prior
figure
for w = 1 : num_words
    subplot(num_words,num_words,num_words*(w-1)+1)
    ca_word_A_giv_A = squeeze(ca_word_stim_prev_mean(:, w, 1, :));
    max_resps_A_giv_A = max(ca_word_A_giv_A,[],2);
    [~, max_resps_A_giv_A_index] = sort(max_resps_A_giv_A, 'descend');
    imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index(1:50), :));
%     imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_A);
    caxis([0 .2]); colorbar;
    title(['P(word 1 to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, 1))])
    ylabel('same neurons across rows')
    
    subplot(num_words,num_words,num_words*(w-1)+2)
    ca_word_A_giv_B = squeeze(ca_word_stim_prev_mean(:, w, 2, :));
    imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index(1:50), :));
%     imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_B);
    caxis([0 .2]); colorbar
    title(['P(word 2 to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, 2))])
    
    subplot(num_words,num_words,num_words*(w-1)+num_words)
    ca_word_A_giv_C = squeeze(ca_word_stim_prev_mean(:, w, num_words, :));
    imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index(1:50), :));
%     imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_C);
    caxis([0 .2]); colorbar
    title(['P(word num_words to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, num_words))])
end

%% separate neuron response by prior word 
% AND separated by FIRST HALF and SECOND HALF of the session
% A -> A, B -> A, C -> A
% A -> B, B -> B, C -> B
% A -> C, B -> C, C -> C
wordA_stim = find(stimInfo.order == 1);
wordB_stim = find(stimInfo.order == 2);
wordC_stim = find(stimInfo.order == num_words);
ca_word_stim_prev_halves = {...
    zeros(num_neur, num_words, sum(wordA_stim<1000), 14),...
    zeros(num_neur, num_words, sum(wordB_stim<1000), 14),...
    zeros(num_neur, num_words, sum(wordC_stim<1000), 14);...
    zeros(num_neur, num_words, sum(wordA_stim>=1000), 14),...
    zeros(num_neur, num_words, sum(wordB_stim>=1000), 14),...
    zeros(num_neur, num_words, sum(wordC_stim>=1000), 14)};
for w = 1 : num_words
    for n = 1 : num_neur
        stim_order = find(stimInfo.order == w);
        for s = stim_order
            if s < 2
                continue
            end
            word_prev = stimInfo.order(s-1);
            which_half = 1;
            if s >= 1000
                which_half = 2;
            end
            ca_word_stim_prev_halves{w, which_half}...
                (n, word_prev, find(stim_order==s), :) = ...
                spikes_stim(n, s, :);
        end
    end
end

%% get averages
ca_word_stim_prev_halves_mean = zeros(num_neur, 2, num_words, num_words, 14);
for n = 1 : num_neur
    for which_half = 1 : 2
        for w = 1 : num_words
            for w_prev = 1 : num_words
                ca_word_stim_prev_halves_mean...
                    (n, which_half, w, w_prev, :) = ...
                    mean(ca_word_stim_prev_halves...
                    {w, which_half}(n, w_prev, :, :), num_words);
            end
        end
    end
end

%% plot neuron response to A by prior
% AND separated by FIRST HALF and SECOND HALF of the session
for half = 1 : 2
    figure(half)
    for w = 1 : num_words
        subplot(num_words,num_words,num_words*(w-1)+1)
        ca_word_A_giv_A = squeeze(...
            ca_word_stim_prev_halves_mean(:, half, w, 1, :));
        max_resps_A_giv_A = max(ca_word_A_giv_A,[],2);
        [~, max_resps_A_giv_A_index] = sort(max_resps_A_giv_A, 'descend');
%         imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index(1:50), :));
        imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index, :));
%         imagesc(ca_word_A_giv_A);
        caxis([0 .2]); colorbar;
        title(['P(word 1 to ' num2str(w) ') = ' ...
            num2str(stimInfo.grammar(w, 1)) ' (' num2str(half) ' half)'])
        ylabel('same neurons across rows')

        subplot(num_words,num_words,num_words*(w-1)+2)
        ca_word_A_giv_B = squeeze(...
            ca_word_stim_prev_halves_mean(:, half, w, 2, :));
%         imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index(1:50), :));
        imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index, :));
%         imagesc(ca_word_A_giv_B);
        caxis([0 .2]); colorbar
        title(['P(word 2 to ' num2str(w) ') = ' ...
            num2str(stimInfo.grammar(w, 2)) ' (' num2str(half) ' half)'])

        subplot(num_words,num_words,num_words*(w-1)+num_words)
        ca_word_A_giv_C = squeeze(...
            ca_word_stim_prev_halves_mean(:, half, w, num_words, :));
%         imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index(1:50), :));
        imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index, :));
%         imagesc(ca_word_A_giv_C);
        caxis([0 .2]); colorbar
        title(['P(word num_words to ' num2str(w) ') = ' ...
            num2str(stimInfo.grammar(w, num_words)) ' (' num2str(half) ' half)'])
    end
end
