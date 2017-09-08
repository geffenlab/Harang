
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

%% find event onset for A, B, and C
eventsAOn = events.eventsOn(stimInfo.order==1);
eventsBOn = events.eventsOn(stimInfo.order==2);
eventsCOn = events.eventsOn(stimInfo.order==3);

%% find mean Ca response
frames_per_sample = 11; % ~0.35s
frames_ahead = 3; % ~0.1s
% frames_ahead = 0; % 0s
ca_word_stim = zeros(316, 2000, frames_per_sample + frames_ahead);
for n = 1 : 316
    for s = 1 : 2000
        onset = events.eventsOn(s) - frames_ahead;
        offset = onset + frames_per_sample + frames_ahead - 1;
        ca_word_stim(n, s, :) = ...
            spikes.raster(n, onset : offset);
%             calcium.npilSubTraces(n, onset : offset);
%             spikes.deconvCaTraces(n, onset : offset);
    end
end
%% find Ca response per word
ca_word = zeros(316, 3, frames_per_sample + frames_ahead);
for n = 1 : 316
    for w = 1 : 3
        ca_word(n, w, :) = ...
            mean(ca_word_stim(n, find(stimInfo.order==w), :), 2);
    end
end
ca_word_mean = [...
    squeeze(mean(ca_word(:,1,:), 3))...
    squeeze(mean(ca_word(:,2,:), 3))...
    squeeze(mean(ca_word(:,3,:), 3))];
% ca_word_std = [...
%     squeeze(std(ca_word(:,1,:), 3))...
%     squeeze(std(ca_word(:,2,:), 3))...
%     squeeze(std(ca_word(:,3,:), 3))];

%% plot neuron Ca response to 3 words
bar(ca_word_mean(1:10,:))

%% plot neuron response across time (one stimuli at a time)
for n = 1 : 316
    plot(squeeze(ca_word(n, 1, :)))
    hold on
    plot(squeeze(ca_word(n, 2, :)))
    plot(squeeze(ca_word(n, 3, :)))
    axis([0 12 0 0.2])
    title(['neuron ' num2str(n)])
    hold off
    waitforbuttonpress
end

%% plot neuron response across time
num_neur_show = 50;
% A
subplot(1,3,1)
ca_word_A = squeeze(ca_word(:, 1, :));
max_resps_A = max(ca_word_A,[],2);
[~, max_resps_A_index] = sort(max_resps_A, 'descend'    );
imagesc(ca_word_A(max_resps_A_index(2 : num_neur_show), :))
% B
subplot(1,3,2)
ca_word_B = squeeze(ca_word(:, 2, :));
max_resps_B = max(ca_word_B,[],2);
[~, max_resps_B_index] = sort(max_resps_B, 'descend');
imagesc(ca_word_B(max_resps_A_index(2 : num_neur_show), :))
% C
subplot(1,3,3)
ca_word_C = squeeze(ca_word(:, 3, :));
max_resps_C = max(ca_word_C,[],2);
[~, max_resps_C_index] = sort(max_resps_C, 'descend');
imagesc(ca_word_C(max_resps_A_index(2 : num_neur_show), :))

%% separate neuron response by prior word
% A -> A, B -> A, C -> A
% A -> B, B -> B, C -> B
% A -> C, B -> C, C -> C
ca_word_stim_prev = {...
    zeros(316, 3, sum(stimInfo.order==1), 14),...
    zeros(316, 3, sum(stimInfo.order==2), 14),...
    zeros(316, 3, sum(stimInfo.order==3), 14)};
% ca_word_stim_prev = zeros(316, 3, 3, 14);
for w = 1 : 3
    for n = 1 : 316
        stim_order = find(stimInfo.order == w);
        for s = stim_order
            if s < 2
                continue
            end
            word_prev = stimInfo.order(s-1);
            ca_word_stim_prev{w}(n, word_prev, find(stim_order==s), :) = ...
                ca_word_stim(n, s, :);
        end
    end
end

%% get averages
ca_word_stim_prev_mean = zeros(316, 3, 3, 14);
for n = 1 : 316
    for w = 1 : 3
        for w_prev = 1 : 3
            ca_word_stim_prev_mean(n, w, w_prev, :) = ...
                mean(ca_word_stim_prev{w}(n, w_prev, :, :), 3);
        end
    end
end

%% plot neuron response to A by prior
for w = 1 : 3
    subplot(3,3,3*(w-1)+1)
    ca_word_A_giv_A = squeeze(ca_word_stim_prev_mean(:, w, 1, :));
    max_resps_A_giv_A = max(ca_word_A_giv_A,[],2);
    [~, max_resps_A_giv_A_index] = sort(max_resps_A_giv_A, 'descend');
%     imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index(1:50), :));
    imagesc(ca_word_A_giv_A(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_A);
    caxis([0 .2]); colorbar;
    title(['P(word 1 to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, 1))])
    ylabel('same neurons across rows')
    
    subplot(3,3,3*(w-1)+2)
    ca_word_A_giv_B = squeeze(ca_word_stim_prev_mean(:, w, 2, :));
%     imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index(1:50), :));
    imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_B);
    caxis([0 .2]); colorbar
    title(['P(word 2 to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, 2))])
    
    subplot(3,3,3*(w-1)+3)
    ca_word_A_giv_C = squeeze(ca_word_stim_prev_mean(:, w, 3, :));
%     imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index(1:50), :));
    imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index, :));
%     imagesc(ca_word_A_giv_C);
    caxis([0 .2]); colorbar
    title(['P(word 3 to ' num2str(w) ') = ' num2str(stimInfo.grammar(w, 3))])
end

%% separate neuron response by prior word 
% AND separated by FIRST HALF and SECOND HALF of the session
% A -> A, B -> A, C -> A
% A -> B, B -> B, C -> B
% A -> C, B -> C, C -> C
wordA_stim = find(stimInfo.order == 1);
wordB_stim = find(stimInfo.order == 2);
wordC_stim = find(stimInfo.order == 3);
ca_word_stim_prev_halves = {...
    zeros(316, 3, sum(wordA_stim<1000), 14),...
    zeros(316, 3, sum(wordB_stim<1000), 14),...
    zeros(316, 3, sum(wordC_stim<1000), 14);...
    zeros(316, 3, sum(wordA_stim>=1000), 14),...
    zeros(316, 3, sum(wordB_stim>=1000), 14),...
    zeros(316, 3, sum(wordC_stim>=1000), 14)};
for w = 1 : 3
    for n = 1 : 316
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
                ca_word_stim(n, s, :);
        end
    end
end

%% get averages
ca_word_stim_prev_halves_mean = zeros(316, 2, 3, 3, 14);
for n = 1 : 316
    for which_half = 1 : 2
        for w = 1 : 3
            for w_prev = 1 : 3
                ca_word_stim_prev_halves_mean...
                    (n, which_half, w, w_prev, :) = ...
                    mean(ca_word_stim_prev_halves...
                    {w, which_half}(n, w_prev, :, :), 3);
            end
        end
    end
end

%% plot neuron response to A by prior
% AND separated by FIRST HALF and SECOND HALF of the session
for half = 1 : 2
    figure(half)
    for w = 1 : 3
        subplot(3,3,3*(w-1)+1)
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

        subplot(3,3,3*(w-1)+2)
        ca_word_A_giv_B = squeeze(...
            ca_word_stim_prev_halves_mean(:, half, w, 2, :));
%         imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index(1:50), :));
        imagesc(ca_word_A_giv_B(max_resps_A_giv_A_index, :));
%         imagesc(ca_word_A_giv_B);
        caxis([0 .2]); colorbar
        title(['P(word 2 to ' num2str(w) ') = ' ...
            num2str(stimInfo.grammar(w, 1)) ' (' num2str(half) ' half)'])

        subplot(3,3,3*(w-1)+3)
        ca_word_A_giv_C = squeeze(...
            ca_word_stim_prev_halves_mean(:, half, w, 3, :));
%         imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index(1:50), :));
        imagesc(ca_word_A_giv_C(max_resps_A_giv_A_index, :));
%         imagesc(ca_word_A_giv_C);
        caxis([0 .2]); colorbar
        title(['P(word 3 to ' num2str(w) ') = ' ...
            num2str(stimInfo.grammar(w, 1)) ' (' num2str(half) ' half)'])
    end
end
