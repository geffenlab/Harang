load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

stim = stimInfo.order; % order of word presentation
uStim = unique(stim); % how many words
e = events.eventsOn; % event times (i.e. stim Onsets)
spikes = spikes.raster; % spike times and size
fr = exptInfo.fr; % frame rate
stimDur = floor((stimInfo.t_dur+stimInfo.ISI)*fr);

% split up the spikes into stim presentations
raster  = zeros(size(spikes,1),stimDur+1,length(e));
for ii = 1:length(e)
    raster(:,:,ii) = spikes(:,e(ii):e(ii)+stimDur);
end

%%
newOrder = [stim(1:end-1)' stim(2:end)']; % order of transitions
% raster = raster(:,:,2:end); % remove fisrt trial as nothng preceded it
uT = unique(newOrder,'rows'); % unique transitions
uT = sortrows(uT,[2,1]);

trast = zeros(size(raster,1),size(raster,2),length(uT)); % make raster of each transition
n = zeros(1,length(uT));
for ii = 1:length(uT)
    rows = find(newOrder(:,1)==uT(ii,1) & newOrder(:,2)==uT(ii,2));
    n(ii) = length(rows);
    trast(:,:,ii) = mean(raster(:,:,rows),3);
end
    
%% Plot each transition raster order by most responsive neurons
figure
grammar_flat = reshape(stimInfo.grammar', [1 9]);
for ii = 1:size(trast,3)
    subplot(3,3,ii)
    [~,index] = sort(max(trast(:,:,ii),[],2),'descend');
    a = trast(index,:,ii);
    imagesc(a(1:50,:), [0 0.5])
    title([ num2str(uT(ii, 1)) '->' num2str(uT(ii,2))])
    xlabel('frames'); ylabel('neurons')
    title(['P(' num2str(uT(ii,1)) '->' ...
        num2str(uT(ii,2)) ')=' num2str(grammar_flat(ii))]);
    colorbar
    plotprefs
end

%% Work out mean response to each transition across neurons
mtr = squeeze(mean(mean(trast,1),2));
mtr = reshape(mtr,[3,3])';
figure
% mtr = mtr-min(mtr(:));
% mtr = reshape(mtr/max(mtr),[3,3])';
g = stimInfo.grammar;
% g = g-min(g(:));
% g = g/max(g(:));
c = corr(g(:),mtr(:));
disp(['corr b/t grammar & mean response = ' num2str(c)])
% imagesc(c)
imagesc(mtr)
colorbar
colormap gray

