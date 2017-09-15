load('F:\HarangData\K070_20170903_artificialGrammar_02.mat')

stim = stimInfo.order; % order of word presentation
uStim = unique(stim); % how many words
e = events.eventsOn; % event times (i.e. stim Onsets)
fr = exptInfo.fr; % frame rate
stimDur = ceil((stimInfo.t_dur+stimInfo.ISI)*fr);
stimPrec = ceil(.1*fr);

% split up the spikes into stim presentations
raster = rasterize(spikes.raster, stimDur, events.eventsOn);



%%
g = stimInfo.grammar;
raster = raster(:,:,2:end); % remove fisrt trial as nothng preceded it
newOrder = [stim(1:end-1)' stim(2:end)']; % order of transitions
uT = unique(newOrder,'rows'); % unique transitions
uT = sortrows(uT,[2,1]);
trials = size(newOrder,1);
tChunks = 2;
bs = 100;
bsSamp = 1000;
chunkSize = floor(trials/tChunks);
trast = zeros(size(raster,1),size(raster,2),length(uT),tChunks);
correl = zeros(1,tChunks);
n = zeros(tChunks, length(uT));
for c = 1:tChunks
    on = (c-1)*chunkSize + 2;
    off = on + chunkSize - 1;
    for t = 1:length(uT)
        rows = find(newOrder(on:off,1)==uT(t,1) & ...
            newOrder(on:off,2)==uT(t,2));
        n(c,t) = length(rows);
%         boot = zeros(bs,size(raster,1),size(raster,2));
%         for b = 1 : bs
%             if length(rows) >= bsSamp
%                 randIndex = randperm(length(rows),bsSamp);
%             else
%                 randIndex = randperm(length(rows));
%             end
%             boot(b,:,:) = mean(raster(:,:,rows(randIndex)), 3);
%         end
        trast(:,:,t,c) = mean(boot,1);
    end
    mtr = squeeze(mean(mean(trast(:,:,:,c),1),2));
    correl(c) = corr(g(:),mtr);
    disp(c)
end
disp(correl)

%% Plot each transition raster order by most responsive neurons
grammar_flat = reshape(stimInfo.grammar', [1 9]);
for c = 1:tChunks
    figure
    for ii = 1:size(trast(:,:,:,c),3)
        subplot(3,3,ii)
        [~,index] = sort(max(trast(:,:,ii,c),[],2),'descend');
        a = trast(index,:,ii);
        imagesc(a(1:50,:), [0 0.5])
        title([ num2str(uT(ii, 1)) '->' num2str(uT(ii,2))])
        xlabel('frames'); ylabel('neurons')
        title(['P(' num2str(uT(ii,1)) '->' ...
            num2str(uT(ii,2)) ')=' num2str(grammar_flat(ii))]);
        colorbar
    end
end

%% Work out mean response to each transition across neurons
mtr = squeeze(mean(mean(trast(:,:,:,1),1),2));
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

