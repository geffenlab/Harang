
raster_m = squeeze(mean(raster, 2));
max_spike = max(max(spikes.raster));

%%
neurons = randperm(size(spikes.raster,1), 100);
plot(raster_m(neurons,:)', '-.')
hold on
plot(mean(raster_m(neurons,:), 1), 'LineWidth', 3)
hold off

%%
for n = 1 : size(raster_m,1)
    plot(spikes.raster(n,:))
    axis([0 22000 0 max_spike])
    waitforbuttonpress
end

%%
chunks = 100;
chunk_size = floor(size(spikes.raster,2)/chunks);
spike_mean_chunk = zeros(size(spikes.raster,1),chunks);
for t = 1 : chunks
    on = (t-1)*chunk_size + 1;
    off = on + chunk_size - 1;
    spike_mean_chunk(:,t) = mean(spikes.raster(:,on:off),2);
end

neurons = randperm(size(spikes.raster,1), 100);
plot(spike_mean_chunk(neurons,:)', '-.')
hold on
plot(mean(spike_mean_chunk(neurons,:), 1), 'LineWidth', 3)
axis([0 100 0 0.5]); xlabel('t'); ylabel('spikes');
plotprefs
hold off

