close all;

% Define parameters
stim_period = [0.1 0.3];
mem_period = [0.3 1];

%% 
stim_ind = ((t >= stim_period(1)) & (t <= stim_period(2)));
mem_ind = ((t >= mem_period(1)) & (t <= mem_period(2)));


%% Plot example neuron's selectivity

figure;
subplot(1,2,1);
plot(squeeze(nanmean(r_rec(1:round(size(r_rec, 1)/8):end, 1, stim_ind, 1:8), 3))');
subplot(1,2,2);
plot(squeeze(nanmean(r_rnd(1:round(size(r_rnd, 1)/8):end, stim_ind, 1:8), 2))');
title('Example neuron responses to individual stimuli into Sensory Network #1');


%% Plot the responses of neurons across stimuli, single stimulus alone into sensory network #1
figure;
subplot(2,2,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, stim_ind, 1:8), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network, Sorted by Preferred Stimulus');

subplot(2,2,3);
temp_resp = squeeze(nanmean(r_rec(:, 2, stim_ind, 1:8), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred Stimulus');

subplot(2,2,[2 4]);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, 1:8), 2));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stimulus');


%% Plot the responses of neurons across stimuli, single stimulus alone into sensory network #2
figure;
subplot(2,2,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, stim_ind, 9:16), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network, Sorted by Preferred Stimulus');

subplot(2,2,3);
temp_resp = squeeze(nanmean(r_rec(:, 2, stim_ind, 9:16), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred Stimulus');

subplot(2,2,[2 4]);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, 9:16), 2));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stimulus');

%% Plot the responses of neurons across stimuli, two stimuli presentations

figure;
subplot(2,3,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, stim_ind, 17:80), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #1, Sorted by Preferred 2-Item Stimulus');
set(gca, 'CLim', [0 80]);

subplot(2,3,4);
temp_resp = squeeze(nanmean(r_rec(:, 2, stim_ind, 17:80), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred 2-Item Stimulus');
set(gca, 'CLim', [0 80]);

subplot(2,3,[2 5]);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, stim_ind, 17:80), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred 2-Item Stimulus');
set(gca, 'CLim', [0 80]);

subplot(2,3,3);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, stim_ind, 1:8), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stim #1 Stimulus');
set(gca, 'CLim', [0 80]);

subplot(2,3,6);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, stim_ind, 9:16), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stim #2 Stimulus');
set(gca, 'CLim', [0 80]);

%% Plot the time course of stimuli to their 'preferred' stimulus

rec_plot_ind = [1:round(size(r_rec, 1)/5):size(r_rec, 1)];
rnd_plot_ind = [1:round(size(r_rnd, 1)/5):size(r_rnd, 1)];

figure;

subplot(1,2,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, stim_ind, :), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
for i = 1:length(rec_plot_ind),
    plot(t, squeeze(r_rec(rec_plot_ind(i), 1, :, max_ind(rec_plot_ind(i))))); hold all;
end
ovr_sum_resp = zeros(length(t), 1);
for i = 1:size(r_rec, 1),
    ovr_sum_resp = ovr_sum_resp + squeeze(r_rec(i, 1, :, max_ind(i)));
end
ovr_sum_resp = ovr_sum_resp./size(r_rec, 1);
plot(t, ovr_sum_resp, 'k-', 'LineWidth', 2);
title('Example (and Average) Sensory Responses to Preferred Stimulus');

subplot(1,2,2);
temp_resp = squeeze(nanmean(r_rnd(:, stim_ind, :), 2));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
for i = 1:length(rnd_plot_ind),
    plot(t, squeeze(r_rnd(rnd_plot_ind(i), :, max_ind(rnd_plot_ind(i))))); hold all;
end
ovr_sum_resp = zeros(1, length(t));
for i = 1:size(r_rnd, 1),
    ovr_sum_resp = ovr_sum_resp + squeeze(r_rnd(i, :, max_ind(i)));
end
ovr_sum_resp = ovr_sum_resp./size(r_rnd, 1);
plot(t, ovr_sum_resp, 'k-', 'LineWidth', 2);
title('Example (and Average) Random Responses to Preferred Stimulus');

%% Plot example neuron's selectivity for memory delay

figure;
subplot(1,2,1);
plot(squeeze(nanmean(r_rec(1:round(size(r_rec, 1)/8):end, 1, mem_ind, 1:8), 3))');
subplot(1,2,2);
plot(squeeze(nanmean(r_rnd(1:round(size(r_rnd, 1)/8):end, mem_ind, 1:8), 2))');
title('Example neuron responses to individual memories into Sensory Network #1');


%% Plot the responses of neurons across stimuli, single stimulus alone into sensory network #1
figure;
subplot(2,2,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, mem_ind, 1:8), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network, Sorted by Preferred Memory');

subplot(2,2,3);
temp_resp = squeeze(nanmean(r_rec(:, 2, mem_ind, 1:8), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred Memory');

subplot(2,2,[2 4]);
temp_resp = squeeze(nanmean(r_rnd(:, mem_ind, 1:8), 2));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Memory');


%% Plot the responses of neurons across stimuli, single stimulus alone into sensory network #2
figure;
subplot(2,2,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, mem_ind, 9:16), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network, Sorted by Preferred Memory');

subplot(2,2,3);
temp_resp = squeeze(nanmean(r_rec(:, 2, mem_ind, 9:16), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred Memory');

subplot(2,2,[2 4]);
temp_resp = squeeze(nanmean(r_rnd(:, mem_ind, 9:16), 2));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Memory');

%% Plot the responses of neurons across stimuli, two stimuli presentations

figure;
subplot(2,3,1);
temp_resp = squeeze(nanmean(r_rec(:, 1, mem_ind, 17:80), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #1, Sorted by Preferred 2-Item Memory');
set(gca, 'CLim', [0 80]); colorbar;

subplot(2,3,4);
temp_resp = squeeze(nanmean(r_rec(:, 2, mem_ind, 17:80), 3));
[~, max_ind] = max(temp_resp, [], 2);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Stimulus Network #2, Sorted by Preferred 2-Item Memory');
set(gca, 'CLim', [0 80]); colorbar;

subplot(2,3,[2 5]);
temp_resp = squeeze(nanmean(r_rnd(:, mem_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, mem_ind, 17:80), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred 2-Item Memory');
set(gca, 'CLim', [0 40]); colorbar;

subplot(2,3,3);
temp_resp = squeeze(nanmean(r_rnd(:, mem_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, mem_ind, 1:8), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stim #1 Memory');
set(gca, 'CLim', [0 40]); colorbar;

subplot(2,3,6);
temp_resp = squeeze(nanmean(r_rnd(:, mem_ind, 17:80), 2));
[~, max_ind] = max(nanmean(r_rnd(:, mem_ind, 9:16), 2), [], 3);
[~, sort_ind] = sort(max_ind);
imagesc(temp_resp(sort_ind, :));
title('Random Network, Sorted by Preferred Stim #2 Memory');
set(gca, 'CLim', [0 40]); colorbar;

%% Plot the responses of neurons comparing the response to two single stimuli and their combination

[i1,i2] = meshgrid([1:8], [1:8]);
rnd_stim = squeeze(nanmean(r_rnd(:, stim_ind, :), 2));
rnd_mem = squeeze(nanmean(r_rnd(:, mem_ind, :), 2));
for i = 1:size(rnd_stim, 1),
    [b_stim(:, i), ~, stats_stim(i)] = glmfit(cat(2, i1(:), i2(:)), rnd_stim(i, 17:80)');
    [b_mem(:, i), ~, stats_mem(i)] = glmfit(cat(2, i1(:), i2(:)), rnd_mem(i, 17:80)');
end

figure;
imagesc(rnd_stim(:, 17:80) - (rnd_stim(:, i1(:)) + rnd_stim(:, i2(:) + 8)));

resid = rnd_stim(:, 17:80) - (rnd_stim(:, i1(:)) + rnd_stim(:, i2(:) + 8));
resid = sum(abs(resid), 2);

num_rnd_perm = 1000;
rnd_resid = NaN*ones(size(resid, 1), num_rnd_perm);
rnd_rnd_stim = rnd_stim(:, 17:80);
for i = 1:num_rnd_perm,
    %Randomly permute the associations
    rnd_rnd_stim = rnd_rnd_stim(:, randperm(size(rnd_rnd_stim, 2)));
    rnd_resid(:, i) = sum(abs(rnd_rnd_stim(:, :) - (rnd_stim(:, i1(:)) + rnd_stim(:, i2(:) + 8))), 2);
end
figure;
plot(resid, 'r-'); hold all;
plot(nanmean(rnd_resid, 2), 'k--');
plot(nanmean(rnd_resid, 2) + nanstd(rnd_resid, [], 2), 'k:');
plot(nanmean(rnd_resid, 2) - nanstd(rnd_resid, [], 2), 'k:');

hist_fig = figure;
[n,x] = hist((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2), 100);
plot(x, n); hold all;
v = axis;
plot(nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2))*[1 1], v(3:4), 'r-');
text(nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2)), v(4), sprintf('%4.3f', nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2))), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');


figure;
imagesc(rnd_mem(:, 17:80) - (rnd_mem(:, i1(:)) + rnd_mem(:, i2(:) + 8)));

resid = rnd_mem(:, 17:80) - (rnd_mem(:, i1(:)) + rnd_mem(:, i2(:) + 8));
resid = sum(abs(resid), 2);

num_rnd_perm = 1000;
rnd_resid = NaN*ones(size(resid, 1), num_rnd_perm);
rnd_rnd_mem = rnd_mem(:, 17:80);
for i = 1:num_rnd_perm,
    %Randomly permute the associations
    rnd_rnd_mem = rnd_rnd_mem(:, randperm(size(rnd_rnd_mem, 2)));
    rnd_resid(:, i) = sum(abs(rnd_rnd_mem(:, :) - (rnd_mem(:, i1(:)) + rnd_mem(:, i2(:) + 8))), 2);
end

figure(hist_fig);
[n,x] = hist((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2), 100);
plot(x, n, '--');
v = axis;
plot(nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2))*[1 1], v(3:4), 'r--');
text(nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2)), v(4), sprintf('%4.3f', nanmean((resid - nanmean(rnd_resid, 2))./nanstd(rnd_resid, [], 2))), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');