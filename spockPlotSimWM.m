function spockPlotSimWM(my_task_id)

%% Move to correct directory
if ispc,
    base_dir = 'B:\Projects\Models of Working Memory\Recurrent-Random Networks';
elseif isunix,
    base_dir = '/jukebox/buschman/Projects/Models of Working Memory/Recurrent-Random Networks';
end
cd(base_dir);

%% Convert task id into parameters setting
sim.RndRec_f = [0:0.1:1];
sim.RndWBaseline = [0 -4];
sim.RndToRecWPositive = [0 10 30];
sim.NumInputs = [1:4];

mysim.RndWBaseline = sim.RndWBaseline(mod(my_task_id, length(sim.RndWBaseline))+1);
my_task_id = floor(my_task_id./length(sim.RndWBaseline));
mysim.RndToRecWPositive = sim.RndToRecWPositive(mod(my_task_id, length(sim.RndToRecWPositive))+1);
my_task_id = floor(my_task_id./length(sim.RndToRecWPositive));
mysim.NumInputs = sim.NumInputs(my_task_id + 1);

fprintf('Plotting recurrent working memory networks with random control units.\n');
fprintf('Parameters:\n');
fn = fieldnames(mysim);
for i = 1:length(fn),
    fprintf('%s: %s\n', fn{i}, mat2str(mysim.(fn{i})));
end

%% Determine the memory estimate over time, sparsity measurement in random network, etc
active_thresh = 3;
fprintf('\nActive neurons are those with a firing rate greater than %d Hz.\n\n', active_thresh);

prct_active_fig = figure;
mem_vect_ang_fig = figure;
mem_vect_len_fig = figure;
mem_peak_height_fig = figure;
mem_peak_ang_fig = figure;
mem_peak_width_fig = figure;
mem_peak_prom_fig = figure;

num_trials = 50;
ovr_vect_ang = NaN*ones(length(sim.RndRec_f), num_trials);
ovr_vect_len = NaN*ones(length(sim.RndRec_f), num_trials);
ovr_peak_height = NaN*ones(length(sim.RndRec_f), num_trials);
ovr_peak_ang = NaN*ones(length(sim.RndRec_f), num_trials);
ovr_peak_width = NaN*ones(length(sim.RndRec_f), num_trials);
ovr_peak_prom = NaN*ones(length(sim.RndRec_f), num_trials);

time_avg = [0.5 1];

for cur_f = sim.RndRec_f,
    
    cur_fn = sprintf('SimWM_f%03.0f_WBase%03.0f_RtRPos%03.0f_Inputs%1.0f.mat', ...
        cur_f*100, mysim.RndWBaseline, mysim.RndToRecWPositive, mysim.NumInputs);
    if ~exist(cur_fn, 'file'),
        fprintf('File %s doesn''t exist. Skipping.\n', cur_fn);
        continue;
    end
    
    %Open up matfile
    try
        data = matfile(cur_fn, 'Writable', true);
    catch err
        fprintf('Ran into error reading file %s. Skipping.\n', cur_fn);
        continue;
    end
    t = data.t;
    opts = data.opts;
    
    %Time average indexes
    time_avg_ind = ((t >= time_avg(1)) & (t <= time_avg(2)));
    
    %Estimate the percent of active cells
    prct_active_rnd = squeeze(nanmean(double(data.r_rnd) >= active_thresh, 1));
    data.Rnd_PrctActive = prct_active_rnd;
    
    figure(prct_active_fig);
    PlotWithStd(gca, t, prct_active_rnd, sprintf('f=%3.2f', cur_f), 0);
    
    %Estimate accuracy of representation
    neuron_angs = [1:size(data, 'r_rec', 1)]'./size(data, 'r_rec', 1)*2*pi;
    mem_ang = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    mem_len = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    mem_peak_ang = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    mem_peak_height = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    mem_peak_width = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    mem_peak_prom = NaN*ones(size(data, 'r_rec', 3), size(data, 'r_rec', 4));
    for i = 1:size(data, 'r_rec', 3),
        temp_r_rec = squeeze(double(data.r_rec(:, 1, i, :)));
        temp_vect = nanmean(temp_r_rec.*repmat(cos(neuron_angs) + sqrt(-1)*sin(neuron_angs), [1 size(data, 'r_rec', 4)]), 1);
        
        mem_ang(i, :) = angle(temp_vect);
        mem_len(i, :) = abs(temp_vect);
        
        for j = 1:size(temp_r_rec, 2),
            %Find peak
            [pks, locs, w, p] = findpeaks(temp_r_rec(:, j), 'SortStr','descend');
            if ~isempty(pks),
                sat_width = find(diff(temp_r_rec(locs(1):end, j) == temp_r_rec(locs(1), j)) == -1, 1, 'first');
                locs = locs(1) + (sat_width - 1)/2;
                
                mem_peak_height(i, j) = pks(1);
                mem_peak_ang(i, j) = neuron_angs(round(locs));
                mem_peak_width(i, j) = w(1);
                mem_peak_prom(i, j) = p(1);
            end
        end %trial loop
    end %loop through time
    
    %Correct for memorandum
    enc_ang = neuron_angs(opts.InputCenter(1));
    mem_ang = mod((mem_ang - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
    mem_peak_ang = mod((mem_peak_ang - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
    
    %Save to file
    data.MemVect_Ang = mem_ang;
    data.MemVect_Len = mem_len;
    data.MemPeakHeight = mem_peak_height;
    data.MemPeakAng = mem_peak_ang;
    data.MemPeakWidth = mem_peak_width;
    data.MemPeakProm = mem_peak_prom;
    
    figure(mem_vect_ang_fig);
    PlotWithStd(gca, t, abs(mem_ang), sprintf('f=%3.2f', cur_f), 0);
    ovr_vect_ang(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(abs(mem_ang(time_avg_ind, :)), 1);
    
    figure(mem_vect_len_fig);
    PlotWithStd(gca, t, mem_len, sprintf('f=%3.2f', cur_f), 0);
    ovr_vect_len(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(mem_len(time_avg_ind, :), 1);
    
    figure(mem_peak_height_fig);
    PlotWithStd(gca, t, mem_peak_height, sprintf('f=%3.2f', cur_f), 0);
    ovr_peak_height(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(mem_peak_height(time_avg_ind, :), 1);
    
    figure(mem_peak_ang_fig);
    PlotWithStd(gca, t, abs(mem_peak_ang), sprintf('f=%3.2f', cur_f), 0);
    ovr_peak_ang(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(abs(mem_peak_ang(time_avg_ind, :)), 1);
    
    figure(mem_peak_width_fig);
    PlotWithStd(gca, t, mem_peak_width, sprintf('f=%3.2f', cur_f), 0);
    ovr_peak_width(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(mem_peak_width(time_avg_ind, :), 1);
    
    figure(mem_peak_prom_fig);
    PlotWithStd(gca, t, mem_peak_prom, sprintf('f=%3.2f', cur_f), 0);
    ovr_peak_prom(find(sim.RndRec_f == cur_f, 1, 'first'), :) = nanmean(mem_peak_prom(time_avg_ind, :), 1);
    
    fprintf('Added f=%3.2f to plots.\n', cur_f);
    
    % Plot an example trial
    figh = PlotSimWM_ExampleTrial(data, 5);
    saveas(figh, sprintf('SimWM_WBase%03.0f_RtRPos%03.0f_Input%03.0f_ExampleTrial.fig', mysim.RndWBaseline, mysim.RndToRecWPositive, mysim.NumInputs));
    close(figh);
    
    clear('data');
end %f loop

figure(mem_vect_ang_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(mem_vect_len_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(mem_peak_height_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(mem_peak_width_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(mem_peak_prom_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(mem_peak_ang_fig);  title(sprintf('%d Inputs', mysim.NumInputs));
figure(prct_active_fig);  title(sprintf('%d Inputs', mysim.NumInputs));

%Save figures to file
fig_pre_str = sprintf('SimWM_WBase%03.0f_RtRPos%03.0f_Input%03.0f', mysim.RndWBaseline, mysim.RndToRecWPositive, mysim.NumInputs);
saveas(prct_active_fig, [fig_pre_str 'RndPrctActive.fig']);
saveas(mem_vect_ang_fig, [fig_pre_str 'MemVectAngle.fig']);
saveas(mem_vect_len_fig, [fig_pre_str 'MemVectLength.fig']);
saveas(mem_peak_ang_fig, [fig_pre_str 'MemPeakAngle.fig']);
saveas(mem_peak_height_fig, [fig_pre_str 'MemPeakHeight.fig']);
saveas(mem_peak_width_fig, [fig_pre_str 'MemPeakWidth.fig']);
saveas(mem_peak_prom_fig, [fig_pre_str 'MemPeakProm.fig']);

%% Do overall plots
figure;
PlotWithStd(gca, sim.RndRec_f, ovr_vect_ang, '', 1);
saveas(gcf, [fig_pre_str 'MemVectAngleByF.fig']);

figure;
PlotWithStd(gca, sim.RndRec_f, ovr_vect_len, '', 1);
saveas(gcf, [fig_pre_str 'MemVectLengthByF.fig']);

figure;
PlotWithStd(gca, sim.RndRec_f, ovr_peak_ang, '', 1);
saveas(gcf, [fig_pre_str 'PeakAngleByF.fig']);

figure;
PlotWithStd(gca, sim.RndRec_f, ovr_peak_height, '', 1);
saveas(gcf, [fig_pre_str 'PeakHeightByF.fig']);

figure;
PlotWithStd(gca, sim.RndRec_f, ovr_peak_width, '', 1);
saveas(gcf, [fig_pre_str 'PeakWidthByF.fig']);

figure;
PlotWithStd(gca, sim.RndRec_f, ovr_peak_prom, '', 1);
saveas(gcf, [fig_pre_str 'PeakPromByF.fig']);


function h = PlotWithStd(cur_axes, t, x, leg_str, inc_std),

axes(cur_axes);
h = plot(t, nanmean(x, 2), '-'); hold all;
if inc_std,
    plot(t, nanmean(x, 2) + nanstd(x, [], 2), ':', 'Color', get(h, 'Color'));
    plot(t, nanmean(x, 2) - nanstd(x, [], 2), ':', 'Color', get(h, 'Color'));
end
legend(cat(2, get(legend, 'String'), {leg_str}));
