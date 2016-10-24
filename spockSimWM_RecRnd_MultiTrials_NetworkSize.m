function spockSimWM_RecRnd_MultiTrials_NetworkSize(my_task_id, num_trials, num_fullsave_trials, overwrite)

%% Initialize random number generator
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

%% Move to correct directory
if ispc,
    base_dir = 'B:\Projects\Models of Working Memory\Recurrent-Random Networks';
elseif isunix,
    base_dir = '/jukebox/buschman/Projects/Models of Working Memory/Recurrent-Random Networks';
end
cd(base_dir);

%% Convert task id into parameters setting
sim.NetworkSize = round(2048*[1 0.75 0.5 0.25 0.1]);
sim.RndRec_f = [0 0.010 0.025 0.050 0.075 0.100 0.150 0.200 0.250 0.300 0.400 0.500 0.750 1.000];
sim.RecToRndW_TargetFR = [0.5];
sim.MaxAvgRndFR = [20];
sim.RecWPositiveWidth = [1];
sim.NumInputs = [1:4];

length(sim.NetworkSize)*length(sim.RndRec_f)*length(sim.RecToRndW_TargetFR)*length(sim.MaxAvgRndFR)*length(sim.RecWPositiveWidth)*length(sim.NumInputs)

mysim.NetworkSize = sim.NetworkSize(mod(my_task_id, length(sim.NetworkSize))+1);
my_task_id = floor(my_task_id./length(sim.NetworkSize));
mysim.RndRec_f = sim.RndRec_f(mod(my_task_id, length(sim.RndRec_f))+1);
my_task_id = floor(my_task_id./length(sim.RndRec_f));
mysim.RecToRndW_TargetFR = sim.RecToRndW_TargetFR(mod(my_task_id, length(sim.RecToRndW_TargetFR))+1);
my_task_id = floor(my_task_id./length(sim.RecToRndW_TargetFR));
mysim.MaxAvgRndFR = sim.MaxAvgRndFR(mod(my_task_id, length(sim.MaxAvgRndFR))+1);
my_task_id = floor(my_task_id./length(sim.MaxAvgRndFR));
mysim.RecWPositiveWidth = sim.RecWPositiveWidth(mod(my_task_id, length(sim.RecWPositiveWidth))+1);
my_task_id = floor(my_task_id./length(sim.RecWPositiveWidth));
mysim.NumInputs = sim.NumInputs(my_task_id + 1);

my_save_fn = sprintf('SimWM_RndN%04.0f_f%04.0f_TargetFR%03.0f_MaxFR%03.0f_RecPosWidth%03.0f_Inputs%1.0f.mat', ...
    mysim.NetworkSize, mysim.RndRec_f*1000, mysim.RecToRndW_TargetFR*100, mysim.MaxAvgRndFR, mysim.RecWPositiveWidth*100, mysim.NumInputs);
	
if exist(my_save_fn, 'file') & ~overwrite,
    fprintf('File %s already exists. Skipping.\n', my_save_fn);
    return;
end

fprintf('Simulating a recurrent working memory network with random control units.\n');
fprintf('Parameters:\n');
fn = fieldnames(mysim);
for i = 1:length(fn),
	fprintf('%s: %s\n', fn{i}, mat2str(mysim.(fn{i})));
end

%% Define my von mises function
vm = @(x, mu, k, amp, base) amp.*exp(k*cos(x - mu))./(2*pi*besseli(0, k)) + base;
fminsearch_opts = optimset('Display', 'off', 'MaxFunEvals', 10^4, 'MaxIter', 10^4);

%% Run my simulation

%Simulate a bunch of full networks, saving as we go
cur_inp_time = repmat([0.1 0.2]', [1 4]); %in seconds
cur_inp_time(:, (mysim.NumInputs+1:end)) = NaN; %only drive the first cur_inp pools

fprintf('\nStarting simulation...\n'); start_time = now;
for cur_trial = 1:num_trials,
    [r_rec, r_rnd, t, opts] = SimWM_MultiRecRnd_new('N_rnd', mysim.NetworkSize, 'RndRec_f', mysim.RndRec_f, 'RecToRndW_TargetFR', mysim.RecToRndW_TargetFR, ...
        'MaxAvgRecFR', mysim.MaxAvgRndFR, 'MaxAvgRndFR', mysim.MaxAvgRndFR, 'RecWPositiveWidth', mysim.RecWPositiveWidth, ...
        'InputTime', cur_inp_time, 'SimTime', 1);
    if cur_trial == 1,
        fprintf('Creating MAT file...\n');
        rec_rnd_outfile = matfile(my_save_fn, 'Writable', true);
        rec_rnd_outfile.r_rec = NaN*ones(size(r_rec, 1), size(r_rec, 2), size(r_rec, 3), num_fullsave_trials);
        rec_rnd_outfile.r_rec(:, :, :, cur_trial) = r_rec;
        rec_rnd_outfile.r_rnd = NaN*ones(size(r_rnd, 1), size(r_rnd, 2), num_fullsave_trials);
        rec_rnd_outfile.r_rnd(:, :, cur_trial) = r_rnd;
        
        rec_rnd_outfile.r_rnd_avg = NaN*ones(length(t), num_trials);
        rec_rnd_outfile.r_rnd_avgFRoverN = NaN*ones(size(r_rnd, 2), num_trials);
        rec_rnd_outfile.r_rnd_avgFRoverT = NaN*ones(size(r_rnd, 1), num_trials);
    
        rec_rnd_outfile.r_rnd_prct_active = NaN*ones(length(t), num_trials);
        
        rec_rnd_outfile.MemVect_Ang = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVect_Len = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVect_AngRel = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        
        rec_rnd_outfile.MemPeakHeight = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakAngle = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakWidth = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakProm = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakAngleRel = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        
        rec_rnd_outfile.MemVonMises_Ang = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVonMises_AngRel = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVonMises_K = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVonMises_Height = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVonMises_Base = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVonMises_CircVar = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        
        rec_rnd_outfile.t = t;
        rec_rnd_outfile.opts = repmat(opts, [1 num_trials]);
    elseif (cur_trial <= num_fullsave_trials),
		rec_rnd_outfile.r_rec(:, :, :, cur_trial) = r_rec;
		rec_rnd_outfile.r_rnd(:, :, cur_trial) = r_rnd;
        rec_rnd_outfile.opts(1, cur_trial) = opts;
    end
	fprintf('Simulated trial %d. (%s)\n', cur_trial, datestr(now-start_time, 'HH:MM:SS'));
    
    %Do some simple calculations on the memories
    neuron_angs = [1:opts.N_rec]'./opts.N_rec*2*pi;
    mem_ang = NaN*ones(length(t), opts.N_rec_pools);
    mem_len = NaN*ones(length(t), opts.N_rec_pools);
    mem_peak_ang = NaN*ones(length(t), opts.N_rec_pools);
    mem_peak_height = NaN*ones(length(t), opts.N_rec_pools);
    mem_peak_width = NaN*ones(length(t), opts.N_rec_pools);
    mem_peak_prom = NaN*ones(length(t), opts.N_rec_pools);
    mem_ang_rel = NaN*ones(length(t), opts.N_rec_pools);
    mem_peak_ang_rel = NaN*ones(length(t), opts.N_rec_pools);
    mem_vm_fit = NaN*ones(5, length(t), opts.N_rec_pools);
    mem_vm_fit_ang_rel = NaN*ones(length(t), opts.N_rec_pools);
    for cur_pool = 1:opts.N_rec_pools,
        temp_vect = nanmean(squeeze(r_rec(:, cur_pool, :)).*repmat(cos(neuron_angs) + sqrt(-1)*sin(neuron_angs), [1 length(t)]), 1);
        mem_ang(:, cur_pool) = angle(temp_vect);
        mem_len(:, cur_pool) = abs(temp_vect);
        
        for cur_t = 1:length(t),
            %Find peak
            [pks, locs, w, p] = findpeaks(double(r_rec(:, cur_pool, cur_t)), 'SortStr','descend');
            if ~isempty(pks),
                sat_width = find(diff(r_rec(locs(1):end, cur_pool, cur_t) == r_rec(locs(1), cur_pool, cur_t)) == -1, 1, 'first');
                locs = locs(1) + (sat_width - 1)/2;
                
                mem_peak_height(cur_t, cur_pool) = pks(1);
                mem_peak_ang(cur_t, cur_pool) = neuron_angs(round(locs));
                mem_peak_width(cur_t, cur_pool) = w(1);
                mem_peak_prom(cur_t, cur_pool) = p(1);
            end
            
            %Fit von mises, use fminsearch with starting values based on MLE fits
            temp_data = squeeze(double(r_rec(:, cur_pool, cur_t)));
            %mem_vm_fit(1:4, cur_t, cur_pool) = fmincon(@(x) sum((vm(neuron_angs, x(1), x(2), x(3), x(4)) - t(:)').^2), [pi 1 70 0], [], [], [], [], [0 0 0 0], [2*pi Inf Inf Inf], [], fmin_opts);
            mem_vm_fit(1:4, cur_t, cur_pool) = fminsearch(@(x) sum((vm(neuron_angs, x(1), x(2), x(3), x(4)) - temp_data(:)).^2), ...
                [mem_ang(cur_t, cur_pool) mem_len(cur_t, cur_pool)^2 max(temp_data) min(temp_data)], fminsearch_opts);            
            
        end %time loop
        
        %If a memorandum was stored, calculate relative angle
        if cur_pool <= mysim.NumInputs,
            enc_ang = neuron_angs(opts.InputCenter(cur_pool));
            mem_ang_rel(:, cur_pool) = mod((mem_ang(:, cur_pool) - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
            mem_peak_ang_rel(:, cur_pool) = mod((mem_peak_ang(:, cur_pool) - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
            mem_vm_fit_ang_rel(:, cur_pool) = mod((squeeze(mem_vm_fit(1, :, cur_pool)) - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
        end
            
    end %pool loop
    
    %Create circular variance statistic
    mem_vm_fit(5, :, :) = 1 - besseli(1, mem_vm_fit(2, :, :))./besseli(0, mem_vm_fit(2, :, :));
    
    %Save to file    
    rec_rnd_outfile.r_rnd_avgFRoverN(:, cur_trial) = nanmean(r_rnd, 1)';
    rec_rnd_outfile.r_rnd_avgFRoverT(:, cur_trial) = nanmean(r_rnd, 2);
    pre_stim_rnd = r_rnd(:, t < 0.1);
    rec_rnd_outfile.r_rnd_prct_active(:, cur_trial) = nanmean(r_rnd >= (nanmean(pre_stim_rnd(:)) + 3*nanstd(pre_stim_rnd(:))), 1)';
    rec_rnd_outfile.MemVect_Ang(:, :, cur_trial) = mem_ang;
    rec_rnd_outfile.MemVect_Len(:, :, cur_trial) = mem_len;
    rec_rnd_outfile.MemPeakHeight(:, :, cur_trial) = mem_peak_height;
    rec_rnd_outfile.MemPeakAngle(:, :, cur_trial) = mem_peak_ang;
    rec_rnd_outfile.MemPeakWidth(:, :, cur_trial) = mem_peak_width;
    rec_rnd_outfile.MemPeakProm(:, :, cur_trial) = mem_peak_prom;
    rec_rnd_outfile.MemVect_AngRel(:, :, cur_trial) = mem_ang_rel;
    rec_rnd_outfile.MemPeakAngleRel(:, :, cur_trial) = mem_peak_ang_rel;
    rec_rnd_outfile.MemVonMises_Ang(:, :, cur_trial) = squeeze(mem_vm_fit(1, :, :));
    rec_rnd_outfile.MemVonMises_AngRel(:, :, cur_trial) = mem_vm_fit_ang_rel;
    rec_rnd_outfile.MemVonMises_K(:, :, cur_trial) = squeeze(mem_vm_fit(2, :, :));
    rec_rnd_outfile.MemVonMises_Height(:, :, cur_trial) = squeeze(mem_vm_fit(3, :, :));
    rec_rnd_outfile.MemVonMises_Base(:, :, cur_trial) = squeeze(mem_vm_fit(4, :, :));
    rec_rnd_outfile.MemVonMises_CircVar(:, :, cur_trial) = squeeze(mem_vm_fit(5, :, :));
    fprintf('Added statistics to file for trial %d. (%s)\n', cur_trial, datestr(now-start_time, 'HH:MM:SS'));

end %trial loop
    
%% Close the file
clear('rec_rnd_outfile');


