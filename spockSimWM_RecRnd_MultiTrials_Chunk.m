function spockSimWM_RecRnd_MultiTrials_Chunk(my_task_id, num_trials, overwrite)

%% Move to correct directory
if ispc,
    base_dir = 'B:\Projects\Models of Working Memory\Recurrent-Random Networks';
elseif isunix,
    base_dir = '/jukebox/buschman/Projects/Models of Working Memory/Recurrent-Random Networks';
end
cd(base_dir);

%% Convert task id into parameters setting
sim.RndRec_f = [0.15 0.2 0.25];
sim.NumInputs = [1:4];
sim.ChunkPercentRemap = [0 0.05 0.10 0.25];
sim.MemoryDistance = [0 pi/16 pi/8 pi/4 pi/2 pi];

length(sim.RndRec_f)*length(sim.ChunkPercentRemap)*length(sim.MemoryDistance)*length(sim.NumInputs)

mysim.RndRec_f = sim.RndRec_f(mod(my_task_id, length(sim.RndRec_f))+1);
my_task_id = floor(my_task_id./length(sim.RndRec_f));
mysim.NumInputs = sim.NumInputs(mod(my_task_id, length(sim.NumInputs))+1);
my_task_id = floor(my_task_id./length(sim.NumInputs));
mysim.ChunkPercentRemap = sim.ChunkPercentRemap(mod(my_task_id, length(sim.ChunkPercentRemap))+1);
my_task_id = floor(my_task_id./length(sim.ChunkPercentRemap));
mysim.MemoryDistance = sim.MemoryDistance(my_task_id + 1);

my_save_fn = sprintf('SimWM_f%03.0f_PrctRemap%03.0f_NearbyMems%03.0f_Inputs%1.0f.mat', ...
    mysim.RndRec_f*100, mysim.ChunkPercentRemap*100, mysim.MemoryDistance*100, mysim.NumInputs);
	
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

%% Run my simulation

%Simulate a bunch of full networks, saving as we go
cur_inp_time = repmat([0.05 0.15]', [1 4]); %in seconds
cur_inp_time(:, (mysim.NumInputs+1:end)) = NaN; %only drive the first cur_inp pools

%Create chunks, and the inputs
curChunkCenters = [3*pi/2 pi/2 NaN NaN]; 
curInputCenter = [mod(curChunkCenters(1) + mysim.MemoryDistance, 2*pi) ...
    mod(curChunkCenters(2) + mysim.MemoryDistance, 2*pi) rand*2*pi rand*2*pi];

fprintf('\nStarting simulation...\n'); start_time = now;
for cur_trial = 1:num_trials,
    [r_rec, r_rnd, t, opts] = SimWM_MultiRecRnd_chunked('RndRec_f', mysim.RndRec_f, 'InputTime', cur_inp_time, ...
        'ChunkCenters', curChunkCenters, 'InputCenter', curInputCenter, 'ChunkPercentRemap', mysim.ChunkPercentRemap);
    if cur_trial == 1,
        fprintf('Creating MAT file...\n');
        rec_rnd_outfile = matfile(my_save_fn, 'Writable', true);
        rec_rnd_outfile.r_rec = NaN*ones(size(r_rec, 1), size(r_rec, 2), size(r_rec, 3), num_trials);
        rec_rnd_outfile.r_rec(:, :, :, cur_trial) = r_rec;
        rec_rnd_outfile.r_rnd = NaN*ones(size(r_rnd, 1), size(r_rnd, 2), num_trials);
        rec_rnd_outfile.r_rnd(:, :, cur_trial) = r_rnd;
        
        rec_rnd_outfile.MemVect_Ang = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVect_Len = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakHeight = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakAngle = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakWidth = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakProm = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemVect_AngRel = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        rec_rnd_outfile.MemPeakAngleRel = NaN*ones(length(t), opts.N_rec_pools, num_trials);
        
        rec_rnd_outfile.t = t;
        rec_rnd_outfile.opts = repmat(opts, [1 num_trials]);
    else
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
        end %time loop
        
        %If a memorandum was stored, calculate relative angle
        if cur_pool <= mysim.NumInputs,
            enc_ang = neuron_angs(opts.InputCenter(cur_pool));
            mem_ang_rel(:, cur_pool) = mod((mem_ang(:, cur_pool) - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
            mem_peak_ang_rel(:, cur_pool) = mod((mem_peak_ang(:, cur_pool) - enc_ang) + pi, 2*pi) - pi; %should now be -pi to pi
        end
            
    end %pool loop
    
    %Save to file    
    rec_rnd_outfile.MemVect_Ang(:, :, cur_trial) = mem_ang;
    rec_rnd_outfile.MemVect_Len(:, :, cur_trial) = mem_len;
    rec_rnd_outfile.MemPeakHeight(:, :, cur_trial) = mem_peak_height;
    rec_rnd_outfile.MemPeakAngle(:, :, cur_trial) = mem_peak_ang;
    rec_rnd_outfile.MemPeakWidth(:, :, cur_trial) = mem_peak_width;
    rec_rnd_outfile.MemPeakProm(:, :, cur_trial) = mem_peak_prom;
    rec_rnd_outfile.MemVect_AngRel(:, :, cur_trial) = mem_ang_rel;
    rec_rnd_outfile.MemPeakAngleRel(:, :, cur_trial) = mem_peak_ang_rel;
    fprintf('Added statistics to file for trial %d. (%s)\n', cur_trial, datestr(now-start_time, 'HH:MM:SS'));

end %trial loop
    
%% Close the file
clear('rec_rnd_outfile');


