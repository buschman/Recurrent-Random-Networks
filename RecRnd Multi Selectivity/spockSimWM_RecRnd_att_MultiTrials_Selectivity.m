function spockSimWM_RecRnd_att_MultiTrials_Selectivity(my_task_id, num_stims, num_trials, overwrite)

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
sim.RndRec_f = [0.150 0.200 0.350];
sim.W_ratio = [3/4 1 4/3];

length(sim.RndRec_f)*length(sim.W_ratio)*num_trials

mysim.RndRec_f = sim.RndRec_f(mod(my_task_id, length(sim.RndRec_f))+1);
my_task_id = floor(my_task_id./length(sim.RndRec_f));
mysim.W_ratio = sim.W_ratio(mod(my_task_id, length(sim.W_ratio))+1);
my_task_id = floor(my_task_id./length(sim.W_ratio));
mysim.TrialNum = my_task_id;

fprintf('Simulating a recurrent working memory network with random control units.\n');
fprintf('Parameters:\n');
fn = fieldnames(mysim);
for i = 1:length(fn),
    fprintf('%s: %s\n', fn{i}, mat2str(mysim.(fn{i})));
end

%% Run my simulation

N_rec_pools = 8;

tot_num_stims = 2*num_stims + num_stims*num_stims;

%Input times
cur_inp_time = repmat([0.1 0.3]', [1 N_rec_pools tot_num_stims]); %in seconds
cur_inp_time(:, setdiff([1:N_rec_pools], 1), 1:num_stims) = NaN; %only drive the first cur_inp pools
cur_inp_time(:, setdiff([1:N_rec_pools], 2), (num_stims + 1):(2*num_stims)) = NaN; %only drive the first cur_inp pools
cur_inp_time(:, setdiff([1:N_rec_pools], [1:2]), (2*num_stims + 1):end) = NaN; %only drive the first cur_inp pools

%Baseline inputs should be structured for these simulations
cur_inp_center = repmat(ceil(rand(1, N_rec_pools)*512), [tot_num_stims 1]);
cur_inp_center(1:num_stims, 1) = ceil([1:num_stims]./num_stims*512);
cur_inp_center((num_stims + 1):(2*num_stims), 2) = ceil([1:num_stims]./num_stims*512);
count = 2*num_stims;
for i = 1:num_stims,
    for j = 1:num_stims,
        count = count + 1;
        cur_inp_center(count, 1) = ceil(i./num_stims*512);
        cur_inp_center(count, 2) = ceil(j./num_stims*512);
    end
end

cur_att_time = NaN*ones(2, N_rec_pools);
cur_att_inp = 3*[7 -1 -1 -1 -1 -1 -1 -1];

%Run the same network simulation repeatedly
SimWM_MultiRecRnd_new_att_multiStim('SimTime', 1, 'RunningPlot', 0, ...
    'RndRec_f', mysim.RndRec_f, 'InputTime', cur_inp_time, 'InputCenter', cur_inp_center, 'N_rnd', 1024, ...
    'RecToRndW_TargetFR', 0.3.*mysim.W_ratio, 'RecToRndW_EIBalance', -1.02, ...
    'RndToRecW_TargetFR', 0.2.*mysim.W_ratio, 'RndToRecW_EIBalance', -1.02, ...
    'AttInputTime', cur_att_time, 'AttInputValue', cur_att_inp, ...
    'FileAppend', sprintf('_Wratio%s_Selectivity_Network%02.0f', strtrim(strrep(rats(mysim.W_ratio), '/', '-')), mysim.TrialNum));
