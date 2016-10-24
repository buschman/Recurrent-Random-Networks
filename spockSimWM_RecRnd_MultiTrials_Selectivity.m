function spockSimWM_RecRnd_MultiTrials_Selectivity(my_task_id, num_stims, num_trials, overwrite)

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
% sim.RndRec_f = [0.150 0.200 0.250];
% sim.RecToRndW_TargetFR = [0.5 0.25];
% sim.MaxAvgRndFR = [10 20];
% sim.RecWPositiveWidth = [1 2];
sim.RndRec_f = [0.150 0.200 0.250];
sim.RecToRndW_TargetFR = [0.5];
sim.MaxAvgRndFR = [20];
sim.RecWPositiveWidth = [1];

length(sim.RndRec_f)*length(sim.RecToRndW_TargetFR)*length(sim.MaxAvgRndFR)*length(sim.RecWPositiveWidth)*num_trials

mysim.RndRec_f = sim.RndRec_f(mod(my_task_id, length(sim.RndRec_f))+1);
my_task_id = floor(my_task_id./length(sim.RndRec_f));
mysim.RecToRndW_TargetFR = sim.RecToRndW_TargetFR(mod(my_task_id, length(sim.RecToRndW_TargetFR))+1);
my_task_id = floor(my_task_id./length(sim.RecToRndW_TargetFR));
mysim.MaxAvgRndFR = sim.MaxAvgRndFR(mod(my_task_id, length(sim.MaxAvgRndFR))+1);
my_task_id = floor(my_task_id./length(sim.MaxAvgRndFR));
mysim.RecWPositiveWidth = sim.RecWPositiveWidth(mod(my_task_id, length(sim.RecWPositiveWidth))+1);
my_task_id = floor(my_task_id./length(sim.RecWPositiveWidth));
mysim.TrialNum = my_task_id + 1;

fprintf('Simulating a recurrent working memory network with random control units.\n');
fprintf('Parameters:\n');
fn = fieldnames(mysim);
for i = 1:length(fn),
    fprintf('%s: %s\n', fn{i}, mat2str(mysim.(fn{i})));
end

%% Run my simulation

tot_num_stims = 2*num_stims + num_stims*num_stims;

%Input times
cur_inp_time = repmat([0.1 0.3]', [1 4 tot_num_stims]); %in seconds
cur_inp_time(:, setdiff([1:4], 1), 1:num_stims) = NaN; %only drive the first cur_inp pools
cur_inp_time(:, setdiff([1:4], 2), (num_stims + 1):(2*num_stims)) = NaN; %only drive the first cur_inp pools
cur_inp_time(:, setdiff([1:4], [1:2]), (2*num_stims + 1):end) = NaN; %only drive the first cur_inp pools

%Baseline inputs should be structured for these simulations
cur_inp_center = repmat(ceil(rand(1, 4)*512), [tot_num_stims 1]);
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

%Run the same network simulation repeatedly
SimWM_MultiRecRnd_new_multiStim('RndRec_f', mysim.RndRec_f, 'RecToRndW_TargetFR', mysim.RecToRndW_TargetFR, ...
    'MaxAvgRecFR', mysim.MaxAvgRndFR, 'MaxAvgRndFR', mysim.MaxAvgRndFR, 'RecWPositiveWidth', mysim.RecWPositiveWidth, ...
    'InputTime', cur_inp_time, 'InputCenter', cur_inp_center, 'SimTime', 1, 'FileAppend', sprintf('_Selectivity_Network%02.0f', mysim.TrialNum));



