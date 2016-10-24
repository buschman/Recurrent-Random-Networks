function spockSimWM_MakeMovies(my_task_id)

cd('RecRnd Multi Simulations');

files = dir('*.mat');
files = {files(:).name};

tok = regexpi(files, 'SimWM_f(\d*)_TargetFR050_MaxFR010_Inputs(\d*).mat');
files = files(~cellfun(@isempty, tok));

SimWM_MakeMovie(files{my_task_id});
cd('..');