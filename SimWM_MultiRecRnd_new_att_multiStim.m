function SimWM_MultiRecRnd_new_att_multiStim(varargin),

% Simulates a recurrent 'working memory' network as well as randomly
% connected coupled network.

%% Parse optional inputs
opts.RunningPlot = 0; %whether or not to plot as we simulate

opts.SimTime = 2; %time to simulate, in seconds
opts.dt = 0.1/1000; %time step, in seconds
opts.Output_dt = 1/1000; %time step, in seconds, of the output
opts.b = 0; %bias in the firing response
opts.RecSynapseTau = 10/1000; %synapse tau, in seconds
opts.RndSynapseTau = 10/1000; %synapse tau, in seconds
opts.fI_Slope = 0.4; %slope of the non-linear f-I function
opts.RecActFxn = [];
opts.RndActFxn = [];
opts.InitSynapseRange = 0.01;

%Divisive normalization constants
opts.UseDivNorm = 0;
opts.MaxAvgRecFR = 20;
opts.MaxAvgRndFR = 20;
opts.DivNormTau = 5/1000; %in seconds

%Recurrent pool information
opts.N_rec_pools = 8; %# of recurrent neuron pools
opts.PoolWBaseline = 0; %baseline weight between recurrent pools (scaled for # of neurons)
opts.PoolWRandom = 0; %+/- range of random weights between recurrent pools (scaled for # of neurons)

opts.N_rec = 512; %# of recurrent neurons per pool
opts.RecWBaseline = 0.28; %baseline of weight matrix for recurrent network
opts.RecWAmp = 2; %amplitude of weight function
opts.RecWPositiveWidth = 1; %width of positive amplification
opts.RecWNegativeWidth = 0.25; %width of negative surround suppression

opts.N_rnd = 1024; %# of randomly-connected neurons
opts.RndWBaseline = 0; %Baseline inhibition between neurons in random network
opts.RndWSelf = 0; %Self-excitation in random network

opts.RndRec_f = 0.2; %likelihood of connection to/from any given recurrent neuron and a random neuron
opts.RecToRndW_NormalizeByNumConnections = 1; %if true, normalizes by the percentage of connections
opts.RecToRndW_TargetFR = 0.5; %firing rate needed in the recurrent network to drive random network to 20 Hz (scales weights)
opts.RecToRndW_EIBalance = -1;
opts.RecToRndW_Baseline = 0; %baseline weight from recurrent network to random network, regardless of connection existing

opts.RndToRecW_NormalizeByNumConnections = 1; %if true, normalizes by the percentage of connections
opts.RndToRecW_TargetFR = 0.25; %firing rate needed in the random network to drive recurrent network to 20 Hz (scales weights)
opts.RndToRecW_EIBalance = -1;
opts.RndToRecW_Baseline = 0; %baseline weight from recurrent network to random network, regardless of connection existing

%Input vector into recurrent network
opts.InputBaseline = 0; %strength of random inputs
opts.InputCenter = ceil(rand(1, opts.N_rec_pools)*opts.N_rec);
opts.InputWidth = round(opts.N_rec/32)*ones(1, opts.N_rec_pools);
opts.InputTime = repmat([0.1 0.2]', [1 opts.N_rec_pools]); %in seconds
opts.InputValue = 10*ones(1, opts.N_rec_pools); %strength of input (current injection)

%Attention input -- set to a random pool for now
opts.AttInputTime = NaN*ones(2, opts.N_rec_pools); opts.AttInputTime(:, randsample(opts.N_rec_pools, 1)) = [1 1.25];
opts.AttInputValue = 10*ones(1, opts.N_rec_pools); %strength of input (current injection)

opts.FileAppend = '_Selectivity';

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    %Check to see if this is an option for the adding to the processing
    %stream
    if isfield(opts, varargin{i}),
        opts.(varargin{i}) = varargin{i+1};
    else
        error(sprintf('Unknown parameter %s passed.', varargin{i}));
    end
end

%% Initialize network

%Define time
t = [0:opts.dt:opts.SimTime];

%Define activation functions
if isempty(opts.RecActFxn),
    %opts.ActFxn = @(g) 0.2/opts.SynapseTau*(1 + tanh(opts.fI_Slope*g - 3));
    opts.RecActFxn = @(g) 0.4/opts.RecSynapseTau*(1 + tanh(opts.fI_Slope*g - 3));
    %opts.RecActFxn = @(g) max(0, opts.fI_Slope*g - 0);
end
if isempty(opts.RndActFxn),
    %opts.RndActFxn = @(g) 2*0.2/opts.RndSynapseTau*heaviside(g - 6);%@(g) 0.2/opts.RndSynapseTau*(1 + tanh(opts.fI_Slope*g - 3));
    opts.RndActFxn = @(g) 0.4/opts.RndSynapseTau*(1 + tanh(opts.fI_Slope*g - 3));
    %opts.RndActFxn = @(g) max(0, opts.fI_Slope*g - 0);
end

%Define recurrent weight matrix, this will be identical for all recurrent pools
w_fun = @(x) opts.RecWBaseline + opts.RecWAmp*exp(opts.RecWPositiveWidth*(cos(x) - 1)) ...
    - opts.RecWAmp*exp(opts.RecWNegativeWidth*(cos(x) - 1));
ang = 2*pi*[1:opts.N_rec]./opts.N_rec;
W_rec = NaN*ones(opts.N_rec, opts.N_rec);
for i = 1:opts.N_rec,
    for j = 1:opts.N_rec,
        W_rec(i,j) = w_fun(ang(i) - ang(j));
        if i == j, W_rec(i,j) = 0; end
    end
end

%Define the inter-pool recurrent weight matrices
W_rec_rec = cell(opts.N_rec_pools, opts.N_rec_pools);
for i = 1:opts.N_rec_pools,
    for j = 1:opts.N_rec_pools,
        if i == j, continue; end
        W_rec_rec{i,j} = opts.PoolWBaseline + opts.PoolWRandom*2*(rand(opts.N_rec, opts.N_rec) - 0.5);
        W_rec_rec{i,j} = W_rec_rec{i,j}./(opts.N_rec.*(opts.N_rec_pools - 1));
    end
end

%Initialize synaptic activation for recurrent network(s)
s_rec = opts.InitSynapseRange*rand(opts.N_rec, opts.N_rec_pools);
a_rec = false(opts.N_rec, opts.N_rec_pools);
r_rec = NaN*ones(opts.N_rec, opts.N_rec_pools, length(t));

%Initialize random weight matrix -- should be symmetric (gets/sends signals
%to same neurons in recurrent network)
W_rec_rnd = opts.RndToRecW_Baseline*ones(opts.N_rec*opts.N_rec_pools, opts.N_rnd);
W_rnd_rec = opts.RecToRndW_Baseline*ones(opts.N_rnd, opts.N_rec*opts.N_rec_pools);
for i = 1:opts.N_rec*opts.N_rec_pools,
    %temp_con = (rand([1 opts.N_rnd]) <= opts.RndRec_f);
    temp_con = false(1, opts.N_rnd);
    temp_con(randsample(opts.N_rnd, round(opts.RndRec_f*opts.N_rnd), 0)) = true;
    
    W_rec_rnd(i, temp_con) = (1/0.1)/0.01*opts.RndToRecW_TargetFR./sum(temp_con); %this should give ~1 spike per input spike
    W_rec_rnd(i, :) = W_rec_rnd(i, :) + opts.RndToRecW_EIBalance*((1/0.1)/0.01*opts.RndToRecW_TargetFR)./opts.N_rnd;
    
    W_rnd_rec(temp_con, i) = (1/0.1)/0.01*opts.RecToRndW_TargetFR./sum(temp_con); %this should give ~1 spike per input spike
    W_rnd_rec(:, i) = W_rnd_rec(:, i) + opts.RecToRndW_EIBalance*((1/0.1)/0.01*opts.RecToRndW_TargetFR)./opts.N_rnd;
end %random neurons

%Recurrence within the random network
W_rnd = opts.RndWBaseline*ones(opts.N_rnd, opts.N_rnd);
W_rnd([0:(opts.N_rnd-1)]*(opts.N_rnd + 1) + 1) = opts.RndWSelf;

%Plot input matrices
if opts.RunningPlot,
    figure;
    subplot(2,2,1);
    imagesc(ang, ang, W_rec'); xlabel('Rec To'); ylabel('Rec From');
    colorbar;
    subplot(2,2,2);
    imagesc(ang, [1:opts.N_rnd], W_rec_rnd');xlabel('Rec To'); ylabel('Rnd From');
    colorbar;
    subplot(2,2,3);
    imagesc([1:opts.N_rnd], [1:opts.N_rnd], W_rnd'); xlabel('Rnd To'); ylabel('Rnd From');
    colorbar;
    subplot(2,2,4);
    imagesc([1:opts.N_rnd], ang, W_rnd_rec'); xlabel('Rnd To'); ylabel('Rec From');
    colorbar;
end

%Initialize synaptic activation for random
s_rnd = opts.InitSynapseRange*rand(opts.N_rnd, 1);
a_rnd = false(opts.N_rnd, 1);
r_rnd = NaN*ones(opts.N_rnd, length(t));

%% Initialize file to save data

opts.NumInputs = squeeze(sum(~any(isnan(opts.InputTime), 1)));
opts.NumStims = size(opts.InputCenter, 1);

my_save_fn = sprintf('SimWM_f%04.0f%s.mat', ...
    opts.RndRec_f*1000, opts.FileAppend);

if exist(my_save_fn, 'file') & ~overwrite,
    fprintf('File %s already exists. Skipping.\n', my_save_fn);
    return;
end

fprintf('Creating MAT file...\n');
out_t = [0:opts.Output_dt:opts.SimTime];
rec_rnd_outfile = matfile(my_save_fn, 'Writable', true);
rec_rnd_outfile.r_rec = NaN*ones(opts.N_rec, opts.N_rec_pools, length(out_t), opts.NumStims);
rec_rnd_outfile.r_rnd = NaN*ones(opts.N_rnd, length(out_t), opts.NumStims);
rec_rnd_outfile.t = out_t;
rec_rnd_outfile.opts = opts;

%% Simulate networks
start_time = now;

for cur_stim = 1:opts.NumStims,
    %% Initialize the networks
    
    %First, initialize synaptic activation for recurrent network(s)
    s_rec = opts.InitSynapseRange*rand(opts.N_rec, opts.N_rec_pools);
    a_rec = false(opts.N_rec, opts.N_rec_pools);
    r_rec = NaN*ones(opts.N_rec, opts.N_rec_pools, length(t));
    
    %Second, initialize synaptic activation for random
    s_rnd = opts.InitSynapseRange*rand(opts.N_rnd, 1);
    a_rnd = false(opts.N_rnd, 1);
    r_rnd = NaN*ones(opts.N_rnd, length(t));
    
    if opts.RunningPlot, line_fig = figure; end
    for t_ind = 1:length(t),
        %Update synaptic activation variables for recurrent networks
        for cp = 1:opts.N_rec_pools,
            dsdt = -s_rec(:, cp)./opts.RecSynapseTau + a_rec(:, cp)./opts.dt;
            s_rec(:, cp) = s_rec(:, cp) + dsdt*opts.dt;
        end
        %Update synaptic activation variables for random network
        dsdt = -s_rnd./opts.RndSynapseTau + a_rnd./opts.dt;
        s_rnd = s_rnd + dsdt*opts.dt;
        
        
        %Current for neurons
        g_rec = zeros(opts.N_rec, opts.N_rec_pools);
        for cp = 1:opts.N_rec_pools,
            %First, integrate my recurrent activity, the input from the random
            %network, and the baseline
            g_rec(:, cp) = W_rec*s_rec(:, cp) + W_rec_rnd((cp-1)*opts.N_rec + [1:opts.N_rec], :)*s_rnd + opts.b;
            
            %Now add the other recurrent networks
            for i = setdiff([1:opts.N_rec_pools], cp)
                g_rec(:, cp) = g_rec(:, cp) + W_rec_rec{cp, i}*s_rec(:, i);
            end
            
            %Finally, check for inputs

            %Sensory inputs
            if ~any(isnan(opts.InputTime(:, cp, cur_stim))) && (t(t_ind) >= opts.InputTime(1, cp, cur_stim)) && (t(t_ind) <= opts.InputTime(2, cp, cur_stim)),
                %Inputs arrive, centered at designated location
                inp_ind = [round(opts.InputCenter(cur_stim, cp) - 3*opts.InputWidth(cp)):round(opts.InputCenter(cur_stim, cp) + 3*opts.InputWidth(cp))];
                inp_scale = normpdf(inp_ind - opts.InputCenter(cur_stim, cp), 0, opts.InputWidth(cp));
                inp_scale = inp_scale./max(inp_scale);
                inp_ind = mod(inp_ind - 1, opts.N_rec) + 1; %wrap around vector
                %Create input vector, this is normalized so net input to
                %network is kept at zero
                inp_vect = zeros(opts.N_rec, 1); inp_vect(inp_ind) = opts.InputValue(cp)*inp_scale;
                inp_vect = inp_vect - sum(inp_vect)./length(inp_vect); %normalize so no net input to network
                
                g_rec(:, cp) = g_rec(:, cp) + inp_vect;
            end
            
            %Attention inputs
            if ~any(isnan(opts.AttInputTime(:, cp))) && (t(t_ind) >= opts.AttInputTime(1, cp)) && (t(t_ind) <= opts.AttInputTime(2, cp)),
                g_rec(:, cp) = g_rec(:, cp) + opts.AttInputValue(cp);
            end
        end
        %Add baseline, random input
        g_rec = g_rec + opts.InputBaseline*rand(size(g_rec));
        
        %Recurrence within random network, input from recurrent networks, and
        %baseline
        g_rnd = W_rnd*s_rnd + W_rnd_rec*s_rec(:) + opts.b;
        
        %Convert current to firing rate
        r_rec(:, :, t_ind) = opts.RecActFxn(g_rec);
        r_rnd(:, t_ind) = opts.RndActFxn(g_rnd);
        
        %Apply divisive normalization
        if opts.UseDivNorm,
            mean_r_rec = mean(r_rec(:, :, t_ind, 1));
            r_rec(:, :, t_ind) = r_rec(:, :, t_ind)./(opts.NormBaseline + repmat(mean_r_rec, [opts.N_rec 1]));
            mean_r_rnd = mean(r_rnd(:, t_ind), 1);
            r_rnd(:, t_ind) = r_rec(:, t_ind)./(opts.NormBaseline + repmat(mean_r_rnd, [opts.N_rnd 1]));
        end
        
        %Determine whether a spike occurred
        a_rec = (poissrnd(r_rec(:, :, t_ind)*opts.dt) >= 1);
        a_rnd = (poissrnd(r_rnd(:, t_ind)*opts.dt) >= 1);
        
        if opts.RunningPlot & mod(t_ind, 200) == 0,
            figure(line_fig);
            subplot(2,2,1); cla;
            plot(ang, r_rec(:, :, t_ind)); axis tight;
            set(gca, 'YLim', [0 80]);
            title(sprintf('Recurrent Network - %4.3f ms', t(t_ind)*1000));
            
            subplot(2,2,2); cla;
            plot([1:opts.N_rnd], r_rnd(:, t_ind)); axis tight;
            set(gca, 'YLim', [0 80]);
            title(sprintf('Random Network - %4.3f ms', t(t_ind)*1000));
            
            subplot(2,2,3); cla;
            plot([1:length(ang)], g_rec); hold all;
            plot(length(ang) + [1:length(s_rnd)], g_rnd);
            axis tight;
            title(sprintf('Recurrent Network - Synaptic Input', t(t_ind)*1000));
            
            subplot(2,2,4); cla;
            bar([1:(opts.N_rec_pools+1)], [squeeze(mean(r_rec(:, :, t_ind), 1)) mean(r_rnd(:, t_ind))]);
            
            drawnow;
            
            %         figure(time_fig);
            %         for i = 1:opts.N_rec_pools,
            %             subplot(2,6,mod(i-1,4)+1+6*floor((i-1)/4));
            %             imagesc(t, [1:size(r_rec, 1)]./size(r_rec, 1)*2*pi, squeeze(r_rec(:, i, :)));
            %             hold all;
            %             v = axis;
            %             plot(opts.InputTime(1)*[1 1], v(3:4), 'w-'); plot(opts.InputTime(2)*[1 1], v(3:4), 'w-');
            %             xlabel('Time'); ylabel('Angle');
            %             title(sprintf('Recurrent Network %d - Activity over Time', i));
            %             colorbar; set(gca, 'CLim', [0 80]);
            %         end
            %
            %         subplot(2,6,[5 6 11 12]);
            %         imagesc(t, [1:size(r_rnd, 1)], r_rnd);
            %         hold all;
            %         v = axis;
            %         plot(opts.InputTime(1)*[1 1], v(3:4), 'w-'); plot(opts.InputTime(2)*[1 1], v(3:4), 'w-');
            %         xlabel('Time'); ylabel('Neuron #');
            %         title('Random Network Activity over Time');
            %         colorbar; set(gca, 'CLim', [0 80]);
        end
    end %time loop
    
    %% Downsample the firing rate to desired output
    out_t = [0:opts.Output_dt:opts.SimTime];
    for i = 1:size(r_rec, 1),
        for j = 1:size(r_rec, 2),
            temp_r_rec(i, j, :) = reshape(interp1(t, squeeze(r_rec(i, j, :)), out_t), [1 1 length(out_t)]);
        end
    end
    r_rec = single(temp_r_rec); clear('temp_r_rec');
    for i = 1:size(r_rnd, 1),
        temp_r_rnd(i, :) = interp1(t, r_rnd(i, :), out_t);
    end
    r_rnd = single(temp_r_rnd); clear('temp_r_rnd');
    t = single(out_t);
    
    %% Write to disk
    rec_rnd_outfile.r_rec(:, :, :, cur_stim) = r_rec;
    rec_rnd_outfile.r_rnd(:, :, cur_stim) = r_rnd;
    
    %% Plot figure (if requested)
    
    if opts.RunningPlot,
        figure;
        for i = 1:opts.N_rec_pools,
            subplot(2,6,mod(i-1,4)+1+6*floor((i-1)/4));
            imagesc(t, [1:size(r_rec, 1)]./size(r_rec, 1)*2*pi, squeeze(r_rec(:, i, :)));
            hold all;
            v = axis;
            if ~any(isnan(opts.InputTime(:, i))),
                plot(opts.InputTime(1, i)*[1 1], v(3:4), 'w-'); plot(opts.InputTime(2, i)*[1 1], v(3:4), 'w-');
            end
            if ~any(isnan(opts.AttInputTime(:, i))),
                plot(opts.AttInputTime(1, i)*[1 1], v(3:4), 'r-'); plot(opts.AttInputTime(2, i)*[1 1], v(3:4), 'r-');
            end
            xlabel('Time'); ylabel('Angle');
            title(sprintf('Recurrent Network %d - Activity over Time', i));
            colorbar; set(gca, 'CLim', [0 80]);
        end
        
        subplot(2,6,[5 6 11 12]);
        imagesc(t, [1:size(r_rnd, 1)], r_rnd);
        hold all;
        v = axis;
        plot(opts.InputTime(1)*[1 1], v(3:4), 'w-'); plot(opts.InputTime(2)*[1 1], v(3:4), 'w-');
        plot(opts.AttInputTime(1)*[1 1], v(3:4), 'r-'); plot(opts.AttInputTime(2)*[1 1], v(3:4), 'r-');
        xlabel('Time'); ylabel('Neuron #');
        title('Random Network Activity over Time');
        colorbar; set(gca, 'CLim', [0 80]);
    end 
    
    fprintf('Simulated stimulation %d. (%s)\n', cur_stim, datestr(now-start_time, 'HH:MM:SS'));
end %stimulation loop