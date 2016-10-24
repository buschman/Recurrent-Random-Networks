function plotSimWM_MultiRecRnd_new_NetworkStructure(varargin),

% Plots the network structure used when simulating a recurrent 'working memory' network as well as randomly
% connected coupled network.

%% Parse optional inputs
opts.RunningPlot = 0; %whether or not to plot as we simulate

opts.SimTime = 0.5; %time to simulate, in seconds
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
opts.MaxAvgRecFR = 10;
opts.MaxAvgRndFR = 10;
opts.DivNormTau = 5/1000; %in seconds

%Recurrent pool information
opts.N_rec_pools = 4; %# of recurrent neuron pools
opts.PoolWBaseline = 0; %baseline weight between recurrent pools (scaled for # of neurons)
opts.PoolWRandom = 0; %+/- range of random weights between recurrent pools (scaled for # of neurons)

opts.N_rec = 512; %# of recurrent neurons per pool
opts.RecWBaseline = 0.3; %baseline of weight matrix for recurrent network
opts.RecWAmp = 2; %amplitude of weight function
opts.RecWPositiveWidth = 1; %width of positive amplification
opts.RecWNegativeWidth = 0.25; %width of negative surround suppression

opts.N_rnd = 2048; %# of randomly-connected neurons
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
    W_rnd_rec(:, i) = W_rnd_rec(:, i) + opts.RecToRndW_EIBalance*((1/0.1)/0.01*opts.RndToRecW_TargetFR)./opts.N_rnd;
end %random neurons

%Recurrence within the random network
W_rnd = opts.RndWBaseline*ones(opts.N_rnd, opts.N_rnd);
W_rnd([0:(opts.N_rnd-1)]*(opts.N_rnd + 1) + 1) = opts.RndWSelf;

%Initialize synaptic activation for random
s_rnd = opts.InitSynapseRange*rand(opts.N_rnd, 1);
a_rnd = false(opts.N_rnd, 1);
r_rnd = NaN*ones(opts.N_rnd, length(t));

%% Make network plots

%Plot the weigth function for a single neuron's connectivity in the
%recurrent, sensory network
figure;
plot([1:opts.N_rec] - round(opts.N_rec/2), W_rec(round(opts.N_rec/2), :));
xlabel('Neuron Position, Relative to Current Neuron');
ylabel('Synaptic Weight');
hold all;
set(gca, 'XLim', [1 opts.N_rec] - round(opts.N_rec/2), 'YLim', [-0.7 0.35]);
v = axis;
plot(v(1:2), [0 0], 'k-');
saveas(gcf, 'SensoryRecurrentNetwork_ConnectivityProfile.fig');
saveas(gcf, 'SensoryRecurrentNetwork_ConnectivityProfile.eps', 'psc2');

%Plot input matrices
figure;
subplot(2,2,1);
imagesc([1:opts.N_rec], [1:opts.N_rec], W_rec'); xlabel('Rec To'); ylabel('Rec From');
colorbar;
subplot(2,2,2);
imagesc([1:opts.N_rec], [1:opts.N_rnd], W_rec_rnd');xlabel('Rec To'); ylabel('Rnd From');
colorbar;
subplot(2,2,3);
imagesc([1:opts.N_rnd], [1:opts.N_rnd], W_rnd'); xlabel('Rnd To'); ylabel('Rnd From');
colorbar;
subplot(2,2,4);
imagesc([1:opts.N_rnd], [1:opts.N_rec], W_rnd_rec'); xlabel('Rnd To'); ylabel('Rec From');
colorbar;
saveas(gcf, 'SensoryRecurrentNetwork_FullConnectivityMatrix.fig');
saveas(gcf, 'SensoryRecurrentNetwork_FullConnectivityMatrix.eps', 'psc2');


%% Plot the recurrent networks and random network
draw_skip = opts.N_rec/16;
neuron_r = 10;
figure;

%Radius of the pools
ovr_r = ((2.5*neuron_r)*floor(opts.N_rec/draw_skip))/(2*pi);

%Calculate position of random neurons
num_rows = round((4*ovr_r + 7*neuron_r)./(2.5*neuron_r));
mid_point_y = (2*ovr_r + 5*neuron_r)/2;
count = 0;
for j = 1:draw_skip:opts.N_rnd,
    count = count + 1;
    rnd_neuron_x_pos(j) = floor((count-1)./num_rows)*(2.5*neuron_r) + 2*(2*ovr_r + 5*neuron_r);
    rnd_neuron_y_pos(j) = (mod((count-1), num_rows) - (num_rows-1)/2)*(2.5*neuron_r) + mid_point_y;
end

angs = [0:(pi/20):(2*pi)];
%Draw each pool in turn
for i = 1:opts.N_rec_pools,
    
    %Center of this pool
    pool_x_pos = floor((i-1)/2)*(2*ovr_r + 5*neuron_r);
    pool_y_pos = mod(i-1, 2)*(2*ovr_r + 5*neuron_r);
    for j = 1:draw_skip:opts.N_rec,
        %Center of this neuron
        neuron_x_pos(j) = pool_x_pos + sin(((j-1)/opts.N_rec)*2*pi)*ovr_r;
        neuron_y_pos(j) = pool_y_pos + cos(((j-1)/opts.N_rec)*2*pi)*ovr_r;
    end
    
    %Draw connections with recurrent, sensory networks
    min_W = min(W_rec(:)); max_W = max(W_rec(:));
    set(gca, 'CLim', [min_W max_W]);
    cmap = colormap(gcf, 'cool');
    for j = 1:draw_skip:opts.N_rec,
        if j == 1, continue; end %this is our source
        plot([neuron_x_pos(1) neuron_x_pos(j)], [neuron_y_pos(1) neuron_y_pos(j)], 'LineWidth', abs(W_rec(1, j)*4), 'Color', cmap(round((W_rec(1, j) - min_W)/(max_W - min_W)*size(cmap, 1))+1, :));
        hold on;
    end
    
    %Draw circles on top of recurrent sensory networks
    for j = 1:draw_skip:opts.N_rec,
        plot(neuron_x_pos(j) + neuron_r*sin(angs), neuron_y_pos(j) + neuron_r*cos(angs), 'Color', [0 0 0]);
        %viscircles([neuron_x_pos(j) neuron_y_pos(j)], neuron_r, 'EdgeColor', [0 0 0], 'DrawBackgroundCircle', true);
    end
    
    %Plot connections to a few random neurons
    for j = 1:draw_skip:opts.N_rec,
        cur_rnd = 1;
        if W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec) > 0,
            plot([rnd_neuron_x_pos(cur_rnd) neuron_x_pos(j)], [rnd_neuron_y_pos(cur_rnd) neuron_y_pos(j)], 'LineWidth', abs(W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec)*4), 'Color', [0.7 0.7 0.7]);
        end
        cur_rnd = (num_rows-1)*draw_skip + 1;
        if W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec) > 0,
            plot([rnd_neuron_x_pos(cur_rnd) neuron_x_pos(j)], [rnd_neuron_y_pos(cur_rnd) neuron_y_pos(j)], 'LineWidth', abs(W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec)*4), 'Color', [0.7 0.7 0.7]);
        end
                cur_rnd = round((num_rows-1)/2*draw_skip + 1);
        if W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec) > 0,
            plot([rnd_neuron_x_pos(cur_rnd) neuron_x_pos(j)], [rnd_neuron_y_pos(cur_rnd) neuron_y_pos(j)], 'LineWidth', abs(W_rnd_rec(cur_rnd, j + (i-1)*opts.N_rec)*4), 'Color', [0.7 0.7 0.7]);
        end
    end
end %recurrent loop

%Draw circles on top of random control network
for j = 1:draw_skip:opts.N_rnd,
    plot(rnd_neuron_x_pos(j) + neuron_r*sin(angs), rnd_neuron_y_pos(j) + neuron_r*cos(angs), 'Color', [0.4 0.4 0.4]);
    %viscircles([rnd_neuron_x_pos(j) rnd_neuron_y_pos(j)], neuron_r, 'EdgeColor', [0.7 0.7 0.7], 'DrawBackgroundCircle', true);
end

set(gca, 'YLim', [-(ovr_r + 2*neuron_r) (3*ovr_r + 7*neuron_r)]);
axis image; axis off;
colorbar ('horiz');
saveas(gcf, 'SensoryRecurrentNetwork_ExampleConnectivity.fig');
saveas(gcf, 'SensoryRecurrentNetwork_ExampleConnectivity.svg');