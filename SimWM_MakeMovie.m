function SimWM_MakeMovie(inp_filestr, varargin),

%% Process optional inputs
opts.DownSample = 1; %whether to downsample in time
opts.SimNumToPlot = 1; %which simulation run to plot

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

%% Loop through files, making movies

if ~iscell(inp_filestr),
    inp_filestr = {inp_filestr};
end

for cur_file = 1:length(inp_filestr),
    
    inp = matfile(inp_filestr{cur_file});
    
    ang = [1:size(inp, 'r_rec', 1)]./size(inp, 'r_rec', 1)*2*pi;
    
    clear F;
    for t_ind = 1:opts.DownSample:size(inp, 'r_rec', 3),
        
        frame_fig = figure;
        for cp = 1:4,
            subplot(2,4,mod(cp-1,2)+1+4*floor((cp-1)/2)); colormap hot;
            cmap = colormap;
            lh(cp) = plot3(sin(ang), cos(ang), squeeze(inp.r_rec(:, cp, t_ind, opts.SimNumToPlot)), 'LineWidth', 1);
            axis off; set(gca, 'ZLim', [0 80]);
            camzoom(1.5)
            drawnow;
            text(0, 0, 0, sprintf('Sensory Network %d', cp), 'HorizontalAlignment', 'center');
        end
        
        subplot(2,4,[3 4 7 8]);
        bh = bar3(reshape(inp.r_rnd(:, t_ind, opts.SimNumToPlot), [32 64]));
        set(gca, 'ZLim', [0 80], 'CLim', [0 40]);
        for k = 1:length(bh)
            zdata = bh(k).ZData;
            bh(k).CData = zdata;
            bh(k).FaceColor = 'interp';
        end
        axis off;
        camzoom(1.4);
        set(gcf, 'Position', [70 20 1200 600]);
        get(gca, 'XLim')
        get(gca, 'YLim')
        text(max(get(gca, 'XLim'))/4, 0, 80, sprintf('Random Control Network\n     [% 6.1f ms]', inp.t(1, t_ind)*1000), 'FontSize', 24);
        %camorbit(35, 0);
        
        for cp = 1:4,
            col = ones(4, size(inp, 'r_rec', 1));
            for i = 1:3,
                col(i, :) = interp1(linspace(0, 80, size(cmap, 1)), cmap(:, i), inp.r_rec(:, cp, t_ind, opts.SimNumToPlot));
            end
            set(lh(cp).Edge, 'ColorType', 'truecolor', 'ColorData', uint8(round(col*255)), 'ColorBinding', 'interpolated');
        end
        
        F(t_ind) = getframe(frame_fig);
        close(frame_fig);
    end %time loop
    
    %Save movie to file
    [pathstr, inp_filename] = fileparts(inp_filestr{1});
    myVideo = VideoWriter(sprintf('%s%s_Run%02.0f.avi', pathstr inp_filename opts.SimNumToPlot));
    myVideo.FrameRate = 30;
    open(myVideo);
    writeVideo(myVideo, F);
    close(myVideo);
    
    clear('myVideo', 'F', 'inp');
    
end %file loop