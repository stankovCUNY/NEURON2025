function [pre_neurons, maze_neurons, post_neurons, cell_type] = saveAllNeuronData(spike_data_fullpath, varargin)
    
    spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    session_name = 'RoyMaze1';
    
    ip = inputParser;
    ip.addRequired('spike_data_fullpath', @isfile);
    ip.addRequired('session_name', @ischar);
    ip.addParameter('combine_epochs', false, @islogical)
    ip.parse(spike_data_fullpath, session_name, varargin{:});

    combine_epochs = ip.Results.combine_epochs;

    %%

    [data_dir, name, ~] = fileparts(spike_data_fullpath);

    load(spike_data_fullpath, 'spikes')

    load(fullfile(data_dir, 'wake-behavior.mat'), 'behavior');
    load(fullfile(data_dir, 'wake-basics.mat'),'basics');
    
    % comment out for now
%     load(fullfile(data_dir, 'sleep-behavior.mat'), 'behavior'); 
%     load(fullfile(data_dir, 'sleep-basics.mat'), 'basics');

    SampleRate = basics.(session_name).SampleRate;
    
    % Make data nicely formatted to work with buzcode
    % NOTE: this next line of code rejects poor quality neurons
    bz_spikes = Hiro_to_bz(spikes.(session_name), session_name);

    % Pull out PRE, MAZE, and POST time limits
    nepochs_plot = 3; % used to keep plots nice!
    if contains(name, 'wake')
        if ~combine_epochs
            nepochs = 3;
            save_name_append = '';
            epoch_names = {'Pre', 'Maze', 'Post'};
            time_list = behavior.(session_name).time/1000; % convert time list to milliseconds

        elseif combine_epochs
            nepochs = 1;
            epoch_names = {'Pre-Maze-Post Combined'};
            time_list = [behavior.(session_name).time(1)/1000, ...
                behavior.(session_name).time(end)/1000];
            save_name_append = '_combined';
        end
%     elseif contains(name, 'sleep')
%         nepochs = 3;
% 
%         % Name epochs nicely (up to 5 maximum!)
%     %     prefixes = {'First', 'Second', 'Third', 'Fourth', 'Fifth'};
%     %     epoch_names = cellfun(@(a) ['Sleep ' a ' ' num2str(1) '/' num2str(nepochs)], ...
%     %         prefixes(1:nepochs), 'UniformOutput', false);
%         epoch_names = arrayfun(@(a) ['Sleep Block ' num2str(a)], 1:nepochs, ...
%             'UniformOutput', false);
%         epoch_times = (0:nepochs)*diff(behavior.(session_name).time/1000)/nepochs + ...
%             behavior.(session_name).time(1)/1000;
%         save_name_append = ['_' num2str(nepochs) 'epochs'];
%         time_list = nan(nepochs,2);
%         for j = 1:nepochs
%             time_list(j,:) = [epoch_times(j), epoch_times(j+1)];
%         end
    end

    
    nneurons = length(spikes.(session_name));
    for j = 1:nepochs
        epoch_bool = bz_spikes.spindices(:,1) >= time_list(j,1) ...
            & bz_spikes.spindices(:,1) <= time_list(j,2); % ID spike times in each epoch
        parse_spikes(j).spindices = bz_spikes.spindices(epoch_bool,:); % parse spikes by epoch into this variable
    end
    % Figure out if pyramidal or inhibitory
    cell_type = repmat('p', 1, length(bz_spikes.quality));
    cell_type(bz_spikes.quality == 8) = 'i';
    
    pre  = parse_spikes(1).spindices;
    maze = parse_spikes(2).spindices;
    post = parse_spikes(3).spindices;
    
    pre_neurons  = {};
    maze_neurons = {};
    post_neurons = {};
    
    for j = 1:length(bz_spikes.UID) 
        
        pre_neurons.spikeTimes{j}  = pre(find(pre(:,2)   == bz_spikes.UID(j)),1);
        maze_neurons.spikeTimes{j} = maze(find(maze(:,2) == bz_spikes.UID(j)),1);
        post_neurons.spikeTimes{j} = post(find(post(:,2) == bz_spikes.UID(j)),1);
    
        pre_neurons.UID(j)  = bz_spikes.UID(j);
        maze_neurons.UID(j) = bz_spikes.UID(j);
        post_neurons.UID(j) = bz_spikes.UID(j);
        
        pre_neurons.shankID(j)  = bz_spikes.shankID(j);
        maze_neurons.shankID(j) = bz_spikes.shankID(j);
        post_neurons.shankID(j) = bz_spikes.shankID(j);
        
    end
    
end