%% Part of pre_v_postCCG_v2.m

%% Send all variables calculated from jitter into data structure.
function save_jitter_vars(jit_var_in, cell_pair, cell_spike_times, row, epoch,data_dir, ...
    epoch_names, session_name, screen_type, jscale, alpha, varargin)
    var_names = {'GSPExc','GSPInh','pvalE','pvalI','ccgR','tR','LSPExc',...
        'LSPInh','JBSIE','JBSII','cell1type','cell2type','cell1shank','cell2shank'};

    jit_var_out = jit_var_in;
    jit_var_out(1, 1).cell_pair = cell_pair;
    
    for var_num = 1:length(var_names)
        jit_var_out(1, 1).(var_names{var_num}) = varargin{var_num};
    end
    jit_var_out(1, 1).epoch = epoch_names{epoch};
    
    jit_var_out.cell1_spike_times = cell_spike_times{1};
    jit_var_out.cell2_spike_times = cell_spike_times{2};
    
    tempFilename = [session_name '_' screen_type '_jscale' num2str(jscale) ...
                    '_alpha' num2str(round(alpha*100)) '_jitter_stats.mat'];
    
    if ~isfile(data_dir + "/" + tempFilename)
        save(data_dir + "/" + tempFilename,'jit_var_out')
    else     
        temp = load(data_dir + "/" + tempFilename);
        jit_var_out = [temp.jit_var_out jit_var_out];
        save([data_dir + "/" + tempFilename],'jit_var_out')
    end
    
end