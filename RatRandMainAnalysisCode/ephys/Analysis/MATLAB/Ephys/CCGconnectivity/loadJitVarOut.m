function [jit_var_out,conn_type_idx] = loadJitVarOut(session_name,conn_type,jscale_name,alpha_name,datapath)

    if strcmp(session_name,'RoyMaze1') | strcmp(session_name,'RoyMaze1NULL')
        filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{1},jscale_name,alpha_name);
        % filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{1},5,alpha_name);
        load([datapath + filename])
        jit_var_out_gap = jit_var_out;

        filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{2},jscale_name,alpha_name);
        load([datapath + filename])
        jit_var_out_exc = jit_var_out;

        filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{3},jscale_name,alpha_name);
        if exist(filename,'file')
            load([datapath + filename])
            jit_var_out_inh = jit_var_out;

            clear jit_var_out
            jit_var_out = [jit_var_out_gap jit_var_out_exc jit_var_out_inh];
            
            conn_type_idx = [  ones(length(jit_var_out_gap),1); ...
                             2*ones(length(jit_var_out_exc),1); ...
                             3*ones(length(jit_var_out_inh),1)];
        else
            jit_var_out = [jit_var_out_gap jit_var_out_exc];
            
            conn_type_idx = [  ones(length(jit_var_out_gap),1); ...
                             2*ones(length(jit_var_out_exc),1)];
        end

        
                     
    elseif strcmp(session_name,'KevinMaze1')
        
        filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{1},jscale_name,alpha_name);
        % filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{1},5,alpha_name);
        load([datapath + filename])
        jit_var_out_gap = jit_var_out;

        filename = sprintf("%s_%s_jscale%d_alpha%d_jitter_stats.mat", session_name,conn_type{2},jscale_name,alpha_name);
        load([datapath + filename])
        jit_var_out_exc = jit_var_out;

        clear jit_var_out
        jit_var_out = [jit_var_out_gap jit_var_out_exc];

        conn_type_idx = [  ones(length(jit_var_out_gap),1); ...
                         2*ones(length(jit_var_out_exc),1)];
    end

end