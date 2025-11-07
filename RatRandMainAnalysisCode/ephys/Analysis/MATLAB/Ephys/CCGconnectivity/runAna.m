clear all; close all;

nullFlag = true;

if ~exist("/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat",'file')
    
    mex -O C:\Users\orest\Documents\GitHub\SuperfineSynchrony\buzcode\analysis\spikes\correlation\CCGHeart.c
    spike_data_fullpath = "C:\Users\orest\Documents\CCNY Related Stuff\SCRP Fellowship\GroupData\wake_new\wake-spikes.mat";

else
        % add functions
    % addpath("/home/nasko/CUNY_Work_Han_Kamran_Nat/buzcode/analysis/spikes/correlation")
    % addpath("/home/nasko/CUNY_Work_Han_Kamran_Nat/ephys/Analysis/MATLAB/Ephys/CCGconnectivity")
    % addpath("/home/nasko/CUNY_Work_Han_Kamran_Nat/ephys/Analysis/MATLAB/Ephys/CCGconnectivity/one_off_functions")
    % addpath("/home/nasko/CUNY_Work_Han_Kamran_Nat/ephys/Analysis/MATLAB/Ephys/CCGconnectivity/scripts")
    % addpath("/home/nasko/CUNY_Work_Han_Kamran_Nat/buzcode/externalPackages/FMAToolbox/Helpers")

    % compile CCG C function
    % mex -O /home/nasko/CUNY_Work_Han_Kamran_Nat/buzcode/analysis/spikes/correlation/CCGHeart.c
    
    if nullFlag
        spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikesNULL.mat";
    else
        spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat";
    end
end

% run wake 
session_name        = 'RoyMaze1';
% session_name        = 'KevinMaze1';
jscaleVal = 1;

pair_types = {'GapPairs','ExcPairs','InhPairs'};

% mode = "conv";
mode = "jitter";

if strcmp(mode,"conv")
    % [pairs_plot_all, pp_pairs] = pre_v_postCCG(spike_data_fullpath, session_name); % jscale 5
    % [pairs_plot_all, pp_pairs] = pre_v_postCCG(spike_data_fullpath, session_name,'jscale',1);
    
    [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal,'conn_type', 'GapPairs','nullFlag',nullFlag);
    [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal,'conn_type', 'ExcPairs','nullFlag',nullFlag);
    [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal,'conn_type', 'InhPairs','nullFlag',nullFlag);
    
    % run sleep 
    % session_name        = 'RoySleep0';
    % spike_data_fullpath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/sleep/sleep-spikes.mat";
    % [pairs_plot_all, pp_pairs] = pre_v_postCCG(spike_data_fullpath, session_name);
elseif strcmp(mode,"jitter")
%     [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal, ... 
%         'plot_jitter',true,'plot_conv',false,'save_jitter',true,'screen_type','two_prong');
     [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal, ... 
        'plot_jitter',true,'plot_conv',false,'save_jitter',true,'conn_type', 'GapPairs','nullFlag',nullFlag);
     [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal, ... 
        'plot_jitter',true,'plot_conv',false,'save_jitter',true,'conn_type', 'ExcPairs','nullFlag',nullFlag);
     [pairs_plot_all, pp_pairs] = pre_v_postCCG_v2(spike_data_fullpath, session_name,'jscale',jscaleVal, ... 
        'plot_jitter',true,'plot_conv',false,'save_jitter',true,'conn_type', 'InhPairs','nullFlag',nullFlag);
end


% If Atanas' Computer and convert ps to pdf
if ~exist("C:\Users\orest\Documents\CCNY Related Stuff\SCRP Fellowship\GroupData\wake_new\",'dir')
    
    for i = 1:3
        datapath = "/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/";
        [status,result] =  system("ps2pdf " + ...
            datapath + sprintf("RoyMaze1_all_%s_jscale%d_coarse_CCGs_%s.ps ",pair_types{i},jscaleVal,mode) + ... 
            datapath + sprintf("RoyMaze1_all_%s_jscale%d_coarse_CCGs_%s.pdf ",pair_types{i},jscaleVal,mode));
        [status,result] =  system("ps2pdf " + ...
            datapath + sprintf("RoyMaze1_all_%s_jscale%d_fine_CCGs_%s.ps ",pair_types{i},jscaleVal,mode) + ... 
            datapath + sprintf("RoyMaze1_all_%s_jscale%d_fine_CCGs_%s.pdf ",pair_types{i},jscaleVal,mode));
    end

end
