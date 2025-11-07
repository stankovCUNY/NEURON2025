function putativeCellType = cellClass(spikeTimes,troughToPeakLength,fs)
    % INPUTS
    % spikeTimes
    % troughToPeak duration
    %
    % OUTPUT
    % putativeCellType  a cell array with assigned cell types
    %
    % You may use this as a template for creating your own classification schemas
    
    % By Peter Petersen
    % petersen.peter@gmail.com
    % Last updated 14-05-2021
    % Modified by Atanas Stankov 2023
    
    % All cells are initially assigned as Pyramidal cells
    putativeCellType = 'p';

    % compute narrow ACGs for CellExplorer
    bins_narrow = 100;
    acg_narrow = CCG(spikeTimes,ones(size(spikeTimes)),'binSize',0.0005,'duration',0.100,'norm','rate','Fs',1/fs);
    
    plotFlag = false;
    fit_params_out = fit_ACG(acg_narrow,plotFlag);
    
    % Cells are reassigned as interneurons by below criteria 
    % Narrow interneuron assigned if troughToPeak <= 0.425 ms (preferences.putativeCellType.troughToPeak_boundary)
    if troughToPeakLength <= 0.425
        putativeCellType = 'i-narrow';
    end
    
    % acg_tau_rise > 6 ms (preferences.putativeCellType.acg_tau_rise_boundary) and troughToPeak > 0.425 ms
    if (troughToPeakLength > 0.425) && (fit_params_out.acg_tau_rise > 6)
        putativeCellType = 'i-wide';
    end
end