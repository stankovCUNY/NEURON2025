function [resFilt] = filterPseudoEMG(res)

    [EMGFromLFP] = bz_EMGFromLFP('/media/nasko/WD_BLACK3/RatU_temp2/');

    %% high correlation times 
    idxCutoff = find((EMGFromLFP.data.^2)>0.3); % r-squared cutoff above 0.6

    offset = EMGFromLFP.timestamps(idxCutoff);
    onset  = EMGFromLFP.timestamps(idxCutoff) - (1/EMGFromLFP.samplingFrequency);
    
    resBinFlags = zeros(length(res),1);
    
    for i = 1:length(idxCutoff)
        
        flagIdx = find((res > onset(i)) & (res < offset(i)));
        
        if isempty(flagIdx)
            continue;
        else
            resBinFlags(flagIdx) = 1; 
        end 
        
    end
    
    resFilt = res; 
    resFilt(find(resBinFlags)) = [];
    
end