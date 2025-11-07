ExcPairsStr = load('RoyMaze1_GapPairs_jscale1_alpha5_jitter_stats.mat');
InhPairsStr = load('RoyMaze1_ExcPairs_jscale1_alpha5_jitter_stats.mat');
GapPairsStr = load('RoyMaze1_InhPairs_jscale1_alpha5_jitter_stats.mat');

ExcPairs = [];
InhPairs = [];
GapPairs = [];

for i = 1:size(ExcPairsStr.jit_var_out,2) 
    ExcPairs = [ExcPairs; ExcPairsStr.jit_var_out(i).cell_pair];
end

for i = 1:size(InhPairsStr.jit_var_out,2)
    InhPairs = [InhPairs; InhPairsStr.jit_var_out(i).cell_pair];
end

for i = 1:size(GapPairsStr.jit_var_out,2)
    GapPairs = [GapPairs; GapPairsStr.jit_var_out(i).cell_pair];
end

% [~,idxUnique,idxNonUnique] = unique(ExcPairs,'rows');
% ExcPairs(idxUnique,:) = [];
% 
% [~,idxUnique,~] = unique(InhPairs,'rows');
% InhPairs(idxUnique,:) = [];
% 
% [~,idxUnique,~] = unique(GapPairs,'rows');
% GapPairs(idxUnique,:) = [];

allPairs = [ExcPairs; InhPairs; GapPairs];
[~,idxUnique,~] = unique(allPairs,'rows');
allPairs = allPairs(idxUnique,:);
