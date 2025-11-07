clear
% 
% load('jitter_stats.mat')
% negLagBound = 76; % -1ms
% posLagBound = 136; % +1ms

load('jitter_stats_coarse.mat')
jit_var_out = jit_var_out(193:2811);
negLagBound = 1; % -1ms
posLagBound = 41; % +1ms


for i = 1:size(jit_var_out,2)
    
    PreIdx(i)  = strcmp(jit_var_out(i).epoch,'Pre');
    MazeIdx(i) = strcmp(jit_var_out(i).epoch,'Maze');
    PostIdx(i) = strcmp(jit_var_out(i).epoch,'Post');
    
    pyr_pyr(i) = strcmp([jit_var_out(i).cell1type jit_var_out(i).cell2type],'pp');
    pyr_int(i) = strcmp([jit_var_out(i).cell1type jit_var_out(i).cell2type],'ip') || ...
                 strcmp([jit_var_out(i).cell1type jit_var_out(i).cell2type],'pi');
    int_int(i) = strcmp([jit_var_out(i).cell1type jit_var_out(i).cell2type],'ii');
    
end

pyr_pyr_pre  = PreIdx  & pyr_pyr;
pyr_pyr_maze = MazeIdx & pyr_pyr;
pyr_pyr_post = PostIdx & pyr_pyr;

pyr_int_pre  = PreIdx  & pyr_int;
pyr_int_maze = MazeIdx & pyr_int;
pyr_int_post = PostIdx & pyr_int;

int_int_pre  = PreIdx  & int_int;
int_int_maze = MazeIdx & int_int;
int_int_post = PostIdx & int_int;

pyr_pyr_pre  = find(pyr_pyr_pre);
pyr_pyr_maze = find(pyr_pyr_maze);
pyr_pyr_post = find(pyr_pyr_post);

pyr_int_pre  = find(pyr_int_pre);
pyr_int_maze = find(pyr_int_maze);
pyr_int_post = find(pyr_int_post);

int_int_pre  = find(int_int_pre);
int_int_maze = find(int_int_maze);
int_int_post = find(int_int_post);


posSigBins  = zeros(size(jit_var_out,2),1);
negSigBins  = zeros(size(jit_var_out,2),1);

for i = 1:size(jit_var_out,2)
    
    posSigBins(i)   = sum(jit_var_out(i).GSPExc(negLagBound:posLagBound));
    negSigBins(i)   = sum(jit_var_out(i).GSPInh(negLagBound:posLagBound));
 
end

pyr_pyr_pre_posBins  = posSigBins(pyr_pyr_pre);  pyr_pyr_pre_negBins  = negSigBins(pyr_pyr_pre);
pyr_pyr_maze_posBins = posSigBins(pyr_pyr_maze); pyr_pyr_maze_negBins = negSigBins(pyr_pyr_maze);
pyr_pyr_post_posBins = posSigBins(pyr_pyr_post); pyr_pyr_post_negBins = negSigBins(pyr_pyr_post);

pyr_int_pre_posBins  = posSigBins(pyr_int_pre);  pyr_int_pre_negBins  = negSigBins(pyr_int_pre);
pyr_int_maze_posBins = posSigBins(pyr_int_maze); pyr_int_maze_negBins = negSigBins(pyr_int_maze);
pyr_int_post_posBins = posSigBins(pyr_int_post); pyr_int_post_negBins = negSigBins(pyr_int_post);

int_int_pre_posBins  = posSigBins(int_int_pre);  int_int_pre_negBins  = negSigBins(int_int_pre);
int_int_maze_posBins = posSigBins(int_int_maze); int_int_maze_negBins = negSigBins(int_int_maze);
int_int_post_posBins = posSigBins(int_int_post); int_int_post_negBins = negSigBins(int_int_post);

%%

subplot(3,3,1)
scatter(pyr_pyr_pre_posBins, pyr_pyr_pre_negBins,'b')
title('Pre Pyr-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,2)
scatter(pyr_int_pre_posBins, pyr_int_pre_negBins,'g')
title('Pre Pyr-Int/Int-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,3)
scatter(int_int_pre_posBins, int_int_pre_negBins,'r')
title('Pre')
title('Pre Int-Int')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

%%

subplot(3,3,4)
scatter(pyr_pyr_maze_posBins, pyr_pyr_maze_negBins,'b')
title('Maze Pyr-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,5)
scatter(pyr_int_maze_posBins, pyr_int_maze_negBins,'g')
title('Maze Pyr-Int/Int-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,6)
scatter(int_int_maze_posBins, int_int_maze_negBins,'r')
title('Maze Int-Int')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

%%

subplot(3,3,7)
scatter(pyr_pyr_post_posBins, pyr_pyr_post_negBins,'b')
title('Post Pyr-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,8)
scatter(pyr_int_post_posBins, pyr_int_post_negBins,'g')
title('Post Pyr-Int/Int-Pyr')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

subplot(3,3,9)
scatter(int_int_post_posBins, int_int_post_negBins,'r')
title('Post Int-Int')
xlabel('neg significant bins')
ylabel('pos significant bins')
xlim([0 10])
ylim([0 10])

%%