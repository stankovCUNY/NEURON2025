load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikes.mat')

spikesNULL = spikes;

for loopUnits = 1:length({spikes.RoyMaze1.time})
    
    spikeTimes = spikes.RoyMaze1(loopUnits).time/1e6; % convert to seconds from microseconds
    
    one_ms = 0.001;
    delta = 5*one_ms;
    spikeTimes = floor(spikeTimes/delta)*delta + rand(size(spikeTimes))*delta;
    
    spikesNULL.RoyMaze1(loopUnits).time = spikeTimes*1e6; % convert to microseconds from seconds
    
end 

save('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/wake-spikesNULL.mat','spikesNULL')

%% sanity check

% jscale         = 1;
% alpha_name     = 5;
% duration       = 0.007;
% fs             = 30000;
% binSize        = 1/fs;
% fig_use        = 102;
% njitter        = 500;
% alpha          = 0.05;
% for_grant      = false;
% plotFlag       = true;
% resolution_use = '-r300'; % >= 300dpi required by Nature Neuro, so that's what I'm using.
% 
% 
% ref = spikes.RoyMaze1(loopUnits).time;
% tar = spikesNULL.RoyMaze1(loopUnits).time;
% 
%  [ccgR, tR] = CCG([ref';tar'],[ones(size(ref'));2*ones(size(tar'))], ...
%                 'binSize', binSize, 'duration', duration, 'Fs', 1/fs,...
%                 'norm', 'counts');
%                  line(tR*1e3,ccgR(:,1,2)/length(ref),'color','k', 'LineWidth',1)