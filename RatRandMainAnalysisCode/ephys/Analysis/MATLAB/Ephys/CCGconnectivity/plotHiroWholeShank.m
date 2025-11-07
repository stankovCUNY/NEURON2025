tWave = (-35:36)*(1/30);

channelVectorMap =      [1,8,9 ,16,17,24,25,32,33,40,41,48,49,56,57,64, ...
                         2,7,10,15,18,23,26,31,34,39,42,47,50,55,58,63, ...
                         3,6,11,14,19,22,27,30,35,38,43,46,51,54,59,62, ...
                         4,5,12,13,20,21,28,29,36,37,44,45,52,53,60,61];

channelSnankVectorMap = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, ...
                         1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, ...
                         1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, ...
                         1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8];
                
tiledlayout(4,16)

for i = 1: length(channelVectorMap)
    
   nexttile(i)
   plot(tWave,spikeAvgMaxChanCell1(channelVectorMap(i),:),'color',"#0072BD")
   hold on 
   plot(tWave,spikeAvgMaxChanCell2(channelVectorMap(i),:),'color',"#D95319")
%    plot(tWave,spikeAvgMaxChanCell1(channelVectorMap(i),:)+spikeAvgMaxChanCell2(channelVectorMap(i),:),'color',"k")
   hold off
   
   if channelVectorMap(i) == chCell1 
       title(['ch:' num2str(channelVectorMap(i)) ', sh: ' num2str(channelSnankVectorMap(i))],'color',"#0072BD")
   elseif channelVectorMap(i) == chCell2 
       title(['ch:' num2str(channelVectorMap(i)) ', sh: ' num2str(channelSnankVectorMap(i))],'color',"#D95319")
   else
       title(['ch:' num2str(channelVectorMap(i)) ', sh: ' num2str(channelSnankVectorMap(i))])
   end
   
   legend([num2str(cellID1) cell1type],[num2str(cellID2) cell2type])
   ylim([-7 2])
   set(gca, 'YDir','reverse')
    
end

tiledlayout(2,8)

nexttile; plot(tWave,spikeAvgMaxChanCell1(1:8,:)');   title('unit 79, \newline ch: 1-8')
nexttile; plot(tWave,spikeAvgMaxChanCell1(9:16,:)');  title('unit 79, \newline ch: 9-16')
nexttile; plot(tWave,spikeAvgMaxChanCell1(17:24,:)'); title('unit 79, \newline ch: 17-24')
nexttile; plot(tWave,spikeAvgMaxChanCell1(25:32,:)'); title('unit 79, \newline ch: 25-32')
nexttile; plot(tWave,spikeAvgMaxChanCell1(33:40,:)'); title('unit 79, \newline ch: 33-40')
nexttile; plot(tWave,spikeAvgMaxChanCell1(41:48,:)'); title('unit 79, \newline ch: 41-48')
nexttile; plot(tWave,spikeAvgMaxChanCell1(49:56,:)'); title('unit 79, \newline ch: 49-56')
nexttile; plot(tWave,spikeAvgMaxChanCell1(57:64,:)'); title('unit 79, \newline ch: 57-64')

nexttile; plot(tWave,spikeAvgMaxChanCell2(1:8,:)');   title('unit 118, \newline ch: 1-8')
nexttile; plot(tWave,spikeAvgMaxChanCell2(9:16,:)');  title('unit 118, \newline ch: 9-16')
nexttile; plot(tWave,spikeAvgMaxChanCell2(17:24,:)'); title('unit 118, \newline ch: 17-24')
nexttile; plot(tWave,spikeAvgMaxChanCell2(25:32,:)'); title('unit 118, \newline ch: 25-32')
nexttile; plot(tWave,spikeAvgMaxChanCell2(33:40,:)'); title('unit 118, \newline ch: 33-40')
nexttile; plot(tWave,spikeAvgMaxChanCell2(41:48,:)'); title('unit 118, \newline ch: 41-48')
nexttile; plot(tWave,spikeAvgMaxChanCell2(49:56,:)'); title('unit 118, \newline ch: 49-56')
nexttile; plot(tWave,spikeAvgMaxChanCell2(57:64,:)'); title('unit 118, \newline ch: 57-64')

%%

tiledlayout(2,8)

tWave = (-35:36)*(1/30);

nexttile; plot(tWave,spikeAvgMaxChanCell1(1:8,:)');   title('unit 20, \newline ch: 1-8')
nexttile; plot(tWave,spikeAvgMaxChanCell1(9:16,:)');  title('unit 20, \newline ch: 9-16')
nexttile; plot(tWave,spikeAvgMaxChanCell1(17:24,:)'); title('unit 20, \newline ch: 17-24')
nexttile; plot(tWave,spikeAvgMaxChanCell1(25:32,:)'); title('unit 20, \newline ch: 25-32')
nexttile; plot(tWave,spikeAvgMaxChanCell1(33:40,:)'); title('unit 20, \newline ch: 33-40')
nexttile; plot(tWave,spikeAvgMaxChanCell1(41:48,:)'); title('unit 20, \newline ch: 41-48')
nexttile; plot(tWave,spikeAvgMaxChanCell1(49:56,:)'); title('unit 20, \newline ch: 49-56')
nexttile; plot(tWave,spikeAvgMaxChanCell1(57:64,:)'); title('unit 20, \newline ch: 57-64')

nexttile; plot(tWave,spikeAvgMaxChanCell2(1:8,:)');   title('unit 45, \newline ch: 1-8')
nexttile; plot(tWave,spikeAvgMaxChanCell2(9:16,:)');  title('unit 45, \newline ch: 9-16')
nexttile; plot(tWave,spikeAvgMaxChanCell2(17:24,:)'); title('unit 45, \newline ch: 17-24')
nexttile; plot(tWave,spikeAvgMaxChanCell2(25:32,:)'); title('unit 45, \newline ch: 25-32')
nexttile; plot(tWave,spikeAvgMaxChanCell2(33:40,:)'); title('unit 45, \newline ch: 33-40')
nexttile; plot(tWave,spikeAvgMaxChanCell2(41:48,:)'); title('unit 45, \newline ch: 41-48')
nexttile; plot(tWave,spikeAvgMaxChanCell2(49:56,:)'); title('unit 45, \newline ch: 49-56')
nexttile; plot(tWave,spikeAvgMaxChanCell2(57:64,:)'); title('unit 45, \newline ch: 57-64')


