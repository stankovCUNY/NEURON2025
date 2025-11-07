load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RatRoyprocessedGroupStats.mat')
 
royAlignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

RoyCCGLatencies  = [];
RoyWaveLatencies = [];

for loopPairs = 1:size(royAlignmentTable,1)
    RoyCCGLatencies  = [RoyCCGLatencies;   royAlignmentTable.CCGAlignedLatencies{loopPairs}];
    RoyWaveLatencies = [RoyWaveLatencies;  royAlignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

AG1CCGLatencies  = [];
AG1WaveLatencies = [];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/RatAGday1processedGroupStats.mat')

AG1AlignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

for loopPairs = 1:size(AG1AlignmentTable,1)
    AG1CCGLatencies  = [AG1CCGLatencies;   AG1AlignmentTable.CCGAlignedLatencies{loopPairs}];
    AG1WaveLatencies = [AG1WaveLatencies;  AG1AlignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

AG2CCGLatencies  = [];
AG2WaveLatencies = [];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/RatAGday2processedGroupStats.mat')

AG2AlignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

for loopPairs = 1:size(AG2AlignmentTable,1)
    AG2CCGLatencies  = [AG2CCGLatencies;   AG2AlignmentTable.CCGAlignedLatencies{loopPairs}];
    AG2WaveLatencies = [AG2WaveLatencies;  AG2AlignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

ratSCCGLatencies  = [];
ratSwaveLatencies = [];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatS/RatSprocessedGroupStats.mat')

ratSalignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

for loopPairs = 1:size(ratSalignmentTable,1)
    ratSCCGLatencies  = [ratSCCGLatencies;   ratSalignmentTable.CCGAlignedLatencies{loopPairs}];
    ratSwaveLatencies = [ratSwaveLatencies;  ratSalignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

ratUCCGLatencies  = [];
ratUwaveLatencies = [];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatU/RatUprocessedGroupStats.mat')

ratUalignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

for loopPairs = 1:size(ratUalignmentTable,1)
    ratUCCGLatencies  = [ratUCCGLatencies;   ratUalignmentTable.CCGAlignedLatencies{loopPairs}];
    ratUwaveLatencies = [ratUwaveLatencies;  ratUalignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

ratNCCGLatencies  = [];
ratNwaveLatencies = [];

load('/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/RatN/RatNprocessedGroupStats.mat')

ratNalignmentTable = pairGroupStatTable(pairGroupStatTable.flagAligned == 1,:);

for loopPairs = 1:size(ratNalignmentTable,1)
    ratNCCGLatencies  = [ratNCCGLatencies;   ratNalignmentTable.CCGAlignedLatencies{loopPairs}];
    ratNwaveLatencies = [ratNwaveLatencies;  ratNalignmentTable.WaveAlignedLatencies{loopPairs}']; 
end

%%

plot(RoyCCGLatencies,RoyWaveLatencies,  'o','LineWidth',3)
hold on
plot([AG1CCGLatencies; AG2CCGLatencies],[AG1WaveLatencies; AG2WaveLatencies],  'o','LineWidth',3)
plot(ratSCCGLatencies,ratSwaveLatencies,'o','LineWidth',3)
plot(ratUCCGLatencies,ratUwaveLatencies,'o','LineWidth',3)
plot(ratNCCGLatencies,ratNwaveLatencies,'o','LineWidth',3)
hold off

legend('Roy','AG','Rat S','Rat U','Rat N','Location','NorthWest')
ylabel('CCG feature latency [\musec]')
xlabel('Waveform feature latency [\musec]')
set(gca,'FontSize',12)
axis equal

%%

histogram([royAlignmentTable.pairDistance
AG1AlignmentTable.pairDistance
AG2AlignmentTable.pairDistance
ratSalignmentTable.pairDistance
ratUalignmentTable.pairDistance
ratNalignmentTable.pairDistance],10);

ylabel('pair count')
xlabel('pair distance [\mum]')
set(gca,'FontSize',12)

