RoySummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/RoySummaryStatsGJ.mat','fracGJ');
AGsummaryStats  = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/AGsummaryStatsGJ.mat','fracGJ','brainRegions','screenPairType');
NSUsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/NSUsummaryStatsGJ.mat','fracGJ');

BOTsummaryStats = load('/home/nasko/CUNY_Work_Han_Kamran_Nat/data/wake_new/BOTsummaryStatsGJ.mat','fracGJ','brainRegions','screenPairType');

brainRegions  = [ "CA1"   
                  "CA3" 
                  "APN"   
                  "CA1"   
                  "CA3"   
                  "DG"    
                  "Eth"   
                  "HPF"   
                  "IntG"  
                  "LGd"   
                  "LGv"   
                  "LP"    
                  "MB"    
                  "MGd"   
                  "MGm"   
                  "MGv"   
                  "MRN"   
                  "NOT"   
                  "PO"    
                  "POL"   
                  "ProS"  
                  "SCig"  
                  "SGN"   
                  "SUB"   
                  "TH"    
                  "VIS"   
                  "VISal" 
                  "VISam" 
                  "VISl"  
                  "VISli" 
                  "VISmma"
                  "VISmmp"
                  "VISp"  
                  "VISpm" 
                  "VISrl" 
                  "VL"    
                  "VPM"   
                  "undef. grey"  ];

hcomb = figure(102);

res_type = 'QHD';
pos = [70 230 1920*0.04 (2/3)*1080/5]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

b = bar(categorical(brainRegions(1:2)),100*[[RoySummaryStats.fracGJ{1,1} + AGsummaryStats.fracGJ{1,1}(1) + NSUsummaryStats.fracGJ{1,1}]/6 ...
                                     [RoySummaryStats.fracGJ{2,1} + AGsummaryStats.fracGJ{2,1}(1) + NSUsummaryStats.fracGJ{2,1}]/6 ...
                                     [RoySummaryStats.fracGJ{4,1} + AGsummaryStats.fracGJ{4,1}(1) + NSUsummaryStats.fracGJ{4,1}]/6; ...
                                     [0 +                            AGsummaryStats.fracGJ{1,1}(2) + 0]/2 ...
                                     [0 +                            AGsummaryStats.fracGJ{2,1}(2) + 0]/2 ...
                                     [0 +                            AGsummaryStats.fracGJ{4,1}(2) + 0]/2]);

ylabel('Percent')
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
box off
                               
close 102

hcomb = figure(102);

%     tiledlayout(1,7,'Padding','none','TileSpacing','compact')
%     nexttile([1,4])

res_type = 'QHD';
pos = [70 230 1920*0.3 (2/3)*1080/4]; a_offset = [0 850 100 900]'; b_offset = [0 0 -100 -100]';
arrayfun(@(a) set(a, 'Position', pos), hcomb(:));

maskNoPairs = find(~(((BOTsummaryStats.fracGJ{1,1} < 0.01) & (BOTsummaryStats.fracGJ{2,1} < 0.01)) & (BOTsummaryStats.fracGJ{3,1} < 0.01)));
[~,tempIdx] = sort(BOTsummaryStats.fracGJ{1,1}(maskNoPairs),'descend');
maskNoPairs = maskNoPairs(tempIdx);

brainRegionsMice = brainRegions(3:end);
brainRegionsMice = categorical(brainRegionsMice(maskNoPairs));
brainRegionsMice = reordercats(brainRegionsMice,cellstr(brainRegionsMice)');

bar(brainRegionsMice,100*[BOTsummaryStats.fracGJ{1,1}(maskNoPairs)/58 BOTsummaryStats.fracGJ{2,1}(maskNoPairs)/58 BOTsummaryStats.fracGJ{4,1}(maskNoPairs)/58]);
legend(AGsummaryStats.screenPairType{1}, AGsummaryStats.screenPairType{2}, ...
[AGsummaryStats.screenPairType{3} ' or ' AGsummaryStats.screenPairType{4}])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
box off

