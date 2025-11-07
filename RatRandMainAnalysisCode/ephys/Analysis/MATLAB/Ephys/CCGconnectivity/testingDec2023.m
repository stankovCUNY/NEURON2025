% load('RatRoyprocessedGroupStats.mat')
% 
% % plot(pairGroupStatTable.pairDistance(1,1)*ones(1,length(pairGroupStatTable.featuresCCGactualLat{1,1})),pairGroupStatTable.featuresCCGactualLat{1,1},'ko')
% % 
% % for loopPairs = 2:size(pairGroupStatTable,1)
% %     hold on
% %     plot(pairGroupStatTable.pairDistance(loopPairs,1)*ones(1,length(pairGroupStatTable.featuresCCGactualLat{loopPairs,1})),pairGroupStatTable.featuresCCGactualLat{loopPairs,1},'ko')
% %     hold off
% % end
% 
% tiledlayout(3,1)
% 
% %%
% nexttile
% 
% neuronIdx = 46;
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% 
% probe = buzsaki64probeLoc;
% 
% coordinates = [probe.x probe.y];
% 
% for loopChans = 1:length(probe.chanNo)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:length(probe.chanNo)
%     plot3((1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           5*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("ref activity map")
% set(gca,'FontSize',12)



% %% 
% nexttile
% waveform = pairGroupStatTable.tarWaveforms{neuronIdx,1};
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% 
% probe = buzsaki64probeLoc;
% 
% coordinates = [probe.x probe.y];
% 
% for loopChans = 1:length(probe.chanNo)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:length(probe.chanNo)
%     plot3((1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           5*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("tar activity map")
% set(gca,'FontSize',12)
% 
% % %%
% nexttile
% % find(pairGroupStatTable.refNeuronID == 79);
% 
% anaCCGidx  = (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 > -0.25) & (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 < 1);
% anaWaveidx = (pairGroupStatTable.tWave{1,1} > -0.25) & (pairGroupStatTable.tWave{1,1} < 1);
% 
% waveform = pairGroupStatTable.tarWaveforms{neuronIdx,1};
% 
% for loopChans = 1:length(probe.chanNo)
%     r_mat = corrcoef(pairGroupStatTable.pairRawCCG{neuronIdx,1}(anaCCGidx),waveform(loopChans,anaWaveidx));
%     r_sq(loopChans) = r_mat(1,2)^2;
% end
% 
% heat_values = r_sq;
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = 'R^2';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% set(gca,'FontSize',12)
% title("ref activity correlated with positive CCG lagtimes")

% %%
% 
% load('RatAGday1processedGroupStats.mat')
%  
% loadProbeInfo
% 
% tiledlayout(3,1)
% 
% nexttile
% 
% neuronIdx = 13;
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find(probe.ChannelNumberComingOutFromProbe == loopChans);
%     coordinates(loopChans,2) = probe.Z(a)';
%     coordinates(loopChans,1) = probe.X(a);
% 
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx);
% end
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           25*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("ref activity map")
% set(gca,'FontSize',12)
% 
% %% 
% 
% nexttile
% 
% waveform = pairGroupStatTable.tarWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find(probe.ChannelNumberComingOutFromProbe == loopChans);
%     coordinates(loopChans,2) = probe.Z(a)';
%     coordinates(loopChans,1) = probe.X(a);
% 
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx);
% end
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           25*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("ref activity map")
% set(gca,'FontSize',12)
% 
% %%
% 
% nexttile
% 
% anaCCGidx  = (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 > -0.25) & (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 < 1);
% anaWaveidx = (pairGroupStatTable.tWave{1,1} > -0.25) & (pairGroupStatTable.tWave{1,1} < 1);
% 
% waveform = pairGroupStatTable.tarWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     r_mat = corrcoef(pairGroupStatTable.pairRawCCG{neuronIdx,1}(anaCCGidx),waveform(loopChans,anaWaveidx));
%     r_sq(loopChans) = r_mat(1,2)^2;
% end
% 
% heat_values = r_sq;
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = 'R^2';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% set(gca,'FontSize',12)
% title("ref activity correlated with positive CCG lagtimes")
% 
% %%
% 
% load('RatNprocessedGroupStats.mat')
% 
% probegroupN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
% 
% loadProbeInfo
% 
% tiledlayout(3,1)
% 
% nexttile
% 
% neuronIdx = 13;
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find((cell2num(probegroupN.channel_id) + 1) == loopChans);
%     coordinates(loopChans,2) = probegroupN.y(a)';
%     coordinates(loopChans,1) = probegroupN.x(a);
% 
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx);
% end
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           25*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("ref activity map")
% set(gca,'FontSize',12)
% 
% %%
% 
% nexttile
% 
% waveform = pairGroupStatTable.tarWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find((cell2num(probegroupN.channel_id) + 1) == loopChans);
%     coordinates(loopChans,2) = probegroupN.y(a)';
%     coordinates(loopChans,1) = probegroupN.x(a);
% 
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx);
% end
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           25*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("tar activity map")
% set(gca,'FontSize',12)
% 
% %%
% 
% nexttile
% 
% anaCCGidx  = (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 > -0.25) & (pairGroupStatTable.CCGbinLagTimes{1,1}*1000 < 1);
% anaWaveidx = (pairGroupStatTable.tWave{1,1} > -0.25) & (pairGroupStatTable.tWave{1,1} < 1);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     r_mat = corrcoef(pairGroupStatTable.pairRawCCG{neuronIdx,1}(anaCCGidx),waveform(loopChans,anaWaveidx));
%     r_sq(loopChans) = r_mat(1,2)^2;
% end
% 
% heat_values = r_sq;
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = 'R^2';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% set(gca,'FontSize',12)
% title("ref activity correlated with positive CCG lagtimes")
% 
% %%
% 
% load('RatNprocessedGroupStats.mat')
% 
% probegroupN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
% 
% loadProbeInfo
% 
% tiledlayout(3,1)
% 
% nexttile
% 
% neuronIdx = 19;
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find((cell2num(probegroupN.channel_id) + 1) == loopChans);
%     coordinates(loopChans,2) = probegroupN.y(a)';
%     coordinates(loopChans,1) = probegroupN.x(a);
% 
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx);
% end
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx)';
% end
% 
% % Create a 3D surface plot
% 
% tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
% trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% %     title('Activity Plot');
% 
% axis equal
% 
% % xlim([0 35])
% % ylim([0 235])
% 
% % Add colorbar
% c = colorbar;
% c.Label.String = '[\muV]';
% colormap parula
% cmp = colormap;
% cmp = flipud(cmp);
% colormap(cmp);
% 
% % Adjust view for better visualization
% view(0, 90);
% 
% hold on
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%     plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
%           25*waveform(loopChans,:) + coordinates(loopChans,2), ...
%           max(waveform(:,tWaveZeroIdx))*ones(size(waveform,2),1),'k','LineWidth',1)
% end
% hold off
% title("ref activity map")
% set(gca,'FontSize',12)


%% Allen activity alignment maps

% clear
% 
% load('AllenInstitute786091066processedGroupStats.mat')
% 
% probe = npxProbeCoords;
% 
% neuronIdx = 254;
% 
% refID      = pairGroupStatTable.refNeuronID(neuronIdx,1);
% tarIdxList = find(pairGroupStatTable.refNeuronID == refID);
% tarChList  = pairGroupStatTable.tarChannel(tarIdxList,1);
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% tCCGZeroIdx  = find(pairGroupStatTable.CCGbinLagTimes{1,1} == 0);
% 
% tarExc = [pairGroupStatTable.GSPExc{tarIdxList,1}];
% tarInh = [pairGroupStatTable.GSPInh{tarIdxList,1}];
% 
% firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find(probe.chanIdx == loopChans);
%     if isempty(a)
%         continue
%     end
%     coordinates(loopChans,2) = probe.Y(a);
%     coordinates(loopChans,1) = probe.X(a);
% end
% 
% for loopTar = 1:length(tarIdxList)
%     a = find(probe.chanIdx == tarChList(loopTar));
%     tarCoordinates(loopTar,2) = probe.Y(a);
%     tarCoordinates(loopTar,1) = probe.X(a);
%     
%     if strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'p') || strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i-wide')
%         tarCellType{loopTar} = 'p';
%     elseif strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i-narrow')
%         tarCellType{loopTar} = 'i';
%     end
% end
% 
% tiledlayout(1,10,'Padding','none','TileSpacing','compact')
% 
% for loopTime = 0:2:18
%     nexttile
%     
% %     for loopTar = 1:length(tarIdxList)
% %         
% %         customColormap = [
% %         0.8, 0.2, 0.2;  % Red
% %         0.2, 0.8, 0.2;  % Green
% %         0.2, 0.2, 0.8   % Blue
% %         ];
% %        
% %         normalizedColorValues = interp1(linspace(min(firingProb(:,loopTar)), max(firingProb(:,loopTar)), size(customColormap, 1)), customColormap, firingProb(:,loopTar));
% %         
% %         scatter3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),1,firingProb(tCCGZeroIdx+loopTime,loopTar),normalizedColorValues(tCCGZeroIdx+loopTime,3),'^','LineWidth',5);
% % %         caxis([min(firingProb(:,loopTar)), max(firingProb(:,loopTar))]);
% %     end
%     
%     for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%         a = find(probe.chanIdx == loopChans);
%         if isempty(a)
%             continue
%         end
%         heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx + loopTime)';
%     end
% 
%     % Create a 3D surface plot
% 
%     tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
%     ytickangle( 90 ) 
%     %     title('Activity Plot');
% 
%     axis equal
% 
%     % xlim([0 35])
%     % ylim([0 235])
% 
%     % Add colorbar
% %     c = colorbar;
% %     c.Label.String = '[\muV]';
% %     colormap parula
% %     cmp = colormap;
% %     cmp = flipud(cmp);
% %     colormap(cmp);
% 
%     % Adjust view for better visualization
%     view(0, 90);
%     
% %     hold on
% %     
% %     for loopTar = 1:length(tarIdxList)
% %         if strcmp(tarCellType(loopTar),'p') 
% %             plot3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),100,'k^','LineWidth',2,'MarkerSize',1000*firingProb(tCCGZeroIdx+loopTime,loopTar));
% %         elseif strcmp(tarCellType(loopTar),'i')
% %             plot3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),100,'ko','LineWidth',2,'MarkerSize',1000*firingProb(tCCGZeroIdx+loopTime,loopTar));
% %         end
% %     end
% %     
% % %     plot3(tarCoordinates(strcmp(tarCellType,'p'),1),tarCoordinates(strcmp(tarCellType,'p'),2),ones(length(find(strcmp(tarCellType,'p'))),1),'k^','LineWidth',2);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'i'),1),tarCoordinates(strcmp(tarCellType,'i'),2),ones(length(find(strcmp(tarCellType,'i'))),1),'ko','LineWidth',2);
% %     
% % %     % excited
% % %     filt = tarExc(tCCGZeroIdx + loopTime,:);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'r^','LineWidth',2);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'ro','LineWidth',2);
% % %         
% % %     % inhited
% % %     filt = tarInh(tCCGZeroIdx + loopTime,:);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'c^','LineWidth',2);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'co','LineWidth',2);
% % %     
% % %     % neither
% % %     filt = ~tarExc(tCCGZeroIdx + loopTime,:) & ~tarInh(tCCGZeroIdx + loopTime,:);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'k^','LineWidth',2);
% % %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'ko','LineWidth',2);
% %     
% % %     for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% % %         plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
% % %               25*waveform(loopChans,:) + coordinates(loopChans,2), ...
% % %               max(waveform(:,idx))*ones(size(waveform,2),1),'k','LineWidth',1)
% % %     end
% %     hold off
%     
%     clim([-60 25])
%     
%     title(['ref act.' newline num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
% 
% end

%%

% clear
% 
% animalsList = [715093703
%                719161530
%                721123822
%                732592105
%                737581020
%                739448407
%                742951821
%                743475441
%                744228101
%                746083955
%                750332458
%                750749662
%                751348571
%                754312389
%                754829445
%                755434585
%                756029989
%                757216464
%                757970808
%                758798717
%                759883607
%                760345702
%                760693773
%                761418226
%                762120172
%                762602078
%                763673393
%                766640955
%                767871931
%                768515987
%                771160300
%                771990200
%                773418906
%                774875821
%                778240327
%                778998620
%                779839471
%                781842082
%                786091066
%                787025148
%                789848216
%                791319847
%                793224716
%                794812542
%                797828357
%                798911424
%                799864342
%                816200189
%                819186360
%                819701982
%                821695405
%                829720705
%                831882777
%                835479236
%                839068429
%                839557629
%                840012044
%                847657808];
% 
% heat_values_combo = [];
% coordinates_combo = [];
% 
% firingProb_combo     = [];
% tarCoordinates_combo = [];
% tarCellType_combo    = [];
% 
% probe = npxProbeCoords;
% 
% for loopAnimal = 1:length(animalsList)
%     
%     display(animalsList(loopAnimal))
%     
%     load(['AllenInstitute' num2str(animalsList(loopAnimal)) 'processedGroupStats.mat'])
%     
% %     pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.refCellExplorerType == 'p'),:);
%     [refList,refIdxList] = unique(pairGroupStatTable.refNeuronID);
% 
%     for loopRefUnits = 1:length(refList)
%         
%         neuronIdx = refIdxList(loopRefUnits);
%                
%         a = find(probe.chanIdx == pairGroupStatTable.refChannel(neuronIdx,1));
% 
%         coorRefAdj(loopRefUnits,2) = probe.Y(a);
%         coorRefAdj(loopRefUnits,1) = probe.X(a);
%         
%     end
% 
%     for loopRefUnits = 1:length(refList)
%         
%         clear firingProb    
%         clear tarCoordinates 
%         clear tarCellType
%         
%         neuronIdx = refIdxList(loopRefUnits);
% 
%         refID      = pairGroupStatTable.refNeuronID(neuronIdx,1);
%         tarIdxList = find(pairGroupStatTable.refNeuronID == refID);
%         tarChList  = pairGroupStatTable.tarChannel(tarIdxList,1);
%         refNorm    = length(pairGroupStatTable.refSpikeTimes{neuronIdx,1});
%         
%         tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
%         tCCGZeroIdx  = find(pairGroupStatTable.CCGbinLagTimes{1,1} == 0);
% 
%         tarExc = [pairGroupStatTable.GSPExc{tarIdxList,1}];
%         tarInh = [pairGroupStatTable.GSPInh{tarIdxList,1}];
% 
% %         firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
%         firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
%         
%         firingProb_combo     = [firingProb_combo firingProb];
%         
%         waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
%         for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%             a = find(probe.chanIdx == loopChans);
% 
%             coordinates(loopChans,2) = probe.Y(a) - coorRefAdj(loopRefUnits,2);
%             coordinates(loopChans,1) = probe.X(a) - coorRefAdj(loopRefUnits,1);
%         end
% 
%         coordinates_combo = [coordinates_combo; coordinates];
% 
%         for loopTar = 1:length(tarIdxList)
%                         
%             a =  find(probe.chanIdx == tarChList(loopTar));
%             
%             tarCoordinates(loopTar,2) = probe.Y(a) - coorRefAdj(loopRefUnits,2);
%             tarCoordinates(loopTar,1) = probe.X(a) - coorRefAdj(loopRefUnits,1);
%     
%             if strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'p') || strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i-wide')
%                 tarCellType{loopTar} = 'p';
%             elseif strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i-narrow')
%                 tarCellType{loopTar} = 'i';
%             end
%         end
%         
%         tarCoordinates_combo = [tarCoordinates_combo; tarCoordinates];
%         tarCellType_combo    = [tarCellType_combo    tarCellType];
%         
%         for loopTime = 0:2:18
%             for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%                 a = find(probe.chanIdx == loopChans);
%                 
%                 heat_values(loopTime/2+1,loopChans) = waveform(loopChans,tWaveZeroIdx + loopTime)';
%             end
%         end
% 
%         heat_values_combo = [heat_values_combo heat_values];
%     end
% end
% 
% tiledlayout(2,10)
% 
% for loopTime = 0:2:18
%     
%     %% activity ref 
%     nexttile
%     
%     % Create a 3D surface plot
%     tri = delaunay(coordinates_combo(:,1), coordinates_combo(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   coordinates_combo(:,1), coordinates_combo(:,2), heat_values_combo(loopTime/2+1,:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
% 
%     axis equal
% 
%     xlim([-50 50])
%     ylim([-200 200])
% 
%     % Add colorbar
%     c = colorbar;
%     c.Label.String = '[\muV]';
%     colormap parula
%     cmp = colormap;
%     cmp = flipud(cmp);
%     colormap(cmp);
%     clim([-10 10])
% 
%     % Adjust view for better visualization
%     view(0, 90);
% 
% %     hold on
% %     for loopTar = 1:length(tarCellType_combo)
% %         if strcmp(tarCellType_combo(loopTar),'p')
% %             plot3(tarCoordinates_combo(loopTar,1),tarCoordinates_combo(loopTar,2),1,'k^','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb_combo(tCCGZeroIdx+loopTime,loopTar));
% %         elseif strcmp(tarCellType_combo(loopTar),'i')
% %             plot3(tarCoordinates_combo(loopTar,1),tarCoordinates_combo(loopTar,2),1,'ko','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb_combo(tCCGZeroIdx+loopTime,loopTar));
% %         end
% %     end
% %     hold off
% 
%     title(['ref activity map ' num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
%     
%     %% synchrony 
%      nexttile
%     
%     % Create a 3D surface plot
%     tri = delaunay(tarCoordinates_combo(:,1), tarCoordinates_combo(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   tarCoordinates_combo(:,1), tarCoordinates_combo(:,2), firingProb_combo(tCCGZeroIdx+loopTime,:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
% 
%     axis equal
% 
%     xlim([-50 50])
%     ylim([-200 200])
% 
%     % Add colorbar
%     c = colorbar;
%     c.Label.String = 'prob.';
%     colormap parula
%     cmp = colormap;
%     cmp = flipud(cmp);
%     colormap(cmp);
%     clim([0 0.01])
% 
%     % Adjust view for better visualization
%     view(0, 90);
% 
%     title(['tar synchrony map ' num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
%     
% end

%% Bapun's rats
% 
% clear
% 
% % load('RatNprocessedGroupStats.mat')
% load('RatSprocessedGroupStats.mat')
% % load('RatUprocessedGroupStats.mat')
% 
% % probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
% probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');
% % probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
% 
% [refList,refIdxList] = unique(pairGroupStatTable.refNeuronID);
% 
% neuronIdx = 127; %refIdxList(6);
% 
% refID      = pairGroupStatTable.refNeuronID(neuronIdx,1);
% tarIdxList = find(pairGroupStatTable.refNeuronID == refID);
% tarChList  = pairGroupStatTable.tarChannel(tarIdxList,1);
% 
% tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
% tCCGZeroIdx  = find(pairGroupStatTable.CCGbinLagTimes{1,1} == 0);
% 
% tarExc = [pairGroupStatTable.GSPExc{tarIdxList,1}];
% tarInh = [pairGroupStatTable.GSPInh{tarIdxList,1}];
% 
% firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
% 
% waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
% for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%     a = find((cell2num(probegroup.channel_id) + 1) == loopChans);
%     if isempty(a)
%         continue
%     end
%     coordinates(loopChans,2) = probegroup.y(a);
%     coordinates(loopChans,1) = probegroup.x(a);
% end
% 
% for loopTar = 1:length(tarIdxList)
%     a =  find((cell2num(probegroup.channel_id) + 1) == tarChList(loopTar));
%     tarCoordinates(loopTar,2) = probegroup.y(a);
%     tarCoordinates(loopTar,1) = probegroup.x(a);
%     
%     if strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'p')
%         tarCellType{loopTar} = 'p';
%     elseif strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i')
%         tarCellType{loopTar} = 'i';
%     end
% end
%     
% tiledlayout(10,1)
% 
% for loopTime = 0:2:18
%     nexttile
%     
% %     for loopTar = 1:length(tarIdxList)
% %         
% %         customColormap = [
% %         0.8, 0.2, 0.2;  % Red
% %         0.2, 0.8, 0.2;  % Green
% %         0.2, 0.2, 0.8   % Blue
% %         ];
% %        
% %         normalizedColorValues = interp1(linspace(min(firingProb(:,loopTar)), max(firingProb(:,loopTar)), size(customColormap, 1)), customColormap, firingProb(:,loopTar));
% %         
% %         scatter3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),1,firingProb(tCCGZeroIdx+loopTime,loopTar),normalizedColorValues(tCCGZeroIdx+loopTime,3),'^','LineWidth',5);
% % %         caxis([min(firingProb(:,loopTar)), max(firingProb(:,loopTar))]);
% %     end
%     
%     for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%         a = find((cell2num(probegroup.channel_id) + 1) == loopChans);
%         if isempty(a)
%             continue
%         end
%         heat_values(loopChans) = waveform(loopChans,tWaveZeroIdx + loopTime)';
%     end
% 
%     % Create a 3D surface plot
% 
%     tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
%     %     title('Activity Plot');
% 
%     axis equal
% 
%     % xlim([0 35])
%     % ylim([0 235])
% 
%     % Add colorbar
%     c = colorbar;
%     c.Label.String = '[\muV]';
%     colormap parula
%     cmp = colormap;
%     cmp = flipud(cmp);
%     colormap(cmp);
% 
%     % Adjust view for better visualization
%     view(0, 90);
%     
%     hold on
%     
%     for loopTar = 1:length(tarIdxList)
%         if strcmp(tarCellType(loopTar),'p')
%             plot3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),1,'k^','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb(tCCGZeroIdx+loopTime,loopTar));
%         elseif strcmp(tarCellType(loopTar),'i')
%             plot3(tarCoordinates(loopTar,1),tarCoordinates(loopTar,2),1,'ko','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb(tCCGZeroIdx+loopTime,loopTar));
%         end
%     end
%     
% %     plot3(tarCoordinates(strcmp(tarCellType,'p'),1),tarCoordinates(strcmp(tarCellType,'p'),2),ones(length(find(strcmp(tarCellType,'p'))),1),'k^','LineWidth',2);
% %     plot3(tarCoordinates(strcmp(tarCellType,'i'),1),tarCoordinates(strcmp(tarCellType,'i'),2),ones(length(find(strcmp(tarCellType,'i'))),1),'ko','LineWidth',2);
%     
% %     % excited
% %     filt = tarExc(tCCGZeroIdx + loopTime,:);
% %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'r^','LineWidth',2);
% %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'ro','LineWidth',2);
% %         
% %     % inhited
% %     filt = tarInh(tCCGZeroIdx + loopTime,:);
% %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'c^','LineWidth',2);
% %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'co','LineWidth',2);
% %     
% %     % neither
% %     filt = ~tarExc(tCCGZeroIdx + loopTime,:) & ~tarInh(tCCGZeroIdx + loopTime,:);
% %     plot3(tarCoordinates(strcmp(tarCellType,'p') & filt,1),tarCoordinates(strcmp(tarCellType,'p') & filt,2),ones(length(find(strcmp(tarCellType,'p') & filt)),1),'k^','LineWidth',2);
% %     plot3(tarCoordinates(strcmp(tarCellType,'i') & filt,1),tarCoordinates(strcmp(tarCellType,'i') & filt,2),ones(length(find(strcmp(tarCellType,'i') & filt)),1),'ko','LineWidth',2);
%     
% %     for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% %         plot3(0.5*(1:size(waveform,2))  + coordinates(loopChans,1) - 25, ...
% %               25*waveform(loopChans,:) + coordinates(loopChans,2), ...
% %               max(waveform(:,idx))*ones(size(waveform,2),1),'k','LineWidth',1)
% %     end
%     hold off
%     
% %     clim([-0.6 0.05])
%     
%     title(['ref activity map ' num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
% 
% end

%%

% clear
% 
% heat_values_combo = [];
% coordinates_combo = [];
% 
% firingProb_combo     = [];
% tarCoordinates_combo = [];
% tarCellType_combo    = [];
% 
% for loopAnimal = 1:3
%     
%     if loopAnimal == 1
%         load('RatNprocessedGroupStats.mat')
%         probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
%     elseif loopAnimal == 2
%         load('RatSprocessedGroupStats.mat')
%         probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatS-Day2NSD-2020-11-27_10-22-29.probegroup.mat');
%     elseif loopAnimal == 3
%         load('RatUprocessedGroupStats.mat')
%         probegroup = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatU_Day2NSD_2021-07-24_08-16-38.probegroup.mat');
%     end
%     
%     pairGroupStatTable = pairGroupStatTable(cell2num(pairGroupStatTable.refCellExplorerType == 'p'),:);
%     [refList,refIdxList] = unique(pairGroupStatTable.refNeuronID);
% 
%     for loopRefUnits = 1:length(refList)
% 
%         neuronIdx = refIdxList(loopRefUnits);
% 
%         a = find((cell2num(probegroup.channel_id) + 1) == pairGroupStatTable.refChannel(neuronIdx,1));
% 
%         coorRefAdj(loopRefUnits,2) = probegroup.y(a);
%         coorRefAdj(loopRefUnits,1) = probegroup.x(a);
% 
%     end
% 
%     for loopRefUnits = 1:length(refList)
%         
%         clear firingProb    
%         clear tarCoordinates 
%         clear tarCellType
%         
%         neuronIdx = refIdxList(loopRefUnits);
% 
%         refID      = pairGroupStatTable.refNeuronID(neuronIdx,1);
%         tarIdxList = find(pairGroupStatTable.refNeuronID == refID);
%         tarChList  = pairGroupStatTable.tarChannel(tarIdxList,1);
%         refNorm    = length(pairGroupStatTable.refSpikeTimes{neuronIdx,1});
%         
%         tWaveZeroIdx = find(pairGroupStatTable.tWave{1,1} == 0);
%         tCCGZeroIdx  = find(pairGroupStatTable.CCGbinLagTimes{1,1} == 0);
% 
%         tarExc = [pairGroupStatTable.GSPExc{tarIdxList,1}];
%         tarInh = [pairGroupStatTable.GSPInh{tarIdxList,1}];
% 
% %         firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
%         firingProb = [pairGroupStatTable.pairRawCCG{tarIdxList,1}]./sum([pairGroupStatTable.pairRawCCG{tarIdxList,1}]);
%         
%         firingProb_combo     = [firingProb_combo firingProb];
%         
%         waveform = pairGroupStatTable.refWaveforms{neuronIdx,1};
% 
%         for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
% 
%             a = find((cell2num(probegroup.channel_id) + 1) == loopChans);
%             if isempty(a)
%                 continue
%             end
%             coordinates(loopChans,2) = probegroup.y(a) - coorRefAdj(loopRefUnits,2);
%             coordinates(loopChans,1) = probegroup.x(a) - coorRefAdj(loopRefUnits,1);
%         end
% 
%         coordinates_combo = [coordinates_combo; coordinates];
% 
%         for loopTar = 1:length(tarIdxList)
%             a =  find((cell2num(probegroup.channel_id) + 1) == tarChList(loopTar));
%             tarCoordinates(loopTar,2) = probegroup.y(a) - coorRefAdj(loopRefUnits,2);
%             tarCoordinates(loopTar,1) = probegroup.x(a) - coorRefAdj(loopRefUnits,1);
%     
%             if strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'p')
%                 tarCellType{loopTar} = 'p';
%             elseif strcmp(pairGroupStatTable.tarCellExplorerType{loopTar,1},'i')
%                 tarCellType{loopTar} = 'i';
%             end
%         end
%         
%         tarCoordinates_combo = [tarCoordinates_combo; tarCoordinates];
%         tarCellType_combo    = [tarCellType_combo    tarCellType];
%         
%         for loopTime = 0:2:18
%             for loopChans = 1:size(pairGroupStatTable.refWaveforms{1,1},1)
%                 a = find((cell2num(probegroup.channel_id) + 1) == loopChans);
%                 if isempty(a)
%                     continue
%                 end
%                 heat_values(loopTime/2+1,loopChans) = waveform(loopChans,tWaveZeroIdx + loopTime)';
%             end
%         end
% 
%         heat_values_combo = [heat_values_combo heat_values];
%     end
% end
% 
% tiledlayout(10,2)
% 
% for loopTime = 0:2:18
%     
%     %% activity ref 
%     nexttile
%     
%     % Create a 3D surface plot
%     tri = delaunay(coordinates_combo(:,1), coordinates_combo(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   coordinates_combo(:,1), coordinates_combo(:,2), heat_values_combo(loopTime/2+1,:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
% 
%     axis equal
% 
%     xlim([-500 500])
%     ylim([-200 200])
% 
%     % Add colorbar
%     c = colorbar;
%     c.Label.String = '[\muV]';
%     colormap parula
%     cmp = colormap;
%     cmp = flipud(cmp);
%     colormap(cmp);
%     clim([-0.05 0.05])
% 
%     % Adjust view for better visualization
%     view(0, 90);
% 
% %     hold on
% %     for loopTar = 1:length(tarCellType_combo)
% %         if strcmp(tarCellType_combo(loopTar),'p')
% %             plot3(tarCoordinates_combo(loopTar,1),tarCoordinates_combo(loopTar,2),1,'k^','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb_combo(tCCGZeroIdx+loopTime,loopTar));
% %         elseif strcmp(tarCellType_combo(loopTar),'i')
% %             plot3(tarCoordinates_combo(loopTar,1),tarCoordinates_combo(loopTar,2),1,'ko','LineWidth',2,'MarkerSize',0.0001 + 1000*firingProb_combo(tCCGZeroIdx+loopTime,loopTar));
% %         end
% %     end
% %     hold off
% 
%     title(['ref activity map ' num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
%     
%     %% synchrony 
%      nexttile
%     
%     % Create a 3D surface plot
%     tri = delaunay(tarCoordinates_combo(:,1), tarCoordinates_combo(:,2)); % Create Delaunay triangulation
%     trisurf(tri,   tarCoordinates_combo(:,1), tarCoordinates_combo(:,2), firingProb_combo(tCCGZeroIdx+loopTime,:), 'FaceColor', 'interp', 'EdgeColor', 'none');
%     xlabel('x [\mum]');
%     ylabel('y [\mum]');
% 
%     axis equal
% 
%     xlim([-500 500])
%     ylim([-200 200])
% 
%     % Add colorbar
%     c = colorbar;
%     c.Label.String = 'prob.';
%     colormap parula
%     cmp = colormap;
%     cmp = flipud(cmp);
%     colormap(cmp);
%     clim([0 0.01])
% 
%     % Adjust view for better visualization
%     view(0, 90);
% 
%     title(['tar synchrony map ' num2str(loopTime/30) ' ms'])
%     set(gca,'FontSize',12)
%     
% end


%% 

clear

fs             = 30000;
binSize        = 1/fs;

duration       = 0.007;
lagLimit       = duration/2;

loopAnimal = 4;

pairIdx = 10;

if loopAnimal == 1
    load('RatNprocessedGroupStats.mat')
    for loopShanks = 1:8
        allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/spike_times.npy']));
        spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/spike_templates.npy']));
        templateChannels{loopShanks} = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/pc_feature_ind.npy']));
        
        probegroupN = load('/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatN_Day2_2019-10-11_03-58-54.probegroup.mat');
        
%         channel_map{loopShanks} = double(cell2num(probegroupN.channel_id(probegroupN.shank_id == (loopShanks - 1))));
        
        channel_map{loopShanks}      = double(readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(loopShanks) '/channel_map.npy']));
        
        % adjust channel_map
        channel_map_adjusted{1} = channel_map{1};
        if loopShanks > 1
            adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
            channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust);
        end
    end
elseif loopAnimal == 2
    load('RatSprocessedGroupStats.mat')
    allSpikeLat      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/spike_times.npy'));
    spikeTemplates   = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/spike_templates.npy'));
    templateChannels = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/pc_feature_ind.npy'));
    channel_map      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/channel_map.npy'));
elseif loopAnimal == 3
    load('RatUprocessedGroupStats.mat')
    allSpikeLat      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/spike_times.npy'));
    spikeTemplates   = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/spike_templates.npy'));
    templateChannels = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/pc_feature_ind.npy'));
    channel_map      = double(readNPY('/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/channel_map.npy'));
elseif loopAnimal == 4
    load('RatAGday1processedGroupStats.mat')
    for loopShanks = 1:6
        allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/spike_times.npy']));
        spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/spike_templates.npy']));
        templateChannels{loopShanks} = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/pc_feature_ind.npy']));
        channel_map{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-23_NSD/Units/s' num2str(loopShanks) '/channel_map.npy']));
        
        % adjust channel_map
        channel_map_adjusted{1} = channel_map{1};
        if loopShanks > 1
            adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
            channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust);
        end
    end
elseif loopAnimal == 5
    load('RatAGday2processedGroupStats.mat')
    for loopShanks = 1:5
        allSpikeLat{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/spike_times.npy']));
        spikeTemplates{loopShanks}   = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/spike_templates.npy']));
        templateChannels{loopShanks} = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/pc_feature_ind.npy']));
        channel_map{loopShanks}      = double(readNPY(['/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/AG_2019-12-27_NSD/Units/s' num2str(loopShanks) '/channel_map.npy']));
        
        % adjust channel_map
        channel_map_adjusted{1} = channel_map{1};
        if loopShanks > 1
            adjust(loopShanks-1) = length(channel_map{loopShanks-1}); 
            channel_map_adjusted{loopShanks} = channel_map{loopShanks} + sum(adjust);
        end
    end
end

R      = pairGroupStatTable.refSpikeTimes{pairIdx,1};
T      = pairGroupStatTable.tarSpikeTimes{pairIdx,1};

refLat = round(R*fs)';
tarLat = round(T*fs)';

tarChannel = pairGroupStatTable.tarChannel(pairIdx,1);

GSPExc = pairGroupStatTable.GSPExc{pairIdx,1};

[CCG,spikeTimesInBin] = CCGtrackedSpikeTimes(R,T,binSize,lagLimit);

Rsync = [];
Tsync = [];
for k = 1:211
    if (GSPExc(k) == 1) 
        if isempty(spikeTimesInBin{k})
            continue
        else
            Rsync = [Rsync; spikeTimesInBin{k}(:,1)];
            Tsync = [Tsync; spikeTimesInBin{k}(:,2)];
        end
    end
end

tarSynchLat = round(Tsync/(1/3e4))';

if sum(loopAnimal == [1,4,5]) == 1
    tarShank = pairGroupStatTable.tarShank(pairIdx);
    
    [~,idx,     ~] = intersect(allSpikeLat{tarShank},tarLat);
    [~,idxSynch,~] = intersect(allSpikeLat{tarShank},tarSynchLat);
    
    templateList      = unique(spikeTemplates{tarShank}(idx));
    templateSynchList = unique(spikeTemplates{tarShank}(idxSynch));
    
     % templateChannels(templateList+1,:);

    tiledlayout(2,1)
    nexttile
    histogram(spikeTemplates{tarShank}(idx),300)
    nexttile
    histogram(spikeTemplates{tarShank}(idxSynch),300)
    xlim([0 300])

    chanIdx = mode(spikeTemplates{tarShank}(idx));
%     mode(spikeTemplates{tarShank}(idxSynch));
    
    clusIdx = find(find((tarChannel - 1) == channel_map_adjusted{tarShank}) == (templateChannels{tarShank}(chanIdx + 1 ,:) ));

else
    [~,idx,     ~] = intersect(allSpikeLat,tarLat);
    [~,idxSynch,~] = intersect(allSpikeLat,tarSynchLat);

    templateList      = unique(spikeTemplates(idx));
    templateSynchList = unique(spikeTemplates(idxSynch));

    % templateChannels(templateList+1,:);

    tiledlayout(2,1)
    nexttile
    histogram(spikeTemplates(idx),300)
    nexttile
    histogram(spikeTemplates(idxSynch),300)
    xlim([0 300])

    chanIdx = mode(spikeTemplates(idx));
    mode(spikeTemplates(idxSynch));

    clusIdx = find(find((tarChannel - 1) == channel_map) == (templateChannels(chanIdx + 1 ,:) ));

end
    
% readNPY("/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/pc_features.npy");

try
    if loopAnimal == 1
        clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatN/RatNDay2NSD/shank' num2str(tarShank) '/pc_features_PCchan' num2str(clusIdx - 1) '.npy']);
    elseif loopAnimal == 2
        clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatS/RatSDay2NSD/RatS-Day2NSD-2020-11-27_10-22-29-1.GUI/pc_features_PCchan' num2str(clusIdx - 1) '.npy']);
    elseif loopAnimal == 3
        clusterVal = readNPY(['/media/nasko/WD_BLACK/BapunRawMultiArray/RatU/RatUDay2NSD/RatU_Day2NSD_2021-07-24_08-16-38-1.GUI/pc_features_PCchan' num2str(clusIdx - 1) '.npy']);
    elseif loopAnimal == 4
        clusterVal = readNPY(['/media/nasko/WD_BLACK/UtkuRawMultiArrayUnpackedNPY/AG_2019-12-23_NSD/Units/s' num2str(tarShank) '/pc_features_PCchan' num2str(clusIdx - 1) '.npy']);
    elseif loopAnimal == 5
        clusterVal = readNPY(['/media/nasko/WD_BLACK/UtkuRawMultiArrayUnpackedNPY/AG_2019-12-27_NSD/Units/s' num2str(tarShank) '/pc_features_PCchan' num2str(clusIdx - 1) '.npy']);
    end

    clusterValRef      = clusterVal(idx,1:3);
    clusterValRefSynch = clusterVal(idxSynch,1:3);

    tiledlayout(1,4)

    tR = (-lagLimit:binSize:lagLimit)*1000;

    nexttile
    plot(tR,CCG/length(R),'k')
    xlim([-1 1])
    ylims = get(gca,'ylim');
    ylabel('Spike Probability')
    xlabel('[ms]')
    %     title([refCelltype '-' tarCelltype ' pair ' num2str(d) '\mum ' region])
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')

    hold on
    scatter(tR(GSPExc == 1),CCG(GSPExc == 1)/length(R),100,'ro','LineWidth',1)
    % scatter(tR(114)*1e3,ccgR(114,1,2)/length(ref),100,'bo','LineWidth',1)
    hold off
    box off

    nexttile
    scatter(clusterValRef(:,1),clusterValRef(:,2))
    hold on
    scatter(clusterValRefSynch(:,1),clusterValRefSynch(:,2))
    hold off

    nexttile
    scatter(clusterValRef(:,1),clusterValRef(:,3))
    hold on
    scatter(clusterValRefSynch(:,1),clusterValRefSynch(:,3))
    hold off

    nexttile
    scatter(clusterValRef(:,2),clusterValRef(:,3))
    hold on
    scatter(clusterValRefSynch(:,2),clusterValRefSynch(:,3))
    hold off
    
catch
        
    display("No Spyking Circus PCA done for this target neuron")
        
end
%%

figure()

if loopAnimal == 1 % rat N
   chanData = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatN/ch' num2str(tarChannel) 'highpass300hz.mat'],'data'); 
elseif loopAnimal == 2 % rat S
   chanData = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatS/ch' num2str(tarChannel) 'highpass300hz.mat'],'data');
elseif loopAnimal == 3 % rat U
   chanData = load(['/media/nasko/WD_BLACK31/BapunFilteredDataPerChannel/RatU/ch' num2str(tarChannel) 'highpass300hz.mat'],'data');
elseif loopAnimal == 4 % rat AG day 1
   chanData = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_23_NSD/shank' num2str(pairGroupStatTable.tarShank(pairIdx)) '/ch' num2str(tarChannel) 'highpass300hz.mat'],'data');
elseif loopAnimal == 5 % rat AG day 1
   chanData = load(['/media/nasko/WD_BLACK31/UtkuFilteredDataPerChannel/AG/AG_2019-12_27_NSD/shank' num2str(pairGroupStatTable.tarShank(pairIdx)) '/ch' num2str(tarChannel) 'highpass300hz.mat'],'data');
end

% chanData       = double(chanData.data);
chanData       = chanData.data;

filterFlag     = false; 
noDemeanFlag   = false;
fpass          = 300;
preLength      = 5*30;
postLength     = 5*30;

tarLatRandPerm = sort(tarLat(randperm(length(tarLat),10000)));

[tarSponWaveMean, tarSponWaveforms]  = waveformAvg(chanData,tarLatRandPerm,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
[tarSynchWaveMean,tarSynchWaveforms] = waveformAvg(chanData,tarSynchLat,   preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);

patchline((-preLength:postLength-1)/(fs/1000), ... 
                 tarSponWaveforms(:,1), ...
                 'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
hold on
for i = 2:size(tarSponWaveforms,2)
    patchline((-preLength:postLength-1)/(fs/1000), ... 
             tarSponWaveforms(:,i), ...
             'linewidth',1,'edgealpha',0.25,'edgecolor',[.7 .7 .7])
end
for i = 1:size(tarSynchWaveforms,2)
    patchline((-preLength:postLength-1)/(fs/1000), ... 
             tarSynchWaveforms(:,i), ...
             'linewidth',1,'edgealpha',0.1,'edgecolor','r')
end
hold off
ylim([-6 2])
xlim([-1,1])
set(gca,'FontSize',5)
set(gca,'FontName','Arial')
xlabel('[ms]');
ylabel('[mV]');
set(gca, 'YDir','reverse')
box off










