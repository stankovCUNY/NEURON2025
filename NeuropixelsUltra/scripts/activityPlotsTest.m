neuronIdx = 6;

load('/home/nasko/CUNY_Work_NPUltraWaveforms/data/groupStatsSteinmetzDataset.mat')
waveform = pairGroupStatTable.refWaveforms{neuronIdx,1} + pairGroupStatTable.tarWaveforms{neuronIdx,1};
% % Parameters
% frame_time = 50; % Specify the time point for the frame in milliseconds
% 
% % Create a figure
% figure;
% 
% % Find the index corresponding to the specified time
% frame_index = find(frame_time >= (1:size(waveform_data, 1)), 1, 'last');
% 
% % Check if the specified time is within the valid range
% if isempty(frame_index)
%     error('Invalid frame_time. It should be within the range of the data.');
% end
% 
% % Create a heatmap
% heatmap(waveform_data(frame_index, :)', 'Colormap', 'viridis');
% 
% % Customize the heatmap as needed
% title(['Neuropixels Heatmap at Time = ' num2str(frame_time) ' ms']);
% xlabel('Channels');
% ylabel('Amplitude');
% 
% % Display colorbar
% colorbar;
% 
% % Pause to keep the figure open
% pause;

idx = 43;

chanLayout = flipud(reshape(1:384,[8,48])');

x = 0:5:5*7;
y = 0:5:5*47;

for loopChans = 1:(48*8)
   
    [a,b] = find(chanLayout == loopChans);
    coordinates(loopChans,2) = y(a);
    coordinates(loopChans,1) = x(b);
    
    heat_values(loopChans) = waveform(idx,loopChans);
end

% Create a 3D surface plot
figure;
tri = delaunay(coordinates(:,1), coordinates(:,2)); % Create Delaunay triangulation
trisurf(tri,   coordinates(:,1), coordinates(:,2), heat_values, 'FaceColor', 'interp', 'EdgeColor', 'none');
xlabel('x [\mum]');
ylabel('y [\mum]');
title('Activity Plot');

axis equal

% xlim([0 35])
% ylim([0 235])

% Add colorbar
c = colorbar;
c.Label.String = '[\muV]';
colormap parula
cmp = colormap;
cmp = flipud(cmp);
colormap(cmp);

% Adjust view for better visualization
view(0, 90);

hold on
for loopChans = 1:(48*8)
    plot3(0.05*(1:size(waveform,1))  + coordinates(loopChans,1) - 2.5, ...
          0.02*waveform(:,loopChans) + coordinates(loopChans,2), ...
          max(waveform(idx,:))*ones(size(waveform,1),1),'k','LineWidth',1)
end
hold off