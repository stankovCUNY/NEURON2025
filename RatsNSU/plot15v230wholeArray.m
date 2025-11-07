tiledlayout(2,12)

tWave = (-35:36)*(1/30);

nexttile; plot(tWave,nonArtifactAvg1(:,1:16));    title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 1-16'])
nexttile; plot(tWave,nonArtifactAvg1(:,17:32));   title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 17-31'])
nexttile; plot(tWave,nonArtifactAvg1(:,33:48));   title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 33-48'])
nexttile; plot(tWave,nonArtifactAvg1(:,49:64));   title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 49-64'])
nexttile; plot(tWave,nonArtifactAvg1(:,65:80));   title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 65-80'])
nexttile; plot(tWave,nonArtifactAvg1(:,81:96));   title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 81-96'])
nexttile; plot(tWave,nonArtifactAvg1(:,97:112));  title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 97-112'])
nexttile; plot(tWave,nonArtifactAvg1(:,113:128)); title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 113-128'])
nexttile; plot(tWave,nonArtifactAvg1(:,129:144)); title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 129-144'])
nexttile; plot(tWave,nonArtifactAvg1(:,145:160)); title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 145-160'])
nexttile; plot(tWave,nonArtifactAvg1(:,161:176)); title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 161-176'])
nexttile; plot(tWave,nonArtifactAvg1(:,177:192)); title(['unit ' num2str(cellID1) ' non-0lag, \newline ch: 177-192'])

nexttile; plot(tWave,artifactAvg1(:,1:16));    title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 1-16'])
nexttile; plot(tWave,artifactAvg1(:,17:32));   title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 17-31'])
nexttile; plot(tWave,artifactAvg1(:,33:48));   title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 33-48'])
nexttile; plot(tWave,artifactAvg1(:,49:64));   title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 49-64'])
nexttile; plot(tWave,artifactAvg1(:,65:80));   title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 65-80'])
nexttile; plot(tWave,artifactAvg1(:,81:96));   title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 81-96'])
nexttile; plot(tWave,artifactAvg1(:,97:112));  title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 97-112'])
nexttile; plot(tWave,artifactAvg1(:,113:128)); title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 113-128'])
nexttile; plot(tWave,artifactAvg1(:,129:144)); title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 129-144'])
nexttile; plot(tWave,artifactAvg1(:,145:160)); title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 145-160'])
nexttile; plot(tWave,artifactAvg1(:,161:176)); title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 161-176'])
nexttile; plot(tWave,artifactAvg1(:,177:192)); title(['unit ' num2str(cellID1) ' 0lag, \newline ch: 177-192'])

%% 

tiledlayout(2,12)

tWave = (-35:36)*(1/30);

nexttile; plot(tWave,nonArtifactAvg2(:,1:16));    title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 1-16'])
nexttile; plot(tWave,nonArtifactAvg2(:,17:32));   title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 17-31'])
nexttile; plot(tWave,nonArtifactAvg2(:,33:48));   title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 33-48'])
nexttile; plot(tWave,nonArtifactAvg2(:,49:64));   title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 49-64'])
nexttile; plot(tWave,nonArtifactAvg2(:,65:80));   title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 65-80'])
nexttile; plot(tWave,nonArtifactAvg2(:,81:96));   title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 81-96'])
nexttile; plot(tWave,nonArtifactAvg2(:,97:112));  title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 97-112'])
nexttile; plot(tWave,nonArtifactAvg2(:,113:128)); title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 113-128'])
nexttile; plot(tWave,nonArtifactAvg2(:,129:144)); title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 129-144'])
nexttile; plot(tWave,nonArtifactAvg2(:,145:160)); title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 145-160'])
nexttile; plot(tWave,nonArtifactAvg2(:,161:176)); title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 161-176'])
nexttile; plot(tWave,nonArtifactAvg2(:,177:192)); title(['unit ' num2str(cellID2) ' non-0lag, \newline ch: 177-192'])


nexttile; plot(tWave,artifactAvg2(:,1:16));    title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 1-16'])
nexttile; plot(tWave,artifactAvg2(:,17:32));   title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 17-31'])
nexttile; plot(tWave,artifactAvg2(:,33:48));   title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 33-48'])
nexttile; plot(tWave,artifactAvg2(:,49:64));   title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 49-64'])
nexttile; plot(tWave,artifactAvg2(:,65:80));   title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 65-80'])
nexttile; plot(tWave,artifactAvg2(:,81:96));   title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 81-96'])
nexttile; plot(tWave,artifactAvg2(:,97:112));  title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 97-112'])
nexttile; plot(tWave,artifactAvg2(:,113:128)); title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 113-128'])
nexttile; plot(tWave,artifactAvg2(:,129:144)); title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 129-144'])
nexttile; plot(tWave,artifactAvg2(:,145:160)); title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 145-160'])
nexttile; plot(tWave,artifactAvg2(:,161:176)); title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 161-176'])
nexttile; plot(tWave,artifactAvg2(:,177:192)); title(['unit ' num2str(cellID2) ' 0lag, \newline ch: 177-192'])



