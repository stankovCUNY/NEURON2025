function [iiPairsAcrossShanks] = ii_conv_list(dayNo)

    fileDir = '/media/nasko/3ce6d8b6-2570-4c8f-941f-9b4b0641301c/home/nasko/CUNY_Work_Han_Kamran_Utku_dataV2/data/';

    %%
    
    if strcmp(dayNo,'1')
        day1 = load([fileDir 'AG_2019-12-23_NSD/' 'UtkuDay1_jscale1_alpha5_pairs.mat']);

        ExcPairs = day1.pairs.ExcPairs(:,1:2);
        InhPairs = day1.pairs.InhPairs(:,1:2);
        GapPairs = day1.pairs.GapPairs(:,1:2);

        allPairs = [ExcPairs; InhPairs; GapPairs];
        [~,idxUnique,~] = unique(allPairs,'rows');
        day1_iiPairs = allPairs(idxUnique,1:2);

        % screen and across shank and for ii pairs
        UtkuData = load([fileDir 'AG_2019-12-23_NSD/' 'AG_2019-12-23_NSD.cell_metrics.cellinfo.mat']);

        iiPairsAcrossShanks = [];

        for i = 1:length(day1_iiPairs)

            shank1    = UtkuData.cell_metrics.shankID(day1_iiPairs(i,1));
            shank2    = UtkuData.cell_metrics.shankID(day1_iiPairs(i,2));

            cell1type = UtkuData.cell_metrics.putativeCellType{1,day1_iiPairs(i,1)};
            cell2type = UtkuData.cell_metrics.putativeCellType{1,day1_iiPairs(i,2)};

            if strcmp(cell1type,'Pyramidal Cell')
                cell1type = 'p';
            elseif strcmp(cell1type,'Narrow Interneuron') || strcmp(cell1type,'Wide Interneuron')
                cell1type = 'i';
            end

            if strcmp(cell2type,'Pyramidal Cell')
                cell2type = 'p';
            elseif strcmp(cell2type,'Narrow Interneuron') || strcmp(cell2type,'Wide Interneuron')
                cell2type = 'i';
            end

            if (strcmp(cell1type,'i') && strcmp(cell2type,'i')) && (shank1 ~= shank2) 
                iiPairsAcrossShanks = [iiPairsAcrossShanks; day1_iiPairs(i,:)];
            end

        end

    %%
    
    elseif strcmp(dayNo,'2')
        day2 = load([fileDir 'AG_2019-12-27_NSD/' 'UtkuDay2_jscale1_alpha5_pairs.mat']);

        ExcPairs = day2.pairs.ExcPairs(:,1:2);
        InhPairs = day2.pairs.InhPairs(:,1:2);
        GapPairs = day2.pairs.GapPairs(:,1:2);

        allPairs = [ExcPairs; InhPairs; GapPairs];
        [~,idxUnique,~] = unique(allPairs,'rows');
        day2_iiPairs = allPairs(idxUnique,1:2);

        % screen and across shank and for ii pairs
        UtkuData = load([fileDir 'AG_2019-12-27_NSD/' 'AG_2019-12-27_NSD.cell_metrics.cellinfo.mat']);

        iiPairsAcrossShanks = [];

        for i = 1:length(day2_iiPairs)

            shank1    = UtkuData.cell_metrics.shankID(day2_iiPairs(i,1));
            shank2    = UtkuData.cell_metrics.shankID(day2_iiPairs(i,2));

            cell1type = UtkuData.cell_metrics.putativeCellType{1,day2_iiPairs(i,1)};
            cell2type = UtkuData.cell_metrics.putativeCellType{1,day2_iiPairs(i,2)};

            if strcmp(cell1type,'Pyramidal Cell')
                cell1type = 'p';
            elseif strcmp(cell1type,'Narrow Interneuron') || strcmp(cell1type,'Wide Interneuron')
                cell1type = 'i';
            end

            if strcmp(cell2type,'Pyramidal Cell')
                cell2type = 'p';
            elseif strcmp(cell2type,'Narrow Interneuron') || strcmp(cell2type,'Wide Interneuron')
                cell2type = 'i';
            end

            if (strcmp(cell1type,'i') && strcmp(cell2type,'i')) && (shank1 ~= shank2) 
                iiPairsAcrossShanks = [iiPairsAcrossShanks; day2_iiPairs(i,:)];
            end

        end
    end
        
end