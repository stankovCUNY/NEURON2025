% dayNo   = {'1','2'};
% shankNo = {'1','2','3','4','5','6'};
% 
% for i = 1:length(dayNo)
%     for j = 1:length(shankNo)
%         if strcmp(dayNo{i},'2') && strcmp(shankNo{j},'6')
%             continue
%         end
%         UtkuCurateNoiseClusters(dayNo{i},shankNo{j});
%     end
% end

% UtkuScreenConvV2('1')
% UtkuScreenConvV2('2')
% UtkuScreenConvV2('B1')
% 
% nullFlag = true;
% UtkuRunCCGsV2('1',nullFlag)
% UtkuRunCCGsV2('2',nullFlag)
% UtkuRunCCGsV2('B1')

% nullFlag = true;
% UtkuRunCCGsExquisiteAndSymmetrical('1',nullFlag)
% UtkuRunCCGsExquisiteAndSymmetrical('2',nullFlag)

% nullFlag = true;
% UtkuRunCCGsExquisiteAndSymmetricalVer2('1',nullFlag)
% UtkuRunCCGsExquisiteAndSymmetricalVer2('2',nullFlag)

% UtkuRunCCGsIIV2ii('1')
% UtkuRunCCGsIIV2ii('2')
% 
% UtkuRunCCGsII('1')
% UtkuRunCCGsII('2')
% 
% UtkuWaveformSnippets('1')
% UtkuWaveformSnippets('2')

nullFlag = true;
UtkuRunCCGputativeGJ('1',nullFlag)
UtkuRunCCGputativeGJ('2',nullFlag)