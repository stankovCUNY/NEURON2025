% just inhibitory one-way

N = 100000;

Rthr     = 1 - 0.05;
Tthr     = 1 - 0.05;
TthrBias = 1 - 0.02; % inhibition bias selection

InhBiasSelTime = 10;

Rbin = zeros(N,1);
Tbin = zeros(N,1);

for i = 1:N
   
    flipR = rand(1);
    flipT = rand(1);
    
    if flipR > Rthr
        Rbin(i) = 1;
    end
    
    if (i > InhBiasSelTime) 
        if (sum(Rbin(i-InhBiasSelTime:i-1)) > 0)
            if flipT > TthrBias
                Tbin(i) = 1;
                
            end
        else
            if flipT > Tthr
                Tbin(i) = 1;
            end
        end
    else
        if flipT > Tthr
            Tbin(i) = 1;
        end
    end
    
    
end

R = find(Rbin)/10000;
T = find(Tbin)/10000;

tsOffsets = crosscorrelogram(R,T,[-0.005 0.005]);
subplot(1,2,1)
histogram(tsOffsets*1000, 100);

%% just inhibitory both ways

% N = 100000;
% 
% Rthr     = 1 - 0.10;
% Tthr     = 1 - 0.10;
% RthrBias = 1 - 0.01; % inhibition bias selection
% TthrBias = 1 - 0.01; % inhibition bias selection
% 
% InhBiasSelTime = 10;
% 
% Rbin = zeros(N,1);
% Tbin = zeros(N,1);
% 
% for i = 1:N
%    
%     flipR = rand(1);
%     flipT = rand(1);
%     
%     if (i > InhBiasSelTime) 
%         if (sum(Tbin(i-InhBiasSelTime:i-4)) > 0)
%             if flipR > RthrBias
%                 Rbin(i) = 1;
%                 
%             end
%         else
%             if flipR > Rthr
%                 Rbin(i) = 1;
%             end
%         end
%     else
%         if flipR > Rthr
%             Rbin(i) = 1;
%         end
%     end
%     
%     if (i > InhBiasSelTime) 
%         if (sum(Rbin(i-InhBiasSelTime:i-4)) > 0)
%             if flipT > TthrBias
%                 Tbin(i) = 1;
%                 
%             end
%         else
%             if flipT > Tthr
%                 Tbin(i) = 1;
%             end
%         end
%     else
%         if flipT > Tthr
%             Tbin(i) = 1;
%         end
%     end
%     
%     
% end
% 
% R = find(Rbin)/10000;
% T = find(Tbin)/10000;
% 
% tsOffsets = crosscorrelogram(R,T,[-0.005 0.005]);
% subplot(1,3,2)
% histogram(tsOffsets*1000, 100);

%%

N = 100000;

Rthr     = 1 - 0.05;
Tthr     = 1 - 0.05;

InhBiasSelTime = 10;

win = gausswin(InhBiasSelTime-1); % create Gaussian kernel
win = 1 - win;

RthrBias = 1 - 0.02*win; % inhibition bias selection
TthrBias = 1 - 0.02*win; % inhibition bias selection

Rbin = zeros(N,1);
Tbin = zeros(N,1);

for i = 1:N
   
    flipR = rand(1);
    flipT = rand(1);
    
    if (i > InhBiasSelTime) 
        if (sum(Tbin(i-InhBiasSelTime:i-4)) > 0)
            idx = find(Tbin(i-InhBiasSelTime:i-4));
            if flipR > RthrBias(idx(1))
                Rbin(i) = 1;
                
            end
        else
            if flipR > Rthr
                Rbin(i) = 1;
            end
        end
    else
        if flipR > Rthr
            Rbin(i) = 1;
        end
    end
    
    if (i > InhBiasSelTime) 
        if (sum(Rbin(i-InhBiasSelTime:i-4)) > 0)
            idx = find(Rbin(i-InhBiasSelTime:i-4));
            if flipT > TthrBias(idx(1))
                Tbin(i) = 1;
                
            end
        else
            if flipT > Tthr
                Tbin(i) = 1;
            end
        end
    else
        if flipT > Tthr
            Tbin(i) = 1;
        end
    end
    
    
end

R = find(Rbin)/10000;
T = find(Tbin)/10000;

tsOffsets = crosscorrelogram(R,T,[-0.005 0.005]);
subplot(1,2,2)
histogram(tsOffsets*1000, 100);