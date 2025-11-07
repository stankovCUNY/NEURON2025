lgmx = 5;

ref=stimTimes(1:1e4)*1000;
tar=spikeTimesTar(1:1e4)*1000;

% ref=res1(1:1e4)*1000;
% tar=res2(1:1e4)*1000;
% 
% ref=res2sync*1000;
% tar=res1sync*1000;

xR = []; yR = []; xT = []; yT = [];  xS = []; yS = [];rSx = []; rSy = [];
for k = 1:numel(ref)
    rtem = ref(ref>=ref(k)-lgmx & ref<=ref(k)+lgmx) - ref(k);
    ttem = tar(tar>=ref(k)-lgmx & tar<=ref(k)+lgmx) - ref(k);
    xR = [xR;rtem];
    yR = [yR;k.*ones(size(rtem))];
    xT = [xT;ttem];
    yT = [yT;k.*ones(size(ttem))];

end

figure

subplot(3,1,1)
scatter(xT,yT,1,'filled','markerfacecolor','r'); hold on;
xlim([-lgmx,lgmx])
title('RASTER: Target spikes aligned with reference spikes')

subplot(3,1,2)
scatter(xR,yR,1,'filled','markerfacecolor','b'); hold on;
xlim([-lgmx,lgmx])
title('RASTER: Ref spikes aligned with ref spikes')

subplot(3,1,3)
scatter(xT,yT,1,'filled','markerfacecolor','r'); hold on;
hold on
scatter(xR,yR,1,'filled','markerfacecolor','b'); hold on;
hold off
xlim([-lgmx,lgmx])