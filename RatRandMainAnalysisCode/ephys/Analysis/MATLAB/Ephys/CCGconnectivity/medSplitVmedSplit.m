metric = 'spikeWidth';
load(['medSplitSets_' metric '_allStates.mat'])

figure(102)

jscale         = 1;
alpha_name     = 5;
duration       = 0.007;
fs             = 30000;
fpass          = 300;
binSize        = 1/fs;
fig_use        = 102;
njitter        = 500;
alpha          = 0.05;
for_grant      = false;
filterFlag     = false;

tiledlayout(2,1)

i = 15;

% for i = 1:size(medSplit.pairs,1) 
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(medSplit.refAboveMed{i},medSplit.tarAboveMed{i},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefAboveMedVtarAboveMed = ccgR(:,1,2)/sum(ccgR(:,1,2));
    
                    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(medSplit.refAboveMed{i},medSplit.tarBelowMed{i},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefAboveMedVtarBelowMed = ccgR(:,1,2)/sum(ccgR(:,1,2));
    
    nexttile
    plot(tR,CCGrefAboveMedVtarBelowMed,'LineWidth',2)
    hold on
    plot(tR,CCGrefAboveMedVtarAboveMed,'LineWidth',2)
    hold off
    legend('CCGrefAboveMedVtarAboveMed','CCGrefAboveMedVtarBelowMed')
                    
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(medSplit.refBelowMed{i},medSplit.tarAboveMed{i},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefBelowMedVtarAboveMed = ccgR(:,1,2)/sum(ccgR(:,1,2));
    
    [GSPExc,GSPInh,pvalE,pvalI,ccgR,tR,LSPExc,LSPInh,JBSIE,JBSII,ccgjm] = ...
            CCG_jitter(medSplit.refBelowMed{i},medSplit.tarBelowMed{i},fs,binSize,duration,'jscale',jscale, ...
                        'plot_flag', false, ...
                        'plot_output', get(fig_use, 'Number'), ...
                        'njitter', njitter, 'alpha', alpha,...
                        'for_grant', for_grant);
    CCGrefBelowMedVtarBelowMed = ccgR(:,1,2)/sum(ccgR(:,1,2));
    
    nexttile
    plot(tR,CCGrefBelowMedVtarBelowMed,'LineWidth',2)
    hold on 
    plot(tR,CCGrefBelowMedVtarAboveMed,'LineWidth',2)
    hold off
    legend('CCGrefBelowMedVtarAboveMed','CCGrefBelowMedVtarBelowMed')
                    
% end