function HodgkinHuxleyEpSohmic(timeMin,DC)
    
    tic

%     close all
    
    time = timeMin*60; % in seconds

    noiseSc = 3; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
%     E_Na = 115 + Vrest; % mV
%     E_K = -6 + Vrest; %mV
%     E_Leak = 10.6 + Vrest; % mV
%     
%     g_Na = 120; % mS/cm^2
%     g_K = 36; % mS/cm^2
%     g_Leak = 0.3; % mS/cm^2
    
    E_Na   = 55; % mV
    E_K    = -90; %mV
    E_Leak = -65; % mV
    
    g_Na = 35; % mS/cm^2
    g_K = 9; % mS/cm^2
    g_Leak = 0.1; % mS/cm^2

    gEpS = 0.5; %EpS coupling constant 
    
    % Vector of timesteps
    t = [0:dt:totalTime];
       
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    % sinusoidal current input
%     F = 150; % Sine wave frequency (hertz)
%     I_current = 10*sin(2*pi*F*t/1000);
    
    pseudoTrial = 30; % in ms

    %%
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_ref(50/dt:end)          = DC; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar((2*60*1000)/dt:end) = DC; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % add noise
    I_current_ref = I_current_ref + frozenNoiseRef;
    I_current_tar = I_current_tar + frozenNoiseTar;
    
    I_current_tar(1:(2*60*1000)/dt) = 0;
    
    % initializing values
    Vref(1) = Vrest; % membrane potential is starting at its resting state
    Vtar(1) = Vrest;
    
    % separate functions to get the alpha and beta values
    [alphaMref, betaMref] = m_equations(Vref(1), Vrest);
    [alphaNref, betaNref] = n_equations(Vref(1), Vrest);
    [alphaHref, betaHref] = h_equations(Vref(1), Vrest);
    
    [alphaMtar, betaMtar] = m_equations(Vtar(1), Vrest);
    [alphaNtar, betaNtar] = n_equations(Vtar(1), Vrest);
    [alphaHtar, betaHtar] = h_equations(Vtar(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    mRef(1) = (alphaMref / (alphaMref + betaMref));
    nRef(1) = (alphaNref / (alphaNref + betaNref));
    hRef(1) = (alphaHref / (alphaHref + betaHref));
    
    mTar(1) = (alphaMtar / (alphaMtar + betaMtar));
    nTar(1) = (alphaNtar / (alphaNtar + betaNtar));
    hTar(1) = (alphaHtar / (alphaHtar + betaHtar));
    
    countSample = 0;
    juncOff     = false;
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % reset like a trial
        if mod(i,30*pseudoTrial) == 0 % 30 samples x # ms
            % initializing values
            Vref(i) = Vrest; % membrane potential is starting at its resting state
            Vtar(i) = Vrest;

            % separate functions to get the alpha and beta values
            [alphaMref, betaMref] = m_equations(Vref(i), Vrest);
            [alphaNref, betaNref] = n_equations(Vref(i), Vrest);
            [alphaHref, betaHref] = h_equations(Vref(i), Vrest);

            [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
            [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
            [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);

            % initializing gating variables to the asymptotic values when membrane potential
            % is set to the membrane resting value based on equation 13
            mRef(i) = (alphaMref / (alphaMref + betaMref));
            nRef(i) = (alphaNref / (alphaNref + betaNref));
            hRef(i) = (alphaHref / (alphaHref + betaHref));

            mTar(i) = (alphaMtar / (alphaMtar + betaMtar));
            nTar(i) = (alphaNtar / (alphaNtar + betaNtar));
            hTar(i) = (alphaHtar / (alphaHtar + betaHtar));
        else
            % calculate new alpha and beta based on last known membrane potenatial
            [alphaMref, betaMref] = m_equations(Vref(i), Vrest);
            [alphaNref, betaNref] = n_equations(Vref(i), Vrest);
            [alphaHref, betaHref] = h_equations(Vref(i), Vrest);

            [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
            [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
            [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);
        end
        
        % for pseudo-trials reset turn junctions off for half of trial
        % length
        if sum(mod(i,30*pseudoTrial) == 0:((30*(pseudoTrial/2))-1)) == 1 
            juncSwitch(i) = 0;
        end 
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_ref(i)  = g_K*(nRef(i)^4);
        conductance_Na_ref(i) = g_Na*(mRef(i)^3)*hRef(i);
        conductance_K_tar(i)  = g_K*(nTar(i)^4);
        conductance_Na_tar(i) = g_Na*(mTar(i)^3)*hTar(i);
        
        % retrieving ionic currents
        I_Na_ref(i)   = conductance_Na_ref(i)*(Vref(i)-E_Na);
        I_K_ref(i)    = conductance_K_ref(i)*(Vref(i)-E_K);
        I_Leak_ref(i) = g_Leak*(Vref(i)-E_Leak);
        I_Na_tar(i)   = conductance_Na_tar(i)*(Vtar(i)-E_Na);
        I_K_tar(i)    = conductance_K_tar(i)*(Vtar(i)-E_K);
        I_Leak_tar(i) = g_Leak*(Vtar(i)-E_Leak);
        
        % Calculating the input
        Input_ref = I_current_ref(i) - (I_Na_ref(i) + I_K_ref(i) + I_Leak_ref(i));
        if i == 1
            Input_tar = I_current_tar(i) - (I_Na_tar(i) + I_K_tar(i) + I_Leak_tar(i));
        else
            Input_tar = I_current_tar(i) - (I_Na_tar(i) + I_K_tar(i) + I_Leak_tar(i)) + juncSwitch(i)*gEpS*Input_ref; % gEpS*(Vref(i)-Vref(i-1)); 
            % Input_tar = I_current_tar(i) - (I_Na_tar(i) + I_K_tar(i) + I_Leak_tar(i)); + ~juncOff*juncSwitch(i)*gEpS*Input_ref; % gEpS*(Vref(i)-Vref(i-1)); 
        end
        
        InputTarCurrent(i) = Input_tar;
        
        % Calculating the new membrane potential
%         Vref(i+1) = Vref(i)*(1 + frozenNoiseRef(i)) + Input_ref* dt*(1/C);
%         Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + Input_tar* dt*(1/C);

        Vref(i+1) = Vref(i) + Input_ref* dt*(1/C);
        Vtar(i+1) = Vtar(i) + Input_tar* dt*(1/C);
        
        % getting new values for the gating variables
        mRef(i+1) = mRef(i) + (alphaMref *(1-mRef(i)) - betaMref * mRef(i))*dt;
        nRef(i+1) = nRef(i) + (alphaNref *(1-nRef(i)) - betaNref * nRef(i))*dt;
        hRef(i+1) = hRef(i) + (alphaHref *(1-hRef(i)) - betaHref * hRef(i))*dt;
        mTar(i+1) = mTar(i) + (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i))*dt;
        nTar(i+1) = nTar(i) + (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i))*dt;
        hTar(i+1) = hTar(i) + (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i))*dt;
        
        % turn off junction after ref spike
%         if (Vref(i+1) > 25) || ~juncOff
%             juncOff = true;
%         end
%         if juncOff
%             countSample = countSample + 1;
%             if countSample == 90
%                 juncOff = false;
%                 countSample = 0;
%             end
%         end
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimesRef] = findpeaks(Vref,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    [~,spikeTimesTar] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
        
    %%
    
    Duration = 0.100;
        
    spikeTimesRef = spikeTimesRef'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
    
    refFR = length(spikeTimesRef)/(spikeTimesRef(end)-spikeTimesRef(1));
    tarFR = length(spikeTimesTar)/(spikeTimesTar(end)-spikeTimesTar(1));
    
    [ccg,tR] = CCG([spikeTimesRef;spikeTimesTar],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
    filterFlag     = false; 
    noDemeanFlag   = true;
    fpass          = 300;
    preLength      = 5*30;
    postLength     = 5*30;
    refAP          = waveformAvg(Vref',            spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);
    refDiffAP      = waveformAvg(diff(Vref'),      spikeTimesRef(1)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);        
    evokedSpikelet = waveformAvg(Vtar',            spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);        
    stimWaveInput  = waveformAvg(InputTarCurrent', spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag); 
    
%     figure
%     tiledlayout(3,1)
    
    nexttile(3)
    plot((-(preLength-1):postLength)*dt,refDiffAP*1000,'k','LineWidth',1)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')
    title('Ref. extrawave')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off

%     nexttile(4)
    nexttile(4)
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    ylims = get(gca,'ylim');
    text(-1.85,0.65*ylims(2),{['ref ' num2str(refFR,'%.1f') 'hz']; ...
                              ['tar ' num2str(tarFR,'%.1f') 'hz']; ...
                              ['dc: ' num2str(DC)];['gain ' num2str(gEpS)]; ...
                              ['trial: ' num2str(pseudoTrial) 'ms']},'FontSize',5)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
    title('CCG of extrawave as input')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
%     plot((-(preLength-1):postLength)*dt,33*refDiffAP*1000,'LineWidth',1) 
%     xlim([-5 5])
%     xlabel('[ms]')
%     title('refDiffAP')
%     ylabel('dVm/dt [mV/ms]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
    
%     nexttile
%     plot(((-(preLength-1):(postLength))*dt)-1,stimWaveInput,'-.','LineWidth',1,'color',"#0072BD")
%     xlim([-1 1])
%     xlabel('[ms]')
%     title('stimWaveInput')
%     ylabel('[mA/ch^2]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off

    
% %     nexttile(9)
%     nexttile
%     plot((-(preLength-1):postLength)*dt,evokedSpikelet*1000,'LineWidth',1) 
%     xlim([-1 1])
%     xlabel('[ms]')
%     ylabel('Vm [mV]')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
%     
% %     nexttile(14)
%     nexttile
%     plot((-(preLength-1):(postLength-1))*dt,33*diff(evokedSpikelet)*1000,'-.','LineWidth',1,'color',"#0072BD")
%     xlim([-1 1])
%     xlabel('[ms]')
%     ylabel('dVm/dt [mV/ms]')  
%     title('dPSP')
% %     title('d(pulse)/dt')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
    
%     
%     figure
%     yyaxis left
%     plot(Vref(1:1e4))
%     hold on
%     plot(Vtar(1:1e4))
%     hold off
% 
%     yyaxis right
%     plot(diff(Vref(1:10001)))
    
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccg(:,1,1))
%     nexttile; plot(tR*1000,ccg(:,1,2))
%     nexttile; plot(tR*1000,ccg(:,2,1))
%     nexttile; plot(tR*1000,ccg(:,2,2))
    
%     plot(tR*1000,ccg(:,1,2),'LineWidth',2)
%     xlim([-5 5])
%     xlabel('[ms]')
%     ylabel('Vm [mV]')
%     title('T v R')
%     set(gca,'FontSize',6)

end

% calculate alpha m and beta m based on Table 2
function [alpha_m, beta_m] = m_equations(V, Vrest)
    alpha_m = (2.5-0.1*(V-Vrest))/(exp(2.5-0.1*(V-Vrest))-1);
    beta_m = 4*exp((Vrest-V)/18);
end

% calculate alpha n and beta n based on Table 2
function [alpha_n, beta_n] = n_equations(V, Vrest)
    alpha_n = (0.1-0.01*(V-Vrest))/(exp(1-0.1*(V-Vrest))-1);
    beta_n = 0.125*exp((Vrest-V)/80);
end

% calculate alp3ha h and beta h based on Table 2
function [alpha_h, beta_h] = h_equations(V, Vrest)
    alpha_h = 0.07*exp((Vrest-V)/20);
    beta_h = 1/(1+exp(3-0.1*(V-Vrest)));
end