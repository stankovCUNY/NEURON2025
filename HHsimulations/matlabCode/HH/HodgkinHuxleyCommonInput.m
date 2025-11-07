function HodgkinHuxleyCommonInput(timeMin,refPSPtype,tarPSPtype,delayRef,delayTar,DC)
    
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

    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % synaptic input constants
    tau = 30*2; % units in samples, example 30 units are 1ms for dt = 1/30 ms
    synapseC = 0.5;
    g_SynapseRef = synapseC;
    g_SynapseTar = synapseC;
%     E_syn = 0;
    
    % define E syn for ref and tar
    if refPSPtype == 1
        E_synRef = 0;
    elseif refPSPtype == -1
        E_synRef = -70;
    end
    if tarPSPtype == 1
        E_synTar = 0;
    elseif tarPSPtype == -1
        E_synTar = -70;
    end
    
    delayRef = delayRef*fs; % delay units are in samples, coefficient units are is msec
    delayTar = delayTar*fs; % delay units are in samples, coefficient units are is msec
    
    % synapse PSP waveform
    for i = 1:(30*30)
        sWave(i) = s_function(i,tau);
    end
    sWave = [0 sWave];
    
    frozenNoiseCI  = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    pseudoTrial = 100; % in ms
    
    %%
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_CI  = ones(1,length(t))*0.0; 
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_CI(50/dt:end)           = DC;
    I_current_ref(50/dt:end)          = DC; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar((2*60*1000)/dt:end) = DC; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    % add noise
    I_current_CI  = I_current_CI  + frozenNoiseCI;
    I_current_ref = I_current_ref + frozenNoiseRef;
    I_current_tar = I_current_tar + frozenNoiseTar;
    
    % initializing values
    Vci(1)  = Vrest; % membrane potential is starting at its resting state
    Vref(1) = Vrest; 
    Vtar(1) = Vrest;
    s(1)    = 0; % initial synaptic gain
    
    % separate functions to get the alpha and beta values
    [alphaMci, betaMci]   = m_equations(Vci(1), Vci);
    [alphaNci, betaNci]   = n_equations(Vci(1), Vci);
    [alphaHci, betaHci]   = h_equations(Vci(1), Vci);
    
    [alphaMref, betaMref] = m_equations(Vref(1), Vrest);
    [alphaNref, betaNref] = n_equations(Vref(1), Vrest);
    [alphaHref, betaHref] = h_equations(Vref(1), Vrest);
    
    [alphaMtar, betaMtar] = m_equations(Vtar(1), Vrest);
    [alphaNtar, betaNtar] = n_equations(Vtar(1), Vrest);
    [alphaHtar, betaHtar] = h_equations(Vtar(1), Vrest);
    
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    mCI(1)  = (alphaMci / (alphaMci + betaMci));
    nCI(1)  = (alphaNci / (alphaNci + betaNci));
    hCI(1)  = (alphaHci / (alphaHci + betaHci));
    
    mRef(1) = (alphaMref / (alphaMref + betaMref));
    nRef(1) = (alphaNref / (alphaNref + betaNref));
    hRef(1) = (alphaHref / (alphaHref + betaHref));
    
    mTar(1) = (alphaMtar / (alphaMtar + betaMtar));
    nTar(1) = (alphaNtar / (alphaNtar + betaNtar));
    hTar(1) = (alphaHtar / (alphaHtar + betaHtar));
    
    % solve for common input spike times
    for i = 1:length(t)
        
        % reset like a trial
        if mod(i,30*pseudoTrial) == 0 % 30 samples x # ms
            % initializing values
            Vci(i)  = Vrest; % membrane potential is starting at its resting state

            % separate functions to get the alpha and beta values
            [alphaMci, betaMci]   = m_equations(Vci(i), Vrest);
            [alphaNci, betaNci]   = n_equations(Vci(i), Vrest);
            [alphaHci, betaHci]   = h_equations(Vci(i), Vrest);

            % initializing gating variables to the asymptotic values when membrane potential
            % is set to the membrane resting value based on equation 13
            mCI(i)  = (alphaMci / (alphaMci + betaMci));
            nCI(i)  = (alphaNci / (alphaNci + betaNci));
            hCI(i)  = (alphaHci / (alphaHci + betaHci));
        else
            % calculate new alpha and beta based on last known membrane potenatial
            [alphaMci, betaMci]   = m_equations(Vci(i), Vrest);
            [alphaNci, betaNci]   = n_equations(Vci(i), Vrest);
            [alphaHci, betaHci]   = h_equations(Vci(i), Vrest);
        end
        
%         for pseudo-trials reset turn junctions off for half of trial
%         length
        if sum(mod(i,30*pseudoTrial) == 0:((30*(pseudoTrial/2))-1)) == 1 
            juncSwitch(i) = 0;
        end 
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMci, betaMci]   = m_equations(Vci(i), Vrest);
        [alphaNci, betaNci]   = n_equations(Vci(i), Vrest);
        [alphaHci, betaHci]   = h_equations(Vci(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_CI(i)   = g_K*(nCI(i)^4);
        conductance_Na_CI(i)  = g_Na*(mCI(i)^3)*hCI(i);
        
        % retrieving ionic currents
        I_Na_CI(i)   = conductance_Na_CI(i)*(Vci(i)-E_Na);
        I_K_CI(i)    = conductance_K_CI(i)*(Vci(i)-E_K);
        I_Leak_CI(i) = g_Leak*(Vci(i)-E_Leak);
        
        % Calculating the input
%         Input_CI  = - (I_Na_CI(i) + I_K_CI(i) + I_Leak_CI(i));
        Input_CI = I_current_CI(i) - (I_Na_CI(i) + I_K_CI(i) + I_Leak_CI(i));
        
        % Calculating the new membrane potential
%         Vci(i+1)  = Vci(i)*(1 + frozenNoiseCI(i))  + Input_CI*  dt*(1/C); 
        Vci(i+1)  = Vci(i) + Input_CI*dt*(1/C); 
        
        % getting new values for the gating variables
        mCI(i+1)  = mCI(i) + (alphaMci *(1-mCI(i)) - betaMci * mCI(i))*dt;
        nCI(i+1)  = nCI(i) + (alphaNci *(1-nCI(i)) - betaNci * nCI(i))*dt;
        hCI(i+1)  = hCI(i) + (alphaHci *(1-hCI(i)) - betaHci * hCI(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimesCI] = findpeaks(Vci,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
    % s function for all timesteps
    sRef = zeros(1,length(t));
    sTar = zeros(1,length(t));
    
    sRef(spikeTimesCI + delayRef) = 1;
    sTar(spikeTimesCI + delayTar) = 1;
    
    sRef = conv(sRef,sWave);
    sRef = sRef(1:length(t));
    
    sTar = conv(sTar,sWave);
    sTar = sTar(1:length(t));
    
    % target ion channels switch
    ionSwitch = ones(1,length(t));
    ionSwitch(1:60*fs) = 0;
    
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % reset like a trial
        if mod(i,30*pseudoTrial) == 0 % 30 samples x # ms
            % initializing values
            % membrane potential is starting at its resting state
            Vref(i) = Vrest; 
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
        
        % common synaptic input currents
        I_CI_ref(i) = g_SynapseRef*sRef(i)*(Vref(i)-E_synRef);
        I_CI_tar(i) = g_SynapseTar*sTar(i)*(Vtar(i)-E_synTar);
        
        % Calculating the input
        Input_ref = I_current_ref(i) - (ionSwitch(i)*(I_Na_ref(i) + I_K_ref(i)) + I_Leak_ref(i) + juncSwitch(i)*I_CI_ref(i));
        Input_tar = I_current_tar(i) - (ionSwitch(i)*(I_Na_tar(i) + I_K_tar(i)) + I_Leak_tar(i) + juncSwitch(i)*I_CI_tar(i));
        
        % Calculating the new membrane potential
%         Vref(i+1) = Vref(i)*(1 + frozenNoiseRef(i)) + Input_ref* dt*(1/C);
%         Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + Input_tar* dt*(1/C);
        Vref(i+1) = Vref(i) + Input_ref*dt*(1/C);
        Vtar(i+1) = Vtar(i) + Input_tar*dt*(1/C);
        
        % getting new values for the gating variables
        mRef(i+1) = mRef(i) + (alphaMref *(1-mRef(i)) - betaMref * mRef(i))*dt;
        nRef(i+1) = nRef(i) + (alphaNref *(1-nRef(i)) - betaNref * nRef(i))*dt;
        hRef(i+1) = hRef(i) + (alphaHref *(1-hRef(i)) - betaHref * hRef(i))*dt;
        
        mTar(i+1) = mTar(i) + (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i))*dt;
        nTar(i+1) = nTar(i) + (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i))*dt;
        hTar(i+1) = hTar(i) + (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 50;   % samples
    [~,spikeTimesRef] = findpeaks(Vref,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    [~,spikeTimesTar] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
        
    %%
    
    Duration = 0.100;
    
    spikeTimesCI  = spikeTimesCI'/fs;
    spikeTimesRef = spikeTimesRef'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
    
    [ccgRT,tR] = CCG([spikeTimesRef;spikeTimesTar],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
    CiFR  = length(spikeTimesCI)/ (spikeTimesCI(end) -spikeTimesCI(1));
    refFR = length(spikeTimesRef)/(spikeTimesRef(end)-spikeTimesRef(1));
    tarFR = length(spikeTimesTar)/(spikeTimesTar(end)-spikeTimesTar(1));
    
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgRT(:,1,1))
%     title("R v R")
%     nexttile; plot(tR*1000,ccgRT(:,1,2))
%     title("R v T")
%     nexttile; plot(tR*1000,ccgRT(:,2,1))
%     title("T v R")
%     nexttile; plot(tR*1000,ccgRT(:,2,2))
%     title("T v T")
    
    [ccgCiR,tR] = CCG([spikeTimesCI;spikeTimesRef],[ones(size(spikeTimesCI));2*ones(size(spikeTimesRef))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgCiR(:,1,1))
%     title("CI v CI")
%     nexttile; plot(tR*1000,ccgCiR(:,1,2))
%     title("CI v R")
%     nexttile; plot(tR*1000,ccgCiR(:,2,1))
%     title("R v CI")
%     nexttile; plot(tR*1000,ccgCiR(:,2,2))
%     title("R v R")
    
    [ccgCiT,tR] = CCG([spikeTimesCI;spikeTimesTar],[ones(size(spikeTimesCI));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
%     figure
%     tiledlayout(2,2)            
%     nexttile; plot(tR*1000,ccgCiT(:,1,1))
%     title('CI v CI')
%     nexttile; plot(tR*1000,ccgCiT(:,1,2))
%     title("CI v T")
%     nexttile; plot(tR*1000,ccgCiT(:,2,1))
%     title("T v CI")
%     nexttile; plot(tR*1000,ccgCiT(:,2,2))
%     title("T v T")
%     
%     figure
%     yyaxis left
%     plot(Vref(1:1e4))
%     hold on
%     plot(Vtar(1:1e4))
%     hold off
% 
%     yyaxis right
%     plot(I_CI_ref(1:1e4))
%     hold on
%     plot(I_CI_tar(1:1e4))
%     hold off
%     
%     title("left ref/tar Vm, right ref/tar EPSP")
    
%     tiledlayout(1,2,'Padding','none','TileSpacing','compact')
    
%     nexttile 
%     plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesCI),'LineWidth',1)
%     hold on
%     plot(tR*1000,ccgCiR(:,1,2)/length(spikeTimesCI),'LineWidth',1)
%     hold off
%     xlim([-10 10])
%     xlabel('[ms]')
%     ylabel('Spike Probability')
% %     title('CI v T and CI v R')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
    
%     nexttile(3)
%     plot(tR*1000,ccgCiT(:,1,2)/length(spikeTimesCI),'k','LineWidth',1)
%     ylims = get(gca,'ylim');
%     text(-2.85,0.65*ylims(2),{['Ci '  num2str(CiFR,'%.1f') 'hz'];         ...
%                               ['ref ' num2str(refFR,'%.1f') 'hz'];        ...
%                               ['tar ' num2str(tarFR,'%.1f') 'hz'];        ...
%                               ['dc: ' num2str(DC)];['gain ' num2str(synapseC)]; ... 
%                               ['trial: ' num2str(pseudoTrial) 'ms']},'FontSize',5)
%     xlim([-3 3])
%     xlabel('[ms]')
%     ylabel('Spike Probability')
% %     title('R v T')
%     title('CCG')
%     set(gca,'FontSize',5)
%     set(gca,'FontName','Arial')
%     box off
    
    nexttile(2)
    plot(tR*1000,ccgRT(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    ylims = get(gca,'ylim');
    text(-1.85,0.65*ylims(2),{['Ci '  num2str(CiFR,'%.1f') 'hz'];         ...
                              ['ref ' num2str(refFR,'%.1f') 'hz'];        ...
                              ['tar ' num2str(tarFR,'%.1f') 'hz'];        ...
                              ['dc: ' num2str(DC)];['gain ' num2str(synapseC)]; ... 
                              ['trial: ' num2str(pseudoTrial) 'ms']},'FontSize',5)
    xlim([-2 2])
    xlabel('[ms]')
    ylabel('Spike Probability')
%     title('R v T')
    title('CCG')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    toc
    
    save(['/media/nasko/WD_BLACK31/SimTemp/HHciMSdc' num2str(DC) '.mat'])
    
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

% calculate alpha h and beta h based on Table 2
function [alpha_h, beta_h] = h_equations(V, Vrest)
    alpha_h = 0.07*exp((Vrest-V)/20);
    beta_h = 1/(1+exp(3-0.1*(V-Vrest)));
end

% synaptic gating alpha function (David Heeger: https://www.cns.nyu.edu/~david/handouts/synapse.pdf)
function s = s_function(t,tau)
    s = (t/tau)*exp(-(t/tau));
end
