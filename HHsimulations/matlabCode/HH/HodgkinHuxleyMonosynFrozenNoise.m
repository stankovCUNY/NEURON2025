function HodgkinHuxleyMonosynFrozenNoise(timeMin)
    
    tic

%     close all
    
    time = timeMin*60; % in seconds

    noiseSc = 0.02; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    C = 1; % uF/cm^2
    
    % constants; values based on Table 1
    E_Na = 115 + Vrest; % mV
    E_K = -6 + Vrest; %mV
    E_Leak = 10.6 + Vrest; % mV
    
    g_Na = 120; % mS/cm^2
    g_K = 36; % mS/cm^2
    g_Leak = 0.3; % mS/cm^2
        
    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % synaptic input constants
    tau = 30*2; % units in samples, example 30 units are 1ms for dt = 1/30 ms
    synapseC = 0.10;
    g_SynapseTar = 1*synapseC;
    E_syn = 0;
    
    delayTar = 0*fs; % delay units are in samples, coeffirefent units are is msec
    
    % synapse PSP waveform
    for i = 1:(30*30)
        sWave(i) = s_function(i,tau);
    end
    sWave = [0 sWave];
    
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    % target ion channels switch
    ionSwitch = ones(1,length(t));
    ionSwitch(1:60*fs) = 0;
    
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    %%
    
    % initializing values
    Vref(1) = Vrest; % membrane potential is starting at its resting state
    Vtar(1) = Vrest;
    s(1)    = 0; % initial synaptic gain
    
    % separate functions to get the alpha and beta values
    [alphaMref, betaMref] = m_equations(Vref(1), Vref);
    [alphaNref, betaNref] = n_equations(Vref(1), Vref);
    [alphaHref, betaHref] = h_equations(Vref(1), Vref);
    
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
    
    % solve for common input spike times
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMref, betaMref]   = m_equations(Vref(i), Vrest);
        [alphaNref, betaNref]   = n_equations(Vref(i), Vrest);
        [alphaHref, betaHref]   = h_equations(Vref(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_Ref(i)   = g_K*(nRef(i)^4);
        conductance_Na_Ref(i)  = g_Na*(mRef(i)^3)*hRef(i);
        
        % retrieving ionic currents
        I_Na_Ref(i)   = conductance_Na_Ref(i)*(Vref(i)-E_Na);
        I_K_Ref(i)    = conductance_K_Ref(i)*(Vref(i)-E_K);
        I_Leak_Ref(i) = g_Leak*(Vref(i)-E_Leak);
        
        % Calculating the input
        Input_Ref  = - (I_Na_Ref(i) + I_K_Ref(i) + I_Leak_Ref(i));
        
        % Calculating the new membrane potential
        Vref(i+1)  = Vref(i)*(1 + frozenNoiseRef(i))  + Input_Ref*  dt*(1/C); 
        
        % getting new values for the gating variables
        mRef(i+1)  = mRef(i) + (alphaMref *(1-mRef(i)) - betaMref * mRef(i))*dt;
        nRef(i+1)  = nRef(i) + (alphaNref *(1-nRef(i)) - betaNref * nRef(i))*dt;
        hRef(i+1)  = hRef(i) + (alphaHref *(1-hRef(i)) - betaHref * hRef(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 10;   % samples
    [~,spikeTimesRef] = findpeaks(Vref,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
    % s function for all timesteps
    sTar = zeros(1,length(t));
    
    sTar(spikeTimesRef + delayTar) = 1;
    
    sTar = conv(sTar,sWave);
    sTar = sTar(1:length(t));
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_tar(30/dt:end) = 0; % Input of 1 microA/cm2 beginning at 30 ms and steady until end of time period.
    
    %% synapse off
    
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
        [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
        [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_tar(i)  = g_K*(nTar(i)^4);
        conductance_Na_tar(i) = g_Na*(mTar(i)^3)*hTar(i);
        
        % retrieving ionic currents
        I_Na_tar(i)   = conductance_Na_tar(i)*(Vtar(i)-E_Na);
        I_K_tar(i)    = conductance_K_tar(i)*(Vtar(i)-E_K);
        I_Leak_tar(i) = g_Leak*(Vtar(i)-E_Leak);
       
        % Calculating the input
        Input_tar =  - (ionSwitch(i)*(I_Na_tar(i) + I_K_tar(i)) + I_Leak_tar(i));
        
        % Calculating the new membrane potential
        Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + Input_tar* dt*(1/C);
        
        % getting new values for the gating variables
        mTar(i+1) = mTar(i) + (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i))*dt;
        nTar(i+1) = nTar(i) + (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i))*dt;
        hTar(i+1) = hTar(i) + (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 10;   % samples
    [~,spikeTimesTarOff] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
        
    %%
    
    for i = 1:length(t)
        
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaMtar, betaMtar] = m_equations(Vtar(i), Vrest);
        [alphaNtar, betaNtar] = n_equations(Vtar(i), Vrest);
        [alphaHtar, betaHtar] = h_equations(Vtar(i), Vrest);
        
        % conductance variables − computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K_tar(i)  = g_K*(nTar(i)^4);
        conductance_Na_tar(i) = g_Na*(mTar(i)^3)*hTar(i);
        
        % retrieving ionic currents
        I_Na_tar(i)   = conductance_Na_tar(i)*(Vtar(i)-E_Na);
        I_K_tar(i)    = conductance_K_tar(i)*(Vtar(i)-E_K);
        I_Leak_tar(i) = g_Leak*(Vtar(i)-E_Leak);
        
        % common synaptic input currents
        I_Ref_tar(i) = g_SynapseTar*sTar(i)*(Vtar(i)-E_syn);
        
        % Calculating the input
        Input_tar = I_current_tar(i) - juncSwitch(i)*I_Ref_tar(i) - (ionSwitch(i)*(I_Na_tar(i) + I_K_tar(i)) + I_Leak_tar(i));
        
        % Calculating the new membrane potential
        Vtar(i+1) = Vtar(i)*(1 + frozenNoiseTar(i)) + Input_tar* dt*(1/C);
        
        % getting new values for the gating variables
        mTar(i+1) = mTar(i) + (alphaMtar *(1-mTar(i)) - betaMtar * mTar(i))*dt;
        nTar(i+1) = nTar(i) + (alphaNtar *(1-nTar(i)) - betaNtar * nTar(i))*dt;
        hTar(i+1) = hTar(i) + (alphaHtar *(1-hTar(i)) - betaHtar * hTar(i))*dt;
        
    end
    
    % spike times
    thresh   = 0;    % volts
    distance = 10;   % samples
    [~,spikeTimesTarOn] = findpeaks(Vtar,'MinPeakHeight',thresh,'MinPeakDistance',distance);
    
    %%
    
    Duration = 0.100;
    
    spikeTimesRef    = spikeTimesRef'/fs;
    spikeTimesTarOff = spikeTimesTarOff'/fs;
    spikeTimesTarOn  = spikeTimesTarOn'/fs;
   
    [ccgOff,tR] = CCG([spikeTimesRef;spikeTimesTarOff],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTarOff))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
    [ccgOn,tR] = CCG([spikeTimesRef;spikeTimesTarOn],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTarOn))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
     
%     filterFlag   = false; 
%     noDemeanFlag = true;
%     fpass        = 300;
%     preLength    = 5*30;
%     postLength   = 5*30;
%     PSP          = waveformAvg(Vtar',spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

    tiledlayout(1,3)
    
    nexttile 
    plot(tR*1e3,ccgOff(:,1,2),'k','LineWidth',2)
    hold on
    plot(tR*1e3,ccgOn(:,1,2),'r','LineWidth',2)
    hold off
    xlim([-20 20])
    xlabel('lag times [sec]')
    ylabel('Counts')
    legend('off','on')
    title('R v T')
%     set(gca,'FontSize',6)
    
    spikeTimesRef = spikeTimesRef(spikeTimesRef > 120);
    spikeTimesTarOn = spikeTimesTarOn(spikeTimesTarOn > 120);
    spikeTimesTarOff = spikeTimesTarOff(spikeTimesTarOff > 120);
    
    nexttile    
    raster(spikeTimesRef,spikeTimesTarOn,spikeTimesTarOff,20)
    
%     nexttile(7)
%     plot((-(preLength-1):postLength)*dt*1000,PSP*1000,'LineWidth',2) 
%     xlim([0 2000])
%     xlabel('[\musec]')
%     ylabel('Vm [mV]')
% 	title('pulse')
%     set(gca,'FontSize',6)
%     set(gca, 'YDir','reverse')
%     
%     nexttile(12)
%     plot((-(preLength-1):(postLength-1))*dt*1000,33*diff(PSP)*1000,'-.','LineWidth',2,'color',"#0072BD")
%     xlim([0 2000])
%     xlabel('[\musec]')
%     ylabel('dVm/dt [mV/\musec]')
%     title('d(pulse)/dt')
%     set(gca,'FontSize',6)
    
    nexttile    
    k = 0.050;
    
    % 
    D = [];
    for j = 1:length(spikeTimesRef)

        Xtemp = spikeTimesTarOn - spikeTimesRef(j);
        Xtemp = Xtemp((Xtemp > 0) & (Xtemp < k));

        Ytemp = spikeTimesTarOff - spikeTimesRef(j) ;
        Ytemp = Ytemp((Ytemp > 0) & (Ytemp < k));

        if isempty(Ytemp) || isempty(Xtemp)
            continue
        end

        D(j) = Ytemp(1) - Xtemp(1);
    end

    D(D == 0) = [];
    histogram(D*1000,'binWidth',1/3e4)
    xlabel('Time Shift [ms]')
    ylabel('Counts')
    title('D = X-Y (OFF-ON)')
    
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

function raster(ref,tarOn,tarOff,lgmx)
    
    % convert to ms
    ref    = ref*1000;
    tarOn  = tarOn*1000;
    tarOff = tarOff*1000;
    
    xR = []; yR = []; xTon = []; yTon = [];  xToff = []; yToff = [];  
    for k = 1:numel(ref)
        rtem    = ref(ref>=ref(k)-lgmx & ref<=ref(k)+lgmx) - ref(k);
        ttemOn  = tarOn(tarOn >=ref(k)-lgmx & tarOn <=ref(k)+lgmx) - ref(k);
        ttemOff = tarOff(tarOff>=ref(k)-lgmx & tarOff<=ref(k)+lgmx) - ref(k);
        xR      = [xR;rtem];
        yR      = [yR;k.*ones(size(rtem))];
        xTon    = [xTon;ttemOn];
        yTon    = [yTon;k.*ones(size(ttemOn))];
        xToff   = [xToff;ttemOff];
        yToff   = [yToff;k.*ones(size(ttemOff))];
    end
    
    scatter(xToff,yToff,1,'filled','markerfacecolor','k'); hold on;
    hold on
    scatter(xTon,yTon,1,'filled','markerfacecolor','r'); hold on;
    hold off
    xlabel('Lag times [sec]')
    ylabel('Trial')
    xlim([-lgmx,lgmx])
    title('RASTER: Target spikes aligned with reference spikes')

end