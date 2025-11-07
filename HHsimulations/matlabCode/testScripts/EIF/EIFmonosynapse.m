function EIFmonosynapse(timeMin)
    
    tic

%     close all
    
    time = timeMin*60;% in seconds
    
    noiseSc = 0.01; % normal noise scale

    Vrest = -65; % mV − change this to −65 if desired
    dt = 1/30; % ms
    fs  = (1/dt)*1000; % hz
    totalTime = time*1e3; % ms
    
    % Neuron parameters
    C = 281;                 % Membrane capacitance (pF)
    gL = 30;                 % Leak conductance (nS)
    EL = -60;              % Leak reversal potential (mV)
    DeltaT = 0.1;              % Slope factor (mV)
    VT = -50.4;              % Spike threshold (mV)
    V_reset = -65;         % Reset potential (mV)
    
    % Initialization
    V = EL;                  % Membrane potential (mV)

    % Vector of timesteps
    t = [0:dt:totalTime];
    
    % synaptic input constants
    tau = 30*2; % units in samples, example 30 units are 1ms for dt = 1/30 ms
    synapseC = 100;
    g_SynapseTar = 1*synapseC;
    E_syn = 0;
    
    delayTar = 0*fs; % delay units are in samples, coeffirefent units are is msec
    
    % synapse PSP waveform
    for i = 1:(30*30)
        sWave(i) = s_function(i,tau);
    end
    sWave = [0 sWave];
    
    % target ion channels switch
    ionSwitch = ones(1,length(t));
    ionSwitch(1:60*fs) = 0;
    
    % junction switch
    juncSwitch = ones(1,length(t));
    juncSwitch(60*fs+1:120*fs) = 0;
    
    frozenNoiseRef = noiseSc*normrnd(0,1,1,length(t));
    frozenNoiseTar = noiseSc*normrnd(0,1,1,length(t));
    
    %%
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_ref = ones(1,length(t))*0.0;
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_ref(50/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    I_current_tar(30/dt:end) = 0; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
    
    spikeTimesRef = [];
    for i = 1:length(t)

        dVdt = (-gL*(V - EL) + gL*DeltaT*exp((V - VT)/DeltaT) + I_current_ref(i))/C;
        V = V*(1+frozenNoiseRef(i)) + dVdt * dt;

        if V >= VT
            spikeTimesRef = [spikeTimesRef i];
            V = V_reset;
        end

        V_trace(i) = V;
    end
    
    %%
    
    % s function for all timesteps
    sTar = zeros(1,length(t));
    
    sTar(spikeTimesRef + delayTar) = 1;
    
    sTar = conv(sTar,sWave);
    sTar = sTar(1:length(t));
    
    % Current input −− change this to see how different inputs affect the neuron
    I_current_tar = ones(1,length(t))*0.0;
    
    I_current_tar(30/dt:end) = 0; % Input of 1 microA/cm2 beginning at 30 ms and steady until end of time period.
    
    spikeTimesTar = [];
    for i = 1:length(t)
        
        % synaptic input currents
        I_Ref_tar(i) = g_SynapseTar*sTar(i)*(V-E_syn);
        
        dVdt = (-gL*(V - EL) + gL*DeltaT*exp((V - VT)/DeltaT) - I_Ref_tar(i))/C;
        V = V*(1+frozenNoiseTar(i)) + dVdt * dt;

        if V >= VT
            spikeTimesTar = [spikeTimesTar i];
            V = V_reset;
        end

        V_trace(i) = V;
    end
    
    %%
    
   Duration = 0.100;
    
    spikeTimesRef = spikeTimesRef'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
   
    [ccg,tR] = CCG([spikeTimesRef;spikeTimesTar],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
    
     
    filterFlag   = false; 
    noDemeanFlag = true;
    fpass        = 300;
    preLength    = 5*30;
    postLength   = 5*30;
    PSP          = waveformAvg(V_trace',spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);            

%     tiledlayout(3,1)
    
    nexttile(2) 
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-5 5])
    xlabel('[ms]')
    ylabel('Prob.')
%     title('R v T')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(8)
    plot((-(preLength-1):postLength)*dt,PSP*1000,'LineWidth',1) 
%     xlim([0 2000])
    xlabel('[ms]')
    ylabel('Vm [mV]')
% 	title('pulse')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
%     set(gca, 'YDir','reverse')
    box off
    
    nexttile(14)
    plot((-(preLength-1):(postLength-1))*dt,33*diff(PSP)*1000,'-.','LineWidth',1,'color',"#0072BD")
%     xlim([0 2000])
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')
%     title('d(pulse)/dt')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
end

% synaptic gating alpha function (David Heeger: https://www.cns.nyu.edu/~david/handouts/synapse.pdf)
function s = s_function(t,tau)
    s = (t/tau)*exp(-(t/tau));
end
