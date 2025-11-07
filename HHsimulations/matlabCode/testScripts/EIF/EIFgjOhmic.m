function EIFgjOhmic(timeMin)
    
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
    
    gC = 25; %GJ ohmic coupling 
    g_GJref = 1*gC;
    g_GJtar = 1*gC;

    % Initialization
    Vref = EL;                  % Membrane potential (mV)
    Vtar = EL;
    
    % Vector of timesteps
    t = [0:dt:totalTime];
    
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
    spikeTimesTar = [];
    
    for i = 1:length(t)
        
        % ref
        if ~isempty(spikeTimesTar) && (spikeTimesTar(end) == i-1)
            dVdt = (-gL*(Vref - EL) + gL*DeltaT*exp((Vref - VT)/DeltaT) + I_current_ref(i) - g_GJtar*(Vref-10))/C;
        else
            dVdt = (-gL*(Vref - EL) + gL*DeltaT*exp((Vref - VT)/DeltaT) + I_current_ref(i) - g_GJtar*(Vref-Vtar))/C;
        end
        Vref = Vref*(1+frozenNoiseRef(i)) + dVdt * dt;

        if Vref >= VT
            spikeTimesRef = [spikeTimesRef i];
            Vref = V_reset;
        end

        Vref_trace(i) = Vref;
        
        % tar
        if ~isempty(spikeTimesRef) && spikeTimesRef(end) == i-1
            dVdt = (-gL*(Vtar - EL) + gL*DeltaT*exp((Vtar - VT)/DeltaT) + I_current_tar(i) - g_GJref*(Vtar-10))/C;
        else
            dVdt = (-gL*(Vtar - EL) + gL*DeltaT*exp((Vtar - VT)/DeltaT) + I_current_tar(i) - g_GJref*(Vtar-Vref))/C;
        end
        Vtar = Vtar*(1+frozenNoiseTar(i)) + dVdt * dt;

        if Vtar >= VT
            spikeTimesTar = [spikeTimesTar i];
            Vtar = V_reset;
        end

        Vtar_trace(i) = Vtar;
        
    end
    
    %%
    
    Duration = 0.100;
        
    spikeTimesRef = spikeTimesRef'/fs;
    spikeTimesTar = spikeTimesTar'/fs;
    
    [ccg,tR] = CCG([spikeTimesRef;spikeTimesTar],[ones(size(spikeTimesRef));2*ones(size(spikeTimesTar))], ...
                    'binSize', 1/fs, 'duration', Duration, 'Fs', 1/fs, ...
                    'norm', 'counts');
                
     
    filterFlag        = false; 
    noDemeanFlag      = true;
    fpass             = 300;
    preLength         = 5*30;
    postLength        = 5*30;
    couplingPotential =  waveformAvg(Vtar_trace',spikeTimesRef(spikeTimesRef < 60)*fs,preLength,postLength,fpass,fs,filterFlag,noDemeanFlag);        
    
    nexttile(3)
    plot(tR*1e3,ccg(:,1,2)/length(spikeTimesRef),'k','LineWidth',1)
    xlim([-5 5])
    xlabel('[ms]')
    ylabel('Prob.')
%     title('R v T')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    nexttile(9)
    plot((-(preLength-1):postLength)*dt,couplingPotential*1000,'LineWidth',1) 
    xlim([-5 5])
    xlabel('[ms]')
    ylabel('Vm [mV]')
    set(gca,'FontSize',5)
    box off
    
    nexttile(15)
    plot((-(preLength-1):(postLength-1))*dt,33*diff(couplingPotential)*1000,'-.','LineWidth',1,'color',"#0072BD")
    xlabel('[ms]')
    ylabel('dVm/dt [mV/ms]')
%     title('d(pulse)/dt')
    set(gca,'FontSize',5)
    set(gca,'FontName','Arial')
    box off
    
    
end
