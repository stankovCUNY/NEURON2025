% CCGwidth = 0.030; % CCG CCGwidth

% filtered all states
runSpWaveAnaMedianSplit('troughAmplitude','allStates')
runSpWaveAnaMedianSplit('peakToTrough',   'allStates')
runSpWaveAnaMedianSplit('halfWidth',      'allStates')
runSpWaveAnaMedianSplit('spikeWidth',     'allStates')

% filtered with no ripple 
runSpWaveAnaMedianSplit('troughAmplitude','noRipple')
runSpWaveAnaMedianSplit('peakToTrough',   'noRipple')
runSpWaveAnaMedianSplit('halfWidth',      'noRipple')
runSpWaveAnaMedianSplit('spikeWidth',     'noRipple')

% ripple only
runSpWaveAnaMedianSplit('troughAmplitude','ripple')
runSpWaveAnaMedianSplit('peakToTrough',   'ripple')
runSpWaveAnaMedianSplit('halfWidth',      'ripple')
runSpWaveAnaMedianSplit('spikeWidth',     'ripple')

% filtered with no theta and gamma 
runSpWaveAnaMedianSplit('troughAmplitude','noThetaAndRipple')
runSpWaveAnaMedianSplit('peakToTrough',   'noThetaAndRipple')
runSpWaveAnaMedianSplit('halfWidth',      'noThetaAndRipple')
runSpWaveAnaMedianSplit('spikeWidth',     'noThetaAndRipple')

% filtered with no theta and gamma 
% runSpWaveAnaMedianSplit('troughAmplitude','noThetaAndGamma')
% runSpWaveAnaMedianSplit('peakToTrough',   'noThetaAndGamma')
% runSpWaveAnaMedianSplit('halfWidth',      'noThetaAndGamma')
% runSpWaveAnaMedianSplit('spikeWidth',     'noThetaAndGamma')

% filtered with theta 
% runSpWaveAnaMedianSplit('troughAmplitude','theta')
% runSpWaveAnaMedianSplit('peakToTrough',   'theta')
% runSpWaveAnaMedianSplit('halfWidth',      'theta')
% runSpWaveAnaMedianSplit('spikeWidth',     'theta')

% filtered with gamma 
runSpWaveAnaMedianSplit('troughAmplitude','gamma')
runSpWaveAnaMedianSplit('peakToTrough',   'gamma')
runSpWaveAnaMedianSplit('halfWidth',      'gamma')
runSpWaveAnaMedianSplit('spikeWidth',     'gamma')


% filtered with ripple 
% runSpWaveAnaMedianSplit('troughAmplitude','ripple')
% runSpWaveAnaMedianSplit('peakToTrough',   'ripple')
% runSpWaveAnaMedianSplit('halfWidth',      'ripple')
% runSpWaveAnaMedianSplit('spikeWidth',     'ripple')

