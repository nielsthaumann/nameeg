function [onsets, ons] = name_ons(audio, varargin)

% Automatic sound onset detection with noise suppression
% 
% Use as: 
% 
% onsets = name_ons(audio)  
% 
% to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'
% 
% The onsets are given as an array of time points in seconds and saved to an output file 
% (by default with the same path and name as the input file added the ending ..._onsets.mat)
% 
% 
% The algorithm includes suppression of noise in the audio medium and vibrato/tremolo 
% and the ability to detect slow attacks
% 
% 
% Optional arguments for time frequency representation (TFR):
% 
% name_ons(..., 'framedur', framedur)   specifies the time frame duration in milliseconds 
%                                       for the TFR analysis (default = 100)
% 
% name_ons(..., 'srtfr', srtfr)         specifies the sampling rate for the TFR in Hz (default = 100)
%
% 
% Optional arguments for audio medium noise suppression:
% 
% name_ons(..., 'ns', ns)         defines whether audio medium noise suppression should be applied, 
%                                 where ns is either true or false (default = true)
% 
% name_ons(..., 'decol', decol)   defines the color of the estimated noise with decrease per frequency, 
%                                 where decol is either 'brown' or 'pink' (default = 'brown')
% 
% name_ons(..., 'debeg', debeg)   the decrease per frequency begins after the defined Hz (default = 1)
% 
% name_ons(..., 'incol', incol)   defines the color of the estimated noise with increase per frequency, 
%                                 where incol is either 'blue' or 'violet' (default = 'violet')
% 
% name_ons(..., 'lowfreq', lowfreq)       defines the low frequency extreme range in Hz for noise
%                                         intensity estimate as [minimum maximum] (default = [0 10])
% 
% name_ons(..., 'highfreq', highfreq)     defines the high frequency extreme range in Hz for noise
%                                         intensity estimate as [minimum maximum] (default = [20000 22050])
% 
% name_ons(..., 'dynbuffer', dynbuffer)   defines added dB white noise as noise floor buffer
%                                         to account for variance in noise over time (default = 0)
% 
% name_ons(..., 'mono', mono)     defines whether noise is estimated for the average TFRs across channels,
%                                 where mono is either true or false (default = true)
% 
% 
% Optional arguments for vibrato and tremolo suppression:
% 
% name_ons(..., 'timesmooth', timesmooth)  defines time-smoothing with lowpass at defined maximum Hz
%                                          to suppress tremolo (default = 20) (cf. Lartillot & Grandjean 2019)
% 
% name_ons(..., 'int', int)       defines minimum intensity increase threshold in dB
%                                 to suppress vibrato/tremolo (default = 1)
% 
% name_ons(..., 'soa', soa)       defines minimum sound onset asynchrony in ms (default = 50)
% 
% name_ons(..., 'dur', dur)       defines minimum sound onset duration in ms (default = 0)
% 
% 
% Optional arguments for the output: 
% 
% name_ons(..., 'outpath', outpath)       output path (default is same as the path with the input audio)
%                                         ( e.g., outpath = 'C:\folder\' )
% 
% name_ons(..., 'overwrite', overwrite)   overwrite any existing files without asking user (true or false)
%                                         (default is false) (e.g., use true for batch processing with no interruptions)
% 
% name_ons(..., 'visualize', visualize)   visualization of the sound onset analysis (true or false) (default is false)
%
% name_ons(..., 'audiosave', audiosave)   save audio file with onsets marked by 50 ms 1000 Hz sine tones (true or false) 
%                                         (default is false)
% 
% name_ons(..., 'outsuffix', outsuffix)   text string added to ending of output files (e.g., for testing specific settings)
%                                         (default is '') ( e.g., outsuffix = '_onsets_auto' )
% 
% name_ons(..., 'echange', echange)       save spectral energy change analyses (true or false) (default is false)
% 
% [onsets, ons] = name_ons(audio)         also stores settings and noise estimates in ons
% 
% 
% Beta version 20230607. Developed at Center for Music in the Brain by Niels Trusbak Haumann. https://musicinthebrain.au.dk/ 
% 
% The onset detection with noise suppression (ONS) is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 

% Parse and check the input arguments
p = inputParser;
addOptional(p, 'framedur', 100) % Time frame duration in milliseconds for the TFR analysis (default = 100)
addOptional(p, 'srtfr', 100) % Sampling rate in Hz for the TFR (default = 100)
addOptional(p, 'ns', true) % Apply audio medium noise suppression, true or false (default = true)
addOptional(p, 'decol', 'brown') % Color of noise with decrease per frequency ('brown' or 'pink') (default = 'brown')
addOptional(p, 'debeg', 1) % For noise with decrease per frequency, decrease begins after defined Hz (default = 1)
addOptional(p, 'incol', 'violet') % Color of noise with decrease per frequency ('blue' or 'violet') (default = 'violet')
addOptional(p, 'lowfreq', [0 10]) % Frequency range in Hz for low frequency extreme ([minimum maximum])
addOptional(p, 'highfreq', [20000 22050]) % Frequency range in Hz for high frequency extreme ([minimum maximum])
addOptional(p, 'dynbuffer', 0) % Added defined dB white noise as noise floor buffer to account for variance in noise over time (default = 0)
addOptional(p, 'mono', true) % (by default the noise is estimated for the average TFRs across channels)
addOptional(p, 'timesmooth', 20) % Time smooting to suppress tremolo (Hz) (lowpass at defined maximum Hz) (cf. Lartillot & Grandjean 2019) (default = 20)
addOptional(p, 'int', 1) % Minimum intensity increase in dB for onset detection above vibrato/tremolo level (default = 1)
addOptional(p, 'soa', 50) % Minimum SOA (sound onset-asynchrony) in ms for note onset detection (default = 50)
addOptional(p, 'dur', 0) % Minimum duration of energy increase in ms until end of increase or new increase (default = 0)
[inputpath, ~, ~] = fileparts(audio);
addOptional(p, 'outpath', [inputpath,filesep]) % Output path (default is same as the path with the input audio)
addOptional(p, 'outsuffix', '') % Text string added to ending of output files (e.g., for testing specific settings) (default is '')
addOptional(p, 'echange', false) % Save spectral energy change analyses (true or false) (default is false)
addOptional(p, 'visualize', false) % Visualization of onset intensities and detected onset times on spectral change curve (true or false) (default is false)
addOptional(p, 'audiosave', false) % Save audio file with onsets marked by 50 ms 1000 Hz sine tones (true or false) (default is false)
addOptional(p, 'overwrite', false) % (by default overwriting existing files is false)
parse(p, varargin{:})
if ~ischar(audio)
    error('audio variable must be a character array (string).')
end
framedur = p.Results.framedur; % Time frame duration in milliseconds for the TFR analysis (default = 100)
if ~isnumeric(framedur)  || length(framedur) ~= 1
    error('framedur variable must be a positive value. The default is 100.')
end
if framedur <= 0
    error('framedur variable must be a positive value. The default is 100.')
end
srtfr = p.Results.srtfr; % Sampling rate in Hz for the TFR (default = 100)
if ~isnumeric(srtfr)  || length(srtfr) ~= 1
    error('srtfr variable must be a positive value. The default is 100.')
end
if srtfr <= 0
    error('srtfr variable must be a positive value. The default is 100.')
end
ns = p.Results.ns; % Apply audio medium noise suppression, true or false (default = true)
if ~islogical(ns) || length(ns) ~= 1
    error('ns variable must be either true or false. The default is true.')
end
overwrite = p.Results.overwrite; 
if ~islogical(overwrite) || length(overwrite) ~= 1
    error('overwrite variable must be either true or false. The default is false.')
end
% (The audio medium noise suppression variables are checked in the name_noise_audmed subfunction.)
decol = p.Results.decol;
debeg = p.Results.debeg;
incol = p.Results.incol;
lowfreq = p.Results.lowfreq;
highfreq = p.Results.highfreq;
dynbuffer = p.Results.dynbuffer;
mono = p.Results.mono;
timesmooth = p.Results.timesmooth; % Time smooting to suppress tremolo (Hz) (lowpass at defined maximum Hz) (cf. Lartillot & Grandjean 2019) (default = 20)
if ~isnumeric(timesmooth)  || length(timesmooth) ~= 1
    error('timesmooth variable must be a positive value. The default is 20.')
end
if timesmooth <= 0
    error('timesmooth variable must be a positive value. The default is 20.')
end
int = p.Results.int; % Minimum intensity increase in dB for onset detection above vibrato/tremolo level (default = 1)
if ~isnumeric(int)  || length(int) ~= 1
    error('int variable must be a positive value. The default is 1.')
end
if int <= 0
    error('int variable must be a positive value. The default is 1.')
end
soa = p.Results.soa; % Minimum SOA (sound onset-asynchrony) in ms for note onset detection (default = 50)
if ~isnumeric(soa)  || length(soa) ~= 1
    error('soa variable must be a positive value. The default is 50.')
end
if soa <= 0
    error('soa variable must be a positive value. The default is 50.')
end
dur = p.Results.dur; % Minimum duration of energy increase in ms until end of increase or new increase (default = 0)
if ~isnumeric(dur)  || length(dur) ~= 1
    error('dur variable must be a positive value. The default is 0.')
end
if dur < 0
    error('dur variable must be a positive value. The default is 0.')
end
outpath = p.Results.outpath; % Output path (default is same as the path with the input audio)
if ~ischar(outpath)
    error('outpath variable must be a character array (string).')
end
outsuffix = p.Results.outsuffix; % Text string added to ending of output files (e.g., for testing specific settings) (default is '')
if ~ischar(outsuffix)
    error('outsuffix variable must be a character array (string).')
end
echange = p.Results.echange; % Save spectral energy change analyses (true or false) (default is false)
if ~islogical(echange) || length(echange) ~= 1
    error('echange variable must be either true or false. The default is false.')
end
visualize = p.Results.visualize; % Visualization of onset intensities and detected onset times on spectral change curve (true or false) (default is false)
if ~islogical(visualize) || length(visualize) ~= 1
    error('visualize variable must be either true or false. The default is false.')
end
audiosave = p.Results.audiosave; % Save audio file with onsets marked by 50 ms 1000 Hz sine tones (true or false) (default is false)
if ~islogical(audiosave) || length(audiosave) ~= 1
    error('audiosave variable must be either true or false. The default is false.')
end


%% Load the audio file

disp(['Loading the audio file ''',audio,'''.'])
[~, outname, outtype] = fileparts(audio);
[ inputWave, srtime ] = audioread( audio );


%%  Time-frequency analysis

disp('Time-frequency analysis with vibrato suppression.')
fftsize = 2*round(framedur/1000*srtime/2); % FFT frame size
% Calculate the TFR
[tfr_a, ~, ~, srtfr, freq, ~] = name_tfr(inputWave, srtime, 'frametime','center', 'srtfr',srtfr, 'framesize',fftsize);


%% Suppress audio medium noise spectrum

if ns
    disp('Suppressing audio medium noise...')
    [tfr_a, noise_audmed] = name_noise_audmed(tfr_a, freq, 'overwrite',overwrite , 'decol',decol, 'debeg', debeg, 'incol', incol, 'lowfreq', lowfreq, 'highfreq', highfreq, 'dynbuffer', dynbuffer, 'mono', mono, 'visualize', visualize);
end


%% Convert the TFR to mono version

disp('Converting the TFR to mono version.')
tfr_a = mean(tfr_a,3); % Convert the TFR to mono version


%% Use time-smoothing to counteract tremolo (cf. Lartillot & Grandjean 2019)

disp('Suppressing tremolo.')
tsmoothsamples = round( ( srtfr / timesmooth)/2 ); % Find half cycle samples for the defined maximum Hz
shiftsize = round(srtime/srtfr); % Time frame shift size in relation to audio time samples
nframes = floor( size(inputWave,1)/shiftsize); % Number of time frames
tfr_a_smooth = zeros(size(tfr_a));
for j=1:nframes % Loop over time frames
    
    % Use the mean of the preceding half cycle time samples for the defined maximum Hz
    if j < tsmoothsamples % (If it's the beginning with less than half cycle available, add zeros first)
        tfr_a_smooth(:,j) = mean( horzcat( zeros(size(tfr_a,1),tsmoothsamples-j), tfr_a(:,1:j) ) , 2 );
    else
        tfr_a_smooth(:,j) = mean( tfr_a(:,j-(tsmoothsamples-1):j) , 2 );
    end
end
tfr_a = tfr_a_smooth; clear('tfr_a_smooth');


%% Visualize the TFR

if visualize
    figure('color','w'), imagesc(20*log10(tfr_a)); set(gca,'ydir','normal'), set(gca,'clim',[-120 max(max(20*log10(tfr_a)))])
    title(['Time frequency representation [',outname,']'],'interpreter','none')
    set(gca,'xtick', 1:0.1*srtfr:size(tfr_a,2))
    set(gca,'xticklabel', 0:0.1:(size(tfr_a,2)-1)/srtfr), xlabel('Time (seconds)')
    set(gca,'ytick', 1:round(1000/(1/(fftsize/srtime))):size(tfr_a,1))
    set(gca,'yticklabel', 0:1000:srtime/2), ylabel('Frequency (Hz)')
end


%% Sound onset detection based on the sum of spectral energy increases

% The analysis of spectra for sound onset detection ensures that
% amplitude increases in the time-domain waveform related to beating
% or roughness are ignored.

disp('Detecting sound onsets based on the sum of spectral energy increases.')

% Finding maximum signal gain for spectral change over variable attack durations
disp('   - Estimating the maximum signal gain for spectral change over variable attack durations.')
attacksamples = 1:round(srtfr*0.050); % Attack durations are tested for up to 50 ms attack duration
spectchangedb_attdurs = []; % Prepare spectral change curve for variable attack durations ( attack duration index , time sample )
spectchangedb_leg = {}; % Prepare figure legend for optional visualization
for j=attacksamples % Loop over attack durations
    spectchangedb_attdurs(j,:) = [zeros(1,j), 20*log10( sum( tfr_a(:,1+j:end) , 1 ) ./ sum( tfr_a(:,1:end-j) , 1 ) )]; % Summed spectral amplitude change across frequency bands in dB
    spectchangedb_leg{j} = ['Spectral change in dB per ',num2str(round(j/srtfr*1000)),' ms']; % Prepare figure legend for optional visualization
end
[spectchangedb_attdurs(size(spectchangedb_attdurs,1)+1,:), attacksamplesid] = max(spectchangedb_attdurs, [], 1); % Add a last entry with the maximum signal gain across attack durations calculated for each time sample...
spectchangedb_attdurs(end,:) = spectchangedb_attdurs(end,:) ./ attacksamples(attacksamplesid); % ... and convert it to dB per time resolution of the TFR
spectchangedb_leg{length(spectchangedb_leg)+1} = ['Maximum signal gain converted to dB per ',num2str(round(1/srtfr*1000)),' ms']; % Add the figure legend for optional visualization of the maximum signal gain estimate
time = 0:1/srtfr:(length(spectchangedb_attdurs)-1)/srtfr; % Calculate the time in seconds for optional visualization and saving of energy change curve
if visualize
    figure('color','w')
    plot(time, spectchangedb_attdurs')
    xlabel('Time s.')
    ylabel('Spectral change (dB)')
    legend(spectchangedb_leg)
    title(['Signal gain for variable attack durations in ''',outname,''''],'interpreter','none')
end
spectchangedb = spectchangedb_attdurs(end,:); % Use only the spectral change estimate based on the attack duration with optimal signal gain calculated for each time sample

% Detect local energy increases
eincr = []; % Energy increase ( increase index , 1 = start time sample and 2 = end time sample)
[~,incrpeaks] = findpeaks( spectchangedb , 'MINPEAKHEIGHT', 0 ); % Detect spectral increase peaks
zcross = find( spectchangedb == 0 | ([0 diff(sign(spectchangedb))] == 2) ); % Detect zero crossing in forwards direction
[~,incrvalleys] = findpeaks( -spectchangedb ); incrvalleys = incrvalleys( spectchangedb(incrvalleys)>0 ); % Detect spectral increase valleys
incrvalleys = unique([1 incrvalleys length(spectchangedb)]); % Add first and last time sample as increase valleys if they are not already detected
eincrcurve = nan(size(spectchangedb)); % Prepare to save the energy increase curve across time samples
for j=1:length(incrpeaks) % Loop over energy increase peaks
    
    % Detect the start time point of the energy increase as the nearest zero-crossing or valley above zero, which precedes or is at the peak in the energy increase
    eincr(j,1) = max( [ zcross( find( incrpeaks(j) - zcross >= 0,1,'last') ) , incrvalleys( find( incrpeaks(j) - incrvalleys >= 0,1,'last') ) ] ); % Energy increase ( increase index , 1 = start time sample )
    
    % Detect the end time point of the energy increase at nearest zero-crossing or valley above zero, which follows after or is at the peak in the energy increase
    eincr(j,2) = min( [ zcross( find(  zcross - incrpeaks(j) >= 0,1,'first') ) , incrvalleys( find( incrvalleys - incrpeaks(j) >= 0,1,'first') ) ] ); % Energy increase ( increase index , 2 = end time sample)
    
    % The energy increase curve is estimated by integrating the spectral change derivative from the start to the end time point
    eincrcurve( eincr(j,1):eincr(j,2) ) = cumsum( spectchangedb( eincr(j,1):eincr(j,2) ) ) - spectchangedb(eincr(j,1)); % Energy increase curve across time samples
    
end
disp(['   - ',num2str(length(eincr)),' local energy increases were found.'])


%% Apply onset detection constraints

% Apply minimum intensity increase per time frame above noise floor
intrem = [];
for j=1:size(eincr,1) % Loop over the detected energy increases
    
    % If the maximum intensity in the energy increase curve is less than the intensity increase treshold, 
    % then assume it's noise and mark the detection to be removed
    if max( eincrcurve( eincr(j,1):eincr(j,2) ) ) < int 
        intrem = vertcat(intrem, j);
    end
end
for j=1:length(intrem) % Loop over the energy increases marked to be removed
    eincrcurve( eincr(intrem(j),1):eincr(intrem(j),2) ) = NaN; % Remove the energy increase curve
end
eincr(intrem,:) = []; % Remove the energy increases marked to be removed
disp(['   - ',num2str(length(eincr)),' local energy increases were retained above ',num2str(int),' dB intensity threshold.'])

% Apply SOA constraint
ei = []; % Prepare to store the maximum intensity for each energy increase
for j=1:size(eincr,1)
    ei(j) = max( eincrcurve( eincr(j,1):eincr(j,2) ) ); % Maximum intensity for each energy increase
end
[~,eiid] = sort(ei,'descend'); % Sort energy increases by their maximum intensity
soarem = [];
for j=1:length(eiid) % Loop over the energy increases, from the highest to the lowest intensity
    if ~(ismember(j,soarem)) % If not it's a weaker energy increase that is already marked to be removed with the SOA constraint...
        % Remove any following energy increases that occur faster than the SOA treshold, assuming that it's noise
        soarem = unique(vertcat(soarem, find( ( eincr(:,1) - eincr(eiid(j),1) < soa/1000*srtfr )  &  ( eincr(:,1) - eincr(eiid(j),1) > 0 ) ))); 
    end
end
for j=1:length(soarem) % Loop over the energy increases marked to be removed
    eincrcurve( eincr(soarem(j),1):eincr(soarem(j),2) ) = NaN; % Remove the energy increase curve
end
eincr(soarem,:) = []; % Remove the energy increases marked to be removed
disp(['   - ',num2str(length(eincr)),' local energy increases were retained above ',num2str(soa),' ms SOA threshold.'])

% Apply energy increase duration constraint
durrem = [];
for j=1:size(eincr,1)-1 % Loop over the energy changes, except the last one
    
    % If the start time between the energy increase and the following energy increase is shorter than the duration treshold, assume it's noise and remove the detection
    if eincr(j+1,1) - eincr(j,1) < dur/1000*srtfr 
        durrem = vertcat(durrem, j); 
    end
end
for j=1:length(durrem) % Loop over the energy increases marked to be removed
    eincrcurve( eincr(durrem(j),1):eincr(durrem(j),2) ) = NaN; % Remove the energy increase curve
end
eincr(durrem,:) = []; % Remove the energy increases marked to be removed
disp(['   - ',num2str(length(eincr)),' local energy increases were retained above ',num2str(dur),' ms duration threshold.'])


%% Save the detected sound onsets

% Collect the detected sound onsets as the start time points of the energy increases in seconds
onsets = (eincr(:,1)-1)/srtfr;
disp('Completed detecting sound onsets.')

if visualize
    
    % Show the detected sound onsets on the energy change curve
    figure('color','w')
    plot(time, spectchangedb)
    hold on
    scatter(onsets, spectchangedb(eincr(:,1)), 'k')
    legend({['Spectral change per time frame (dB/',num2str(round(1000/srtfr)),' ms)'], 'Detected sound onset time points'})
    xlabel('Time (seconds)')
    ylabel('Spectral change (dB)')
    title(['Detected sound onset time points in ''',outname,''''],'interpreter','none')
    
    % Show the estimated energy increase curve
    figure('color','w')
    plot(time, eincrcurve, 'r')
    hold on
    plot(time, int,'k')
    legend({'Sound intensity curve (dB)','Sound intensity threshold (dB)'})
    xlabel('Time (seconds)')
    ylabel('Sound intensity (dB)')
    title(['Sound intensity increase for each detected sound onset in ''',outname,''''],'interpreter','none')
end

% Save an audio file with the detected sound onsets marked with audio feedback
if audiosave
    timesamples = (0:round(0.050*srtime)-1); % Time samples for 50 ms sine wave
    sinesamples = repmat( sin( 1000*2*pi * timesamples/srtime)', [1, size(inputWave,2)]); % 50 ms sine wave
    outputWave = inputWave; % Copy the input audio waveform to the output audio waveform
    for j=1:length(onsets) % Looop over detected sound onsets
        
        if round(onsets(j)*srtime) +1 + timesamples <= size(inputWave,1) % As long as the 50 ms audio feedback does not exceed the duration of the audio waveform
            
            amplitude = max(abs(inputWave(round(onsets(j)*srtime)+1 + timesamples))); % Adjust the amplitude of the onset marker to the sound intensity of the input audio
            amplitude = 10^(6/20)*amplitude; % Increase the onset marker intensity by 6 dB
            
            outputWave( round(onsets(j)*srtime)+1 + timesamples, : ) = inputWave( round(onsets(j)*srtime) +1 + timesamples , :) + amplitude*sinesamples; % Add the audio feedback to the output waveform
        else
            
            remsamples = round(onsets(j)*srtime):size(inputWave,1); % Use the remaining time samples in the audio waveform
            
            amplitude = max(abs(inputWave(remsamples))); % Adjust the amplitude of the onset marker to the sound intensity of the input audio
            amplitude = 10^(6/20)*amplitude; % Increase the onset marker intensity by 6 dB
            
            outputWave( remsamples, : ) = inputWave( remsamples, : ) + amplitude*sinesamples(1:length(remsamples), :); % Add the audio feedback to the remaining time samples in the output waveform
        end
        
    end
    if max(abs(outputWave)) > 1 % If audio will be clipped ( <-1 or >+1 ), lower the amplitude to the maximum tolerated ( -1 or +1 )
        outputWave = outputWave/max(abs(outputWave));
    end
    % If overwriting files is not permitted in the options and the output file already exists, 
    % ask the user to confirm overwriting or changing the output file
    if overwrite == false && exist(fullfile([outpath, outname, outsuffix, '_onsets', outtype]),'file') == 2
        warning(['File ',fullfile([outpath, outname, outsuffix, '_onsets', outtype]),' already exists.'])
        [outname, outpath] = uiputfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Save audio with onsets as...', [outpath, outname, outsuffix, '_onsets', outtype]);
        if outname==0
            error('Please select an output file.')
        end
        [outpath,outname,filetype] = fileparts(fullfile([outpath,outname])); % Reconstruct the output path, file name, and file type
        outpath = [outpath,filesep]; % Restore the system-specific sign between the path and file name
        outsuffix = ''; % Apply the user specied file name with no additional user predefined suffix
        suffix = ''; % Apply the user specied file name without the default suffix
    else
        suffix = '_onsets'; % Apply the default suffix
    end
    audiowrite( fullfile([outpath, outname, outsuffix, suffix, outtype]), outputWave, srtime );
    disp(['Audio with onsets were saved to ',fullfile([outpath, outname, outsuffix, suffix, outtype])])
end

ons = p.Results; % Save all the applied options in ons structure
ons.audio = audio; % Add the input audio path and file name
ons.noise_audmed = noise_audmed; % Add the audio medium noise amplitude estimates in dB (1 = Brownian (red), 2 = Pink, 3 = Blue, 4 = Violet, 5 = White)

% If chosen, save Matlab vector with the detected sound onsets, energy change analyses, and the settings
if echange
    ons.time = time; % Time in seconds
    ons.spectchangedb = spectchangedb; % Spectral change in dB
    ons.eincr = eincr; % Energy increase ( increase index , 1 = start time sample and 2 = end time sample)
    ons.eincrcurve = eincrcurve; % Estimated energy increase curves in dB over time
else
    % Else, only save Matlab vector with the detected sound onset time points in seconds and the settings
end

% If overwriting files is not permitted in the options and the output file already exists,
% ask the user to confirm overwriting or changing the output file
if overwrite == false && exist(fullfile([outpath, outname, outsuffix,'.mat']),'file') == 2
    warning(['File ',fullfile([outpath, outname, outsuffix]),' already exists.'])
    [outname, outpath] = uiputfile('*.mat', 'Save onsets as...', [outpath, outname, outsuffix,'.mat']);
    if outname==0
        error('Please select an output file.')
    end
    [outpath,outname,filetype] = fileparts(fullfile([outpath,outname])); % Reconstruct the output path, file name, and file type
    outpath = [outpath,filesep]; % Restore the system-specific sign between the path and file name
    outsuffix = ''; % Apply the user specied file name with no additional user predefined suffix
end
save(fullfile([outpath, outname, outsuffix]),'onsets','ons');
if ons.echange
    disp(['Onsets in seconds, energy change analyses, and detection settings were saved to ''',fullfile([outpath, outname, outsuffix]),'.mat''.'])
else
    disp(['Onsets in seconds and detection settings were saved to ''',fullfile([outpath, outname, outsuffix]),'.mat''.'])
end

end