function [onsets, intensity, ons] = name_ons(audio, varargin)

% Automatic sound onset detection with noise suppression
% 
% Use as: 
% 
% onsets = name_ons(audio)  
% 
% to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'
% 
% The onsets are given as an array of time points in seconds.  
% 
% The algorithm includes suppression of noise in the audio medium. 
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
% name_ons(..., 'tsm', tsm)       defines whether time-smoothing (17 Hz lowpass filter) should be applied
%                                 to suppress vibrato/tremolo (Lartillot & Grandjean 2019), 
%                                 where tsm is either true or false (default = true)
% 
% name_ons(..., 'cthr', cthr)     defines minimum intensity increase threshold 
%                                 to suppress vibrato/tremolo in normalized values between 0-1 
%                                 (default = 0.01)
% 
% name_ons(..., 'thr', thr)       defines minimum intensity threshold 
%                                 to suppress background noise in normalized values between 0-1 
%                                 (default = 0)
% 
% 
% Optional arguments for the output: 
% 
% name_ons(..., 'onsave', onsave)         save onsets time points to MAT-file (true or false) 
%                                         (default is false) (e.g., use true for batch processing with no interruptions)
% 
% name_ons(..., 'outpath', outpath)       output path (default is same as the path with the input audio)
%                                         ( e.g., outpath = 'C:\folder' )
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
% [onsets, intensity, ons] = name_ons(audio)    also stores the peak intensity (dB) at each sound onset 
%                                               and the settings and noise estimates in ons
% 
% 
% Please make sure that a recent version of the MIRtoolbox is installed. 
% (This function was tested with MIRtoolbox 1.8.1.)
% E.g, visit: https://github.com/olivierlar/mirtoolbox
% 
% Please cite this paper when applying name_ons, which calls MIRtoolbox functions: 
% Lartillot, O., & Toiviainen, P. (2007). A Matlab Toolbox for Musical
% Feature Extraction From Audio. International Conference on Digital
% Audio Effects. 
% 
% Please also cite this paper when using the time-smooting (tsm) option: 
% Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A., &
% Grandjean, D. (2013). Estimating tempo and metrical features by
% tracking the whole metrical hierarchy. Paper presented at the The 3rd
% International Conference on Music & Emotion, Jyväskylä, Finland,
% June 11-15, 2013. 
% 
% Beta version 20250605. Developed at Center for Music in the Brain by Niels Trusbak Haumann. https://musicinthebrain.au.dk/ 
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
addOptional(p, 'tsm', true) % Time smooting to suppress tremolo, true or false (Lartillot & Grandjean 2019) (default = true)
addOptional(p, 'cthr', 0.01) % Minimum intensity increase threshold to suppress vibrato/tremolo in normalized values between 0-1 (default = 0.01)
addOptional(p, 'thr', 0) % Minimum intensity threshold to suppress background noise in normalized values between 0-1 (default = 0)
addOptional(p, 'onsave', false) % (by default saving onsets to MAT file is false)
[inputpath, ~, ~] = fileparts(audio);
addOptional(p, 'outpath', inputpath) % Output path (default is same as the path with the input audio)
addOptional(p, 'outsuffix', '') % Text string added to ending of output files (e.g., for testing specific settings) (default is '')
addOptional(p, 'visualize', false) % Visualization of onset intensities and detected onset times on spectral change curve (true or false) (default is false)
addOptional(p, 'audiosave', false) % Save audio file with onsets marked by 50 ms 1000 Hz sine tones (true or false) (default is false)
addOptional(p, 'overwrite', false) % (by default overwriting existing files is false)
parse(p, varargin{:})
if ~ischar(audio)
    error('audio variable must be a character array (string).')
end
[outpath, outname, outtype] = fileparts(audio); % Store audio filname and filetype
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
onsave = p.Results.onsave; 
if ~islogical(onsave) || length(onsave) ~= 1
    error('onsave variable must be either true or false. The default is false.')
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
tsm = p.Results.tsm; % Time smooting to suppress tremolo, true or false (Lartillot & Grandjean 2019) (default = true)
if ~islogical(tsm) || length(tsm) ~= 1
    error('tsm variable must be either true or false. The default is true.')
end
cthr = p.Results.cthr; % Minimum intensity increase threshold to suppress vibrato/tremolo in normalized values between 0-1 (default = 0.01)
if ~isnumeric(cthr)  || length(cthr) ~= 1
    error('cthr variable must be a positive value between 0 and 1. The default is 0.01.')
end
if cthr < 0 || cthr > 1 
    error('cthr variable must be a positive value between 0 and 1. The default is 0.01.')
end
thr = p.Results.thr; % Minimum intensity threshold to suppress background noise in normalized values between 0-1 (default = 0)
if ~isnumeric(thr)  || length(thr) ~= 1
    error('thr variable must be a positive value between 0 and 1. The default is 0.01.')
end
if thr < 0 || thr > 1 
    error('thr variable must be a positive value between 0 and 1. The default is 0.01.')
end
outpath = p.Results.outpath; % Output path (default is same as the path with the input audio)
if ~ischar(outpath)
    error('outpath variable must be a character array (string).')
end
outsuffix = p.Results.outsuffix; % Text string added to ending of output files (e.g., for testing specific settings) (default is '')
if ~ischar(outsuffix)
    error('outsuffix variable must be a character array (string).')
end
visualize = p.Results.visualize; % Visualization of onset intensities and detected onset times on spectral change curve (true or false) (default is false)
if ~islogical(visualize) || length(visualize) ~= 1
    error('visualize variable must be either true or false. The default is false.')
end
audiosave = p.Results.audiosave; % Save audio file with onsets marked by 50 ms 1000 Hz sine tones (true or false) (default is false)
if ~islogical(audiosave) || length(audiosave) ~= 1
    error('audiosave variable must be either true or false. The default is false.')
end


%% Suppress audio medium noise spectrum

if ns
    if overwrite == false && exist([outpath, filesep, outname,'_ns',outtype],'file')
        % Skip if not allowed to overwrite existing files and the noise suppressed audio file already exists
        
        disp(['Using already existing audio file with noise suppression: ',outpath, filesep, outname,'_ns',outtype])
        
    else
        disp('Suppressing audio medium noise...')
        noise_audmed = name_ns(audio, 'framedur', framedur, 'srtfr', srtfr, 'decol',decol, 'debeg', debeg, 'incol', incol, 'lowfreq', lowfreq, 'highfreq', highfreq, 'dynbuffer', dynbuffer, 'mono', mono, 'visualize', visualize, 'outpath', outpath, 'overwrite', overwrite);
    end
end


%% Sound onset detection using finetuned functions from the Music Information Retrieval (MIR) toolbox

if isempty(which('miraudio')) || isempty(which('mirspectrum')) || isempty(which('mirflux')) || isempty(which('mirevents')) 
    error('Please make sure that a recent version of the MIRtoolbox (Lartillot and Toiviai) is installed. (This function was tested with MIRtoolbox 1.8.1.) E.g, visit: https://github.com/olivierlar/mirtoolbox')
end

disp('Detecting sound onsets using finetuned functions from the Music Information Retrieval (MIR) toolbox.')
if ns
    a = miraudio([outpath, filesep, outname,'_ns',outtype]);
else
    a = miraudio(audio);
end


disp('Obtaining short-time Fourier transform (STFT) magnitude spectra with Terhardt model...')
% For features based on short-time Fourier spectra,
% use the Blackmann-Harris window function to suppress scalloping loss.
framedurs = framedur/1000; % Time frame duration conversion from milliseconds to seconds (the effective frame duration will be shorter due to the window function)
feature_spectrum = mirspectrum( a , 'Window', 'blackmanharris', 'Terhardt', 'Frame', framedurs, 's', 1/srtfr, 's' );

if tsm

    disp('Estimating spectral flux using the ''emerge'' (or ''smoothgate'') option...')
    feature_spectral_flux = mirflux(feature_spectrum, 'Inc','BackSmooth','Lartillot','Dist','Gate', 'Frame', framedurs, 's', 1/srtfr, 's');

else
    
    disp('Estimating spectral flux...')
    feature_spectral_flux = mirflux(feature_spectrum, 'Frame', framedurs, 's', 1/srtfr, 's');
    
end

disp('Detecting sound onsets...')
o = mirevents(feature_spectral_flux, 'Contrast', cthr, 'Threshold', thr);
onsets = mirgetdata(o) + framedurs/2;

% Estimating attacks
v = mirpeaks(feature_spectral_flux, 'Valleys', 'Contrast', 0.01, 'Threshold', 1); % Detect valleys in the spectral flux curve...
valleys = sort(mirgetdata(v) + framedurs/2);
attacks = zeros(length(onsets), 2); % Initialize storing attack time points (onset#, begin and end)
attacks(:,1) = onsets;
for i=1:length(onsets)
    attacks(i,2) = valleys( find(valleys > onsets(i), 1, 'first') ); % Identify the attack end as the first valley after the attack begin
end

disp(['   -> ',num2str(length(onsets)),' sound onsets were detected.'])

clear('a','feature_spectrum','feature_spectral_flux','o'); % Cleanup memroy


%% Save the detected sound onsets

% Store the peak intensity (dB) at each sound onset
[ inputWave, srtime ] = audioread(audio); % Read input audio waveform
intensity = zeros(size(onsets)); 
for i=1:length(onsets)
    
    timesamples = round(attacks(i,1)*srtime):round(attacks(i,2)*srtime); % Attack time samples in the audio
    timesamples( timesamples < 1 | timesamples > size(inputWave,1) ) = []; % Ensure the onset time samples are between the first and last audio sample
    intensity(i,1) = 20*log10( max(max(abs( inputWave( timesamples , :) ))) ); % Maximum intensity (dB) across onset time samples and channels

end

if visualize
    
    disp('Visualizing the results...')

    % Show the detected sound onsets
    time = 0:1/srtime:(size(inputWave,1)-1)/srtime;
    figure('color','w')
    plot(time, inputWave)
    hold on
    plot(time, NaN*ones(size(inputWave)),'k')
    for i=1:length(onsets)
        rectangle('Position',[attacks(i,1), -10.^(intensity(i)/20), attacks(i,2)-attacks(i,1), 2*10.^(intensity(i)/20)])
    end
    drawnow
    legend({'Audio waveform','Detected sound onsets (attack range)'})
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    title(['Detected sound onset time points in ''',outname,''''],'interpreter','none')
    
end

% Save an audio file with the detected sound onsets marked with audio feedback
if audiosave
    
    timesamples = (0:round(0.050*srtime)-1); % Time samples for 50 ms sine wave
    sinesamples = repmat( sin( 1000*2*pi * timesamples/srtime)', [1, size(inputWave,2)]); % 50 ms sine wave
    outputWave = inputWave; % Copy the input audio waveform to the output audio waveform
    for i=1:length(onsets) % Looop over detected sound onsets
        
        if round(onsets(i)*srtime) +1 + timesamples <= size(inputWave,1) % As long as the 50 ms audio feedback does not exceed the duration of the audio waveform
            
            amplitude = max(abs(inputWave(round(onsets(i)*srtime)+1 + timesamples))); % Adjust the amplitude of the onset marker to the sound intensity of the input audio
            amplitude = 10^(6/20)*amplitude; % Increase the onset marker intensity by 6 dB
            
            outputWave( round(onsets(i)*srtime)+1 + timesamples, : ) = inputWave( round(onsets(i)*srtime) +1 + timesamples , :) + amplitude*sinesamples; % Add the audio feedback to the output waveform
        else
            
            remsamples = round(onsets(i)*srtime):size(inputWave,1); % Use the remaining time samples in the audio waveform
            
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
    if overwrite == false && exist([outpath, filesep, outname, outsuffix, '_onsets', outtype],'file') == 2
        warning(['File ',[outpath, filesep, outname, outsuffix, '_onsets', outtype],' already exists.'])
        [outname, outpath] = uiputfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Save audio with onsets as...', [outpath, outname, outsuffix, '_onsets', outtype]);
        if outname==0
            error('Please select an output file.')
        end
        [outpath,outname,filetype] = fileparts([outpath, filesep, outname]); % Reconstruct the output path, file name, and file type
        outsuffix = ''; % Apply the user specied file name with no additional user predefined suffix
        suffix = ''; % Apply the user specied file name without the default suffix
    else
        suffix = '_onsets'; % Apply the default suffix
    end
    audiowrite( [outpath, filesep, outname, outsuffix, suffix, outtype], outputWave, srtime );
    disp(['Audio with onsets were saved to ',[outpath, filesep, outname, outsuffix, suffix, outtype]])

end

ons = p.Results; % Save all the applied options in ons structure
ons.audio = audio; % Add the input audio path and file name
if ns & exist('noise_audmed','var')
    ons.noise_audmed = noise_audmed; % Add the audio medium noise amplitude estimates in dB (1 = Brownian (red), 2 = Pink, 3 = Blue, 4 = Violet, 5 = White)
end

if onsave % If chosen, save Matlab vector with the detected sound onsets, energy change analyses, and the settings

    % If overwriting files is not permitted in the options and the output file already exists,
    % ask the user to confirm overwriting or changing the output file
    if overwrite == false && exist([outpath, outname, outsuffix,'.mat'],'file') == 2
        warning(['File ',[outpath, outname, outsuffix],' already exists.'])
        [outname, outpath] = uiputfile('*.mat', 'Save onsets as...', [outpath, outname, outsuffix,'.mat']);
        if outname==0
            error('Please select an output file.')
        end
        [outpath,outname,filetype] = fileparts([outpath, filesep, outname]); % Reconstruct the output path, file name, and file type
        outsuffix = ''; % Apply the user specied file name with no additional user predefined suffix
    end
    save([outpath, filesep, outname, outsuffix],'onsets','ons');
    disp(['Onsets in seconds and detection settings were saved to ''',[outpath, filesep, outname, outsuffix],'.mat''.'])
    
end


end
