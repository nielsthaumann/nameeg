function noise_audmed = name_ns(audio, varargin)

% Audio medium noise suppression
% 
% Use as: 
% 
% noise_audmed = name_ns(audio)  
% 
% to suppress noise in an audio file, where, e.g., audio = 'C:\folder\audio.wav'
% 
% The audio with noise suppression is saved to an output file 
% (by default with the same path and name as the input file added the ending ..._ns.(filetype))
% 
% 
% Optional arguments for time frequency representation (TFR):
% 
% name_ns(..., 'framedur', framedur)   specifies the time frame duration in milliseconds 
%                                       for the TFR analysis (default = 100)
% 
% name_ns(..., 'srtfr', srtfr)         specifies the sampling rate for the TFR in Hz (default = 100)
%
% 
% Optional arguments for audio medium noise suppression:
% 
% name_ns(..., 'decol', decol)   defines the color of the estimated noise with decrease per frequency, 
%                                where decol is either 'brown' or 'pink' (default = 'brown')
% 
% name_ns(..., 'debeg', debeg)   the decrease per frequency begins after the defined Hz (default = 1)
% 
% name_ns(..., 'incol', incol)   defines the color of the estimated noise with increase per frequency, 
%                                where incol is either 'blue' or 'violet' (default = 'violet')
% 
% name_ns(..., 'lowfreq', lowfreq)       defines the low frequency extreme range in Hz for noise
%                                        intensity estimate as [minimum maximum] (default = [0 10])
% 
% name_ns(..., 'highfreq', highfreq)     defines the high frequency extreme range in Hz for noise
%                                        intensity estimate as [minimum maximum] (default = [20000 22050])
% 
% name_ns(..., 'dynbuffer', dynbuffer)   defines added dB white noise as noise floor buffer
%                                        to account for variance in noise over time (default = 0)
% 
% name_ns(..., 'mono', mono)     defines whether noise is estimated for the average TFRs across channels,
%                                where mono is either true or false (default = true)
% 
% 
% Optional arguments for the output: 
% 
% name_ns(..., 'visualize', visualize)   visualization of the noise analysis (true or false) (default is false)
%
% name_ns(..., 'outpath', outpath)       output path (default is same as the path with the input audio)
%                                        ( e.g., outpath = 'C:\folder\' )
% 
% name_ns(..., 'overwrite', overwrite)   overwrite any existing files without asking user (true or false)
%                                        (default is false) (e.g., use true for batch processing with no interruptions)
% 
% name_ns(..., 'duration', duration)     save noise spectrum with defined duration in seconds for perceptual validation
%                                        (default is 0 for no saving of noise spectrum audio)
% 
% name_ns(..., 'bitdepth', bitdepth)     valid audio bit depths are 8, 16, 24, 32, or 64. The default is 16.
% 
% noise_audmed = name_ns(audio)          noise_audmed shows the audio medium noise amplitude estimates in dB 
%                                        (for 1 = Brownian (red), 2 = pink, 3 = blue, 4 = violet, and 5 = white noise)
% 
% 
% Beta version 20230607. Developed at Center for Music in the Brain by Niels Trusbak Haumann. https://musicinthebrain.au.dk/ 
% 
% The noise suppression (NS) is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 


% Parse and check the input arguments
p = inputParser;
addOptional(p, 'overwrite', false) % (by default overwriting existing files is false)
addOptional(p, 'framedur', 100) % Time frame duration in milliseconds for the TFR analysis (default = 100)
addOptional(p, 'srtfr', 100) % Sampling rate in Hz for the TFR (default = 100)
addOptional(p, 'decol', 'brown') % Color of noise with decrease per frequency ('brown' or 'pink') (default = 'brown')
addOptional(p, 'debeg', 1) % For noise with decrease per frequency, decrease begins after defined Hz (default = 1)
addOptional(p, 'incol', 'violet') % Color of noise with decrease per frequency ('blue' or 'violet') (default = 'violet')
addOptional(p, 'lowfreq', [0 10]) % Frequency range in Hz for low frequency extreme ([minimum maximum])
addOptional(p, 'highfreq', [20000 22050]) % Frequency range in Hz for high frequency extreme ([minimum maximum])
addOptional(p, 'dynbuffer', 0) % Added defined dB white noise as noise floor buffer to account for variance in noise over time (default = 0)
addOptional(p, 'mono', true) % (by default the noise is estimated for the average TFRs across channels)
addOptional(p, 'visualize', false) % (by default visualization of noise modeling is false)
addOptional(p, 'duration', 0) % Save noise spectrum with defined duration in seconds for perceptual validation (default is 0 for no saving of noise spectrum audio)
addOptional(p, 'tfr_p', []) % (by default TFR phase is not provided)
addOptional(p, 'win', []) % (by default window function is not provided)
addOptional(p, 'srtime', []) % (by default the sampling rate of the time series is not provided)
addOptional(p, 'bitdepth', 16) % (by default the audio bit depth is 16)
[inputpath, ~, ~] = fileparts(audio);
addOptional(p, 'outpath', [inputpath,filesep]) % Output path (default is same as the path with the input audio)
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
visualize = p.Results.visualize; % Visualize the noise modeling (false or true)
if ~islogical(visualize) || length(visualize) ~= 1
    error('visualize variable must be either true or false. The default is false.')
end
% (The audio medium noise suppression variables are checked in the name_noise_audmed subfunction.)
overwrite = p.Results.overwrite; 
decol = p.Results.decol;
debeg = p.Results.debeg;
incol = p.Results.incol;
lowfreq = p.Results.lowfreq;
highfreq = p.Results.highfreq;
dynbuffer = p.Results.dynbuffer;
mono = p.Results.mono;
duration = p.Results.duration; % Save noise spectrum with defined duration in seconds for perceptual validation (default is 0 seconds)
tfr_p = p.Results.tfr_p; % TFR phase in radians
win = p.Results.win; % Window function
srtime = p.Results.srtime; % Sampling rate of the time series
bitdepth = p.Results.bitdepth; % Audio bit depth
outpath = p.Results.outpath; % Output path (default is same as the path with the input audio)
if ~ischar(outpath)
    error('outpath variable must be a character array (string).')
end


%% Load the audio file

disp(['Loading the audio file ''',audio,'''.'])
[ inputWave, srtime ] = audioread( audio );


%%  Time-frequency analysis

disp('Time-frequency analysis with vibrato suppression.')
fftsize = 2*round(framedur/1000*srtime/2); % FFT frame size
[tfr_a, tfr_p, win, srtfr, freq, ~] = name_tfr(inputWave, srtime, 'frametime','center', 'srtfr',srtfr, 'framesize',fftsize);


%% Suppress audio medium noise and save audio with suppressed noise

disp('Suppressing audio medium noise and saving audio with supressed noise...')
[~, noise_audmed] = name_noise_audmed(tfr_a, freq, 'overwrite', overwrite, 'decol',decol, 'debeg', debeg, 'incol', incol, 'lowfreq', lowfreq, 'highfreq', highfreq, 'dynbuffer', dynbuffer, 'mono', mono, 'duration', duration, 'visualize', visualize, 'audio', audio, 'frametime','center', 'tfr_p', tfr_p, 'win', win, 'srtfr', srtfr, 'srtime',srtime, 'bitdepth', bitdepth, 'outpath',outpath);


end