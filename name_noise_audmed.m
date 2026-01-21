function [tfr_a_ns, noise_audmed] = name_noise_audmed(tfr_a, freq, varargin)

% Suppress noise from audio medium
% 
% Use as: 
% 
% [tfr_a_ns, noise_audmed] = name_noise_audmed(tfr_a, freq) 
% 
% Required inputs are: 
% 
% tfr_a         amplitude spectra over time in a numeric frequency by time by channel matrix 
%               with the frequnecy bins along the first dimension, time frames
%               along the second dimension, and channels along the third dimension
% 
% freq          frequencies in Hz
%
% 
% Optional arguments for audio medium noise suppression:
% 
% [...] = name_noise_audmed(..., 'decol', decol)   defines the color of the estimated noise with decrease per frequency, 
%                                                  where decol is either 'brown' or 'pink' (default = 'brown')
% 
% [...] = name_noise_audmed(..., 'debeg', debeg)   the decrease per frequency begins after the defined Hz (default = 1)
% 
% [...] = name_noise_audmed(..., 'incol', incol)   defines the color of the estimated noise with increase per frequency, 
%                                                  where incol is either 'blue' or 'violet' (default = 'violet')
% 
% [...] = name_noise_audmed(..., 'lowfreq', lowfreq)       defines the low frequency extreme range in Hz for noise
%                                                          intensity estimate as [minimum maximum] (default = [0 10])
% 
% [...] = name_noise_audmed(..., 'highfreq', highfreq)     defines the high frequency extreme range in Hz for noise
%                                                          intensity estimate as [minimum maximum] (default = [20000 22050])
% 
% [...] = name_noise_audmed(..., 'dynbuffer', dynbuffer)   defines added dB white noise as noise floor buffer
%                                                          to account for variance in noise over time (default = 0)
% 
% [...] = name_noise_audmed(..., 'mono', mono)             defines whether noise is estimated for the average TFRs across channels,
%                                                          where mono is either true or false (default = true)
% 
% [...] = name_noise_audmed(..., 'static', static)         by default suppression of the static part of noise is false (it can be set to true)
%                                                          (For sound onset detection the static noise spectrum is returned to avoid 
%                                                          infinitely large dB increases between 0 and the noise floor amplitude. 
%                                                          The static noise supression option only affects the tfr_a_ns output
%                                                          and not the audio output with noise suppression.)
% 
% Optional arguments for the output: 
% 
% [...] = name_noise_audmed(..., 'visualize', visualize)   visualization of the noise analysis (true or false) (default is false)
%
% [...] = name_noise_audmed(..., 'audio', audio)           Save audio with noise suppression to the defined audio path and file name, 
%                                                          e.g., audio = 'C:\folder\audio_ns.wav'
%                                                          (by default saving to audio file name is disabled)
% 
% [...] = name_noise_audmed(..., 'outpath', outpath)       output path (default is current folder)
%                                                          ( e.g., outpath = 'C:\folder\' )
% 
% [...] = name_noise_audmed(..., 'overwrite', overwrite)   overwrite any existing files without asking user (true or false)
%                                                          (default is false) (e.g., use true for batch processing with no interruptions)
% 
% [...] = name_noise_audmed(..., 'duration', duration)     save noise spectrum with defined duration in seconds for perceptual validation
%                                                          (default is 0 for no saving of noise spectrum audio)
% 
% [...] = name_noise_audmed(..., 'bitdepth', bitdepth)     valid audio bit depths are 8, 16, 24, 32, or 64. The default is 16.
% 
% 
% NB: When using the audio option the following input variables tfr_p, win, srtfr, srtime, and frametime are required: 
% 
% [...] = name_noise_audmed(..., 'tfr_p', tfr_p)             phase spectra over time in radians (between -pi and +pi) in a numeric frequency 
%                                                            by time by channel matrix with the frequency bins along the first dimension, 
%                                                            time frames along the second dimension, and channels along the third dimension
% 
% [...] = name_noise_audmed(..., 'win', win)                 Blackman-Harris window
% 
% [...] = name_noise_audmed(..., 'srtfr', srtfr)             sampling rate of the TFR in Hz
% 
% [...] = name_noise_audmed(..., 'srtime', srtime)           sampling rate of the time series in Hz
% 
% [...] = name_noise_audmed(..., 'frametime', frametime)     frame time reference, either 'center' or 'end'
% 
%  
% Outputs are: 
% 
% tfr_a_ns              amplitude spectra over time with noise suppression, 
%                       a numeric frequency by time by channel matrix 
%                       with the frequnecy bins along the first dimension, time frames
%                       along the second dimension, and channels along the third dimension
% 
% noise_audmed          noise_audmed shows the audio medium noise amplitude estimates in dB 
%                       (for 1 = Brownian (red), 2 = pink, 3 = blue, 4 = violet, and 5 = white noise)
% 
% Beta version 20230607. 
% 
% name_noise_audmed is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 


% Parse and check the input arguments
p = inputParser; 
addOptional(p, 'decol', 'brown') % Color of noise with decrease per frequency ('brown' or 'pink') (default = 'brown')
addOptional(p, 'debeg', 1) % For noise with decrease per frequency, decrease begins after defined Hz (default = 1 Hz)
addOptional(p, 'incol', 'violet') % Color of noise with decrease per frequency ('blue' or 'violet') (default = 'violet')
addOptional(p, 'lowfreq', [0 10]) % Frequency range in Hz for low frequency extreme ([minimum maximum])
addOptional(p, 'highfreq', [20000 22050]) % Frequency range in Hz for high frequency extreme ([minimum maximum])
addOptional(p, 'dynbuffer', 0) % Added defined dB white noise as noise floor buffer to account for variance in noise over time (default = 0)
addOptional(p, 'mono', true) % (by default the noise is estimated for the average TFRs across channels)
addOptional(p, 'static', false) % (by default suppression of static part of noise is false)
addOptional(p, 'visualize', false) % (by default visualization of noise modeling is false)
addOptional(p, 'audio', []) % (by default saving to audio file name is disabled)
addOptional(p, 'overwrite', false) % (by default overwriting existing files is false)
addOptional(p, 'duration', 0) % Save noise spectrum with defined duration in seconds for perceptual validation (default is 0 for no saving of noise spectrum audio)
addOptional(p, 'tfr_p', []) % (by default TFR phase is not provided)
addOptional(p, 'win', []) % (by default window function is not provided)
addOptional(p, 'srtfr', []) % (by default the sampling rate of the TFR is not provided)
addOptional(p, 'srtime', []) % (by default the sampling rate of the time series is not provided)
addOptional(p, 'frametime', []) % (by default the frame time reference is is not provided)
addOptional(p, 'bitdepth', 16) % (by default the audio bit depth is 16)
addOptional(p, 'outpath', pwd) % Output path (default is the current folder)
parse(p, varargin{:})
decol = p.Results.decol; % Color of noise with decrease per frequency ('brown' or 'pink') (default = 'brown')
if ~ischar(decol)
    error('decol variable must be either ''brown'' or ''pink''. The default is ''brown''.')
end
if ~ismember({'brown','pink'}, decol)
    error('decol variable must be either ''brown'' or ''pink''. The default is ''brown''.')
end
debeg = p.Results.debeg; % For noise with decrease per frequency, decrease begins after defined Hz (default = 1 Hz)
if ~isnumeric(debeg) || length(debeg) ~= 1
    error(['debeg variable must be a value within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is 1.'])
end
if debeg < freq(1) || debeg > freq(end)
    error(['debeg variable must be a value within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is 1.'])
end
incol = p.Results.incol; % Color of noise with decrease per frequency ('blue' or 'violet') (default = 'violet')
if ~ischar(incol)
    error('incol variable must be either ''blue'' or ''violet''. The default is ''violet''.')
end
if ~ismember({'blue','violet'}, incol)
    error('incol variable must be either ''blue'' or ''violet''.')
end
lowfreq = p.Results.lowfreq; % Frequency range in Hz for low frequency extreme ([minimum maximum])
if ~isnumeric(lowfreq)  || length(lowfreq) ~= 2
    error(['lowfreq variable must be provided as [min max] with values within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is [0 10].'])
end
if lowfreq(1) < freq(1) || lowfreq(1) > freq(end) || lowfreq(2) < freq(1) || lowfreq(2) > freq(end)
    error(['lowfreq variable must be provided as [min max] with values within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is [0 10].'])
end
highfreq = p.Results.highfreq; % Frequency range in Hz for high frequency extreme ([minimum maximum])
if ~isnumeric(highfreq)  || length(highfreq) ~= 2
    error(['highfreq variable must be provided as [min max] with values within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is [20000 22050].'])
end
if highfreq(1) < freq(1) || highfreq(1) > freq(end) || highfreq(2) < freq(1) || highfreq(2) > freq(end) || sum( size(highfreq) ~= [1, 2] )>0
    error(['highfreq variable must be provided as [min max] with values within the available frequency range between ',num2str(freq(1)),' and ',num2str(freq(end)),' Hz. The default is [20000 22050].'])
end
dynbuffer = p.Results.dynbuffer; % Added defined dB white noise as noise floor buffer to account for variance in noise over time (default = 0)
if ~isnumeric(dynbuffer) || length(dynbuffer) ~= 1
    error('dynbuffer variable must be a dB value. The default is 0.')
end
mono = p.Results.mono; % Estimate noise based on the average TFRs across channels (false or true)
if ~islogical(mono) || length(mono) ~= 1
    error('mono variable must be either true or false. The default is true.')
end
static_ns = p.Results.static; % Suppression of the static part of noise (false or true)
if ~islogical(static_ns) || length(static_ns) ~= 1
    error('static variable must be either true or false. The default is false.')
end
visualize = p.Results.visualize; % Visualize the noise modeling (false or true)
if ~islogical(visualize) || length(visualize) ~= 1
    error('visualize variable must be either true or false. The default is false.')
end
audiofile = p.Results.audio; % Save to audio file name (empty [] or string)
if ~ischar(audiofile) && ~isempty(audiofile)
    error('The variable audio must be provided as a string. E.g., ''C:\folder\audio_ns.wav''')
end
overwrite = p.Results.overwrite; % (by default overwriting existing files is false)
duration = p.Results.duration; % Save noise spectrum with defined duration in seconds for perceptual validation (default is 10 seconds)
tfr_p = p.Results.tfr_p; % TFR phase in radians
win = p.Results.win; % Window function
srtfr = p.Results.srtfr; % Sampling rate of the TFR
srtime = p.Results.srtime; % Sampling rate of the time series
frametime = p.Results.frametime; % Frame time reference
bitdepth = p.Results.bitdepth; % Audio bit depth
outpath = p.Results.outpath; % Output path (default is current folder)
if ~isempty(audiofile) % Verify that the required variables are provided when using the audio option...
    if ~islogical(overwrite) || length(overwrite) ~= 1
        error('overwrite variable must be either true or false. The default is false.')
    end
    if ~isnumeric(duration) || length(duration) ~= 1
        error('The variable duration must be provided as a positive value when using the audio option. The default is 10.')
    end
    if duration < 0
        error('The variable duration must be provided as a positive value when using the audio option. The default is 0.')
    end
    if ~isnumeric(tfr_p)
        error('The variable tfr_p of same size as tfr_a must be provided as a numeric matrix when using the audio option.')
    end
    if ~isequal(size(tfr_p), size(tfr_a))
        error('The variable tfr_p of same size as tfr_a must be provided as a numeric matrix when using the audio option.')
    end
    if ~isnumeric(win)
        error(['The variable win of size [',num2str(2*(size(tfr_a,1)-1)),', 1] must be provided as a numeric array when using the audio option.'])
    end
    if ~isequal(size(win), [2*(size(tfr_a,1)-1) 1])
        error(['The variable win of size [',num2str(2*(size(tfr_a,1)-1)),', 1] must be provided as a numeric array when using the audio option.'])
    end
    if ~isnumeric(srtfr) || length(srtfr) ~= 1
        error('The variable srtfr must be provided as a value when using the audio option.')
    end
    if ~isnumeric(srtime) || length(srtime) ~= 1
        error('The variable srtime must be provided as a value when using the audio option.')
    end
    if ~ischar(frametime)
        error('The variable frametime must be provided as either ''center'' or ''end'' when using the audio option.')
    end
    if ~ismember({'center','end'}, frametime)
        error('The variable frametime must be provided as either ''center'' or ''end'' when using the audio option.')
    end
    [~,~,filetype] = fileparts(audiofile);
    if strcmpi(filetype,'.flac')
        if ~ismember([8,16,24], bitdepth)
            error('Valid audio bit depths for FLAC audio are 8, 16, or 24. The default is 16.')
        end
    else
        if ~ismember([8,16,24,32,64], bitdepth)
            error('Valid audio bit depths are 8, 16, 24, 32, or 64. The default is 16.')
        end
    end
    if ~ischar(outpath)
        error('outpath variable must be a character array (string).')
    end
end


%% Simulate Brownian, pink, blue, violet, and white noise spectra

disp('Simulating noise spectra.')
noisespect = []; % Simulated audio medium noise amplitude spectra ( frequency bin , 1 = Brownian (red) - 2 = Pink - 3 = Blue - 4 = Violet - 5 = White)
colors = []; % Show colors for optional visualization
% #1 Brownian (red) noise
beta = 1; % Brownian noise beta = 1 for amplitude spectrum (=2 for power spectrum)
noisespect(:,1) = 1./(freq.^beta); % Amplitude per frequency bin
noisespect( freq < debeg , 1) = 1/(debeg^beta); % Constrain beginning of decrease until defined Hz
colors(:,1) =  [1, 0, 0]; % Red color
% #2 Pink noise
beta = 0.5; % Pink noise beta = 0.5 for amplitude spectrum (=1 for power spectrum)
noisespect(:,2) = 1./(freq.^beta); % Amplitude per frequency bin
noisespect( freq < debeg , 2) = 1/(debeg^beta); % Constrain beginning of decrease until defined Hz
colors(:,2) =  [0.9600, 0.7600, 0.7600]; % Pink color
% #3 Blue noise
beta = 0.5; % Blue noise beta = 0.5 for amplitude spectrum (=1 for power spectrum)
noisespect(:,3) = freq.^beta; % Amplitude per frequency bin
colors(:,3) =  [0, 0, 1]; % Blue color
% #4 Violet noise
beta = 1; % Violet noise beta = 1 for amplitude spectrum (=2 for power spectrum)
noisespect(:,4) = freq.^beta; % Amplitude per frequency bin
colors(:,4) =  [0.4196, 0.2980, 0.6039]; % Violet color
% #5 White noise
beta = 0; % White noise beta = 0 for amplitude spectrum (=0 for power spectrum)
noisespect(:,5) = freq.^beta; % Amplitude per frequency bin
colors(:,5) =  [1, 1, 1]; % White color
if visualize
    figure('color','k')
    hold on
    for j=[4,3,5,2,1]
        plot(freq, 20*log10(noisespect(:,j)) , 'color', colors(:,j))
    end
    set(gca,'Color','k')
    set(gca,'xscale','log')
    title('Simulated noise spectra','color','w')
    legend({'Violet', 'Blue', 'White', 'Pink', 'Brownian'}, 'TextColor','w','location','SouthWest')
    xlabel('Frequency (Hz)','color','w')
    ylabel('dB','color','w')
    set(gca,'xcolor','w')
    set(gca,'ycolor','w')
end


%% Estimate initial fit of the simulated noise spectra to the complete average audio spectrum across time

if mono
    disp('Estimating noise based on a mono average across channels.')
else
    disp('Estimating noise for each channel separately.')
end

disp('Estimating initital fit of the simulated noise spectra to the complete audio spectrum.')
noisecmask = zeros(length(freq),5); % Prepare noise color noise mask
if strcmp(decol,'brown')
    disp(['Using Brownian (red) noise with offset at ',num2str(debeg),' Hz for fitting.'])
    noisecmask(:,1) = 1; % Pass Brownian (red) noise in noise mask
elseif strcmp(decol,'pink')
    disp('Using pink noise with offset at ',num2str(debeg),' Hz for fitting.')
    noisecmask(:,2) = 1; % Pass pink noise in noise mask
end
if strcmp(incol,'blue')
    disp('Using blue noise for fitting.')
    noisecmask(:,3) = 1; % Pass blue noise in noise mask
elseif strcmp(incol,'violet')
    disp('Using violet noise for fitting.')
    noisecmask(:,4) = 1; % Pass violet noise in noise mask
end
noisecmask(:,5) = 1; % Pass white noise in noise mask

if mono
    % Estimate the noise amplitude with non-negative least squares multiple regression for the mean across channels
    noise_amplitude = lsqnonneg( noisecmask .* noisespect , mean(mean(tfr_a,3),2) , optimset('TolX',eps) )'; 
    
else % multichannel
    noise_amplitude = []; 
    for c=1:size(tfr_a,3) % Loop over channels
        % Estimate the noise amplitude with non-negative least squares multiple regression for each channel separately
        noise_amplitude(c,:) = lsqnonneg( noisecmask .* noisespect , mean(tfr_a(:,:,c),2) , optimset('TolX',eps) )'; 
    end
end


%%  Adjust the estimated noise spectrum to the low and high frequency extremes

disp(['Adjusting the estimated noise spectrum to the low (',num2str(lowfreq(1)),'-',num2str(lowfreq(2)),' Hz) and high (',num2str(highfreq(1)),'-',num2str(highfreq(2)),' Hz) frequency extremes.'])

freqbins_low = find(freq >= lowfreq(1) & freq <= lowfreq(2)); % Find the low frequency extreme bins
freqbins_high = find(freq >= highfreq(1) & freq <= highfreq(2)); % Find the high frequency extreme bins

if mono
    [~,freqbins_lowmaxid] = max( mean(mean( tfr_a(freqbins_low,:,:) , 3),2) , [], 1 ); % Find the maximum audio amplitude in relation to the low frequency extreme bins
    freqbins_lowmaxid = freqbins_low(freqbins_lowmaxid); % Find the frequency bin with maximum audio amplitude at the low frequency extreme
    [~,freqbins_highmaxid] = max( mean(mean( tfr_a(freqbins_high,:,:) , 3),2) , [], 1 ); % Find the maximum audio amplitude in relation to the high frequency extreme bins
    freqbins_highmaxid = freqbins_high(freqbins_highmaxid); % Find the frequency bin with maximum audio amplitude at the high frequency extreme
    sum_colored_noise = sum( noisespect([freqbins_lowmaxid, freqbins_highmaxid], 1:4) .* repmat(noise_amplitude(1:4), [2, 1]) , 2 );  % Sum the colored noise at the found maxima bins
    b = lsqnonneg( [ sum_colored_noise , ones(size(sum_colored_noise)) ] , mean(mean(tfr_a([freqbins_lowmaxid, freqbins_highmaxid], :,:),3),2) ,optimset('TolX',eps) )'; % Adjust the noise amplitude with non-negative least squares multiple regression for the mean across channels
    noise_amplitude(1:4) = b(1) * noise_amplitude(1:4); % Adjust the colored part of the noise spectrum
    noise_amplitude(5) = b(2); % Adjust the white noise part of the noise spectrum
    noisespect_est = sum( repmat(noise_amplitude, [length(freq), 1]) .* noisespect , 2) * 10^(dynbuffer/20); % Estimated noise spectrum
    disp('   - Estimated noise amplitude per frequency band (f Hz): ')
    disp(['      - Brownian noise = 20*log10(   1/f   ) ',num2str(round(20*log10(noise_amplitude(1)))),' dB'])
    disp(['      - Pink noise     = 20*log10( 1/f^0.5 ) ',num2str(round(20*log10(noise_amplitude(2)))),' dB'])
    disp(['      - Blue noise     = 20*log10(  f^0.5  ) ',num2str(round(20*log10(noise_amplitude(3)))),' dB'])
    disp(['      - Violet noise   = 20*log10(    f    ) ',num2str(round(20*log10(noise_amplitude(4)))),' dB'])
    disp(['      - White noise    =                     ',num2str(round(20*log10(noise_amplitude(5)))),' dB'])
    
else % multichannel
    [~,freqbins_lowmaxid] = max( squeeze(mean(tfr_a(freqbins_low,:,:),2)) , [], 1 ); % Find the maximum audio amplitude in relation to low frequency extreme bins
    freqbins_lowmaxid = freqbins_low(freqbins_lowmaxid); % Find the frequency bin with the maximum audio amplitude at the low frequency extreme
    [~,freqbins_highmaxid] = max( squeeze(mean(tfr_a(freqbins_high,:,:),2)) , [], 1 ); % Find the maximum audio amplitude in relation to the high frequency extreme bins
    freqbins_highmaxid = freqbins_high(freqbins_highmaxid); % Fidn the frequency bin with the maximum audio amplitude at the high frequency extreme
    sum_colored_noise = [];
    noisespect_est = []; 
    for c=1:size(tfr_a,3)
        sum_colored_noise(c,:) = sum( noisespect([freqbins_lowmaxid(c), freqbins_highmaxid(c)], 1:4) .* repmat(noise_amplitude(c,1:4), [2, 1]) , 2 );  % Sum of colored noise at the found maximum bins
        b = lsqnonneg( [ sum_colored_noise(c,:)' , ones(size(sum_colored_noise(c,:)))' ] , mean(tfr_a([freqbins_lowmaxid(c), freqbins_highmaxid(c)], :,c),2) ,optimset('TolX',eps) )'; % Adjust the noise amplitude with non-negative least squares multiple regression for each channel seperately
        noise_amplitude(c,1:4) = b(1) * noise_amplitude(c,1:4); % Adjust the colored part of the noise spectrum
        noise_amplitude(c,5) = b(2); % Adjust the white noise part of the spectrum
        noisespect_est(c,:) = sum( repmat(noise_amplitude(c,:), [length(freq), 1]) .* noisespect , 2) * 10^(dynbuffer/20); % Estimated noise spectrum
        disp(['   - Estimated noise amplitude per frequency band (f Hz) in channel ',num2str(c),': '])
        disp(['      - Brownian noise = 20*log10(   1/f   ) ',num2str(round(20*log10(noise_amplitude(c,1)))),' dB'])
        disp(['      - Pink noise     = 20*log10( 1/f^0.5 ) ',num2str(round(20*log10(noise_amplitude(c,2)))),' dB'])
        disp(['      - Blue noise     = 20*log10(  f^0.5  ) ',num2str(round(20*log10(noise_amplitude(c,3)))),' dB'])
        disp(['      - Violet noise   = 20*log10(    f    ) ',num2str(round(20*log10(noise_amplitude(c,4)))),' dB'])
        disp(['      - White noise    =                     ',num2str(round(20*log10(noise_amplitude(c,5)))),' dB'])
    end
end
if visualize
    if mono
        figure('color','w')
        plot(freq, 20*log10(mean(mean(tfr_a,3),2))) % Mono
        set(gca,'xscale','log')
        hold on
        plot(freq, 20*log10( noisespect_est ) ,'k') % Mono
        set(gca,'xscale','log')
        title('Simulated noise spectrum fitted to audio spectrum (mono)','interpreter','none')
        xlabel('Frequency (Hz)')
        ylabel('dB')
        grid on
        legend({'Audio spectrum','Fitted noise spectrum (threshold)'})
        ylims = get(gca,'ylim'); ylim([ylims(1) 0])
    else % multichannel
        for c=1:size(tfr_a,3)
            figure('color','w')
            plot(freq, 20*log10(mean(tfr_a(:,:,c),2))) % Multichannel
            set(gca,'xscale','log')
            hold on
            plot(freq, 20*log10( noisespect_est(c,:) ) ,'k') % Multichannel
            set(gca,'xscale','log')
            title(['Simulated noise spectrum fitted to audio spectrum (channel ',num2str(c),' of ',num2str(size(tfr_a,3)),')'],'interpreter','none')
            xlabel('Frequency (Hz)')
            ylabel('dB')
            grid on
            legend({'Audio spectrum','Fitted noise spectrum (threshold)'})
            ylims = get(gca,'ylim'); ylim([ylims(1) 0])
        end
    end
end


%% Suppress the estimated audio medium noise

disp('Suppressing the estimated dynamic and static noise of the audio medium in the TFR.')
if mono
    tfr_a_ns = double( tfr_a .*  ( tfr_a > repmat( noisespect_est , [1, size(tfr_a,2), size(tfr_a,3)])) ); % Suppress the dynamic and static noise spectrum
else % multichannel
    tfr_a_ns = zeros(size(tfr_a)); 
    for c=1:size(tfr_a,3) % Loop over channels
        tfr_a_ns(:,:,c) = double( tfr_a(:,:,c) .*  ( tfr_a(:,:,c) > repmat( noisespect_est(c,:)' , [1, size(tfr_a,2)])) ); % Suppress the dynamic and static noise spectrum
    end
end
if ~static_ns
    disp('Returning the estimated static noise of the audio medium to the TFR.')
    if mono
        tfr_a_ns = tfr_a_ns + repmat( noisespect_est , [1, size(tfr_a,2), size(tfr_a,3)]); % Return the noise spectrum as static version
    else % multichannel
        for c=1:size(tfr_a,3) % Loop over channels
            tfr_a_ns(:,:,c) = tfr_a_ns(:,:,c) + repmat( noisespect_est(c,:)' , [1, size(tfr_a,2)]); % Return the noise spectrum as static version
        end
    end
end
disp('Completed suppressing audio medium noise in the TFR.')

if ~isempty(audiofile)
    
    % Save audio with suppressed noise
    disp(' ')
    [~,filename,filetype] = fileparts(audiofile);
    % If overwriting files is not permitted in the options and the output file already exists, 
    % ask the user to confirm overwriting or changing the output file
    if overwrite == false && exist(fullfile([outpath, filesep, filename,'_ns', filetype]),'file') == 2
        warning(['File ',fullfile([outpath, filesep, filename,'_ns', filetype]),' already exists.'])
        [filename, outpath] = uiputfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Save audio with noise suppression as...', [outpath,filename,'_ns',filetype]);
        if filename==0
            error('Please select an output file.')
        end
        [outpath,filename,filetype] = fileparts(fullfile([outpath,filename])); % Reconstruct the output path, file name, and file type
        suffix = ''; % Apply the user specied file name without the default suffix
    else
        suffix = '_ns'; % Apply the default suffix
    end
    disp(['Saving audio with noise suppression to ''',fullfile([outpath, filesep, filename, suffix, filetype]),'''...'])
    if mono
        outputSpectrogram = tfr_a - repmat( noisespect_est , [1, size(tfr_a,2), size(tfr_a,3)]); % Suppress noise by subtraction of noise spectrum
    else % multichannel
        outputSpectrogram = zeros(size(tfr_a));
        for c=1:size(tfr_a,3) % Loop over channels
            outputSpectrogram(:,:,c) = tfr_a(:,:,c) - repmat( noisespect_est(c,:)' , [1, size(tfr_a,2)]); % Suppress noise by subtraction of noise spectrum
        end
    end
    outputSpectrogram( outputSpectrogram<0 ) = 0; % Ensure all values are non-negative
    outputWave = name_invtfr(outputSpectrogram, tfr_p, win, srtfr, srtime, frametime); % Apply inverse TFR
    disp(['Adding white noise dither at ',num2str(20*log10(1/2^(bitdepth-1))), ' dB to mask quantization errors in audio sampled at ',num2str(bitdepth),' bit depth.'])
    outputWave = outputWave + 1/2^(bitdepth-1)*repmat(2*rand(size(outputWave,1),1)-1, [1, size(outputWave,2)]); % Add white noise dither at relevant bit depth amplitude limit
    if strcmpi(filetype,'.wav') || strcmpi(filetype,'.flac')
        audiowrite( fullfile([outpath, filesep, filename, suffix, filetype]) , outputWave , srtime, 'BitsPerSample', bitdepth)
    else % Else if is lossy format that does not support definition of BitsPerSample
        audiowrite( fullfile([outpath, filesep, filename, suffix, filetype]) , outputWave , srtime)
    end
    disp('Completed saving audio with noise suppression.')
    
    if duration >0
        % If chosen, save defined seconds audio with estimated noise
        disp(' ')
        [~,filename,filetype] = fileparts(audiofile);
        % If overwriting files is not permitted in the options and the output file already exists,
        % ask the user to confirm overwriting or changing the output file
        if overwrite == false && exist(fullfile([outpath, filesep, filename,'_noise', filetype]),'file') == 2
            warning(['File ',fullfile([outpath, filesep, filename,'_noise', filetype]),' already exists.'])
            [filename, outpath] = uiputfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Save noise as...', [outpath,filename,'_noise',filetype]);
            if filename==0
                error('Please select an output file.')
            end
            [outpath,filename,filetype] = fileparts(fullfile([outpath,filename])); % Reconstruct the output path, file name, and file type
            suffix = ''; % Apply the user specied file name without the default suffix
        else
            suffix = '_noise'; % Apply the default suffix
        end
        disp(['Saving ',num2str(duration),' seconds audio with estimated noise to ''',fullfile([outpath, filesep, filename,suffix, filetype]),'''...'])
        if mono
            outputSpectrogram = repmat( noisespect_est , [1, duration*srtfr]); % Estimated noise spectrum
        else % mutichannel
            outputSpectrogram = zeros(size(tfr_a,1), duration*srtfr, size(tfr_a,3));
            for c=1:size(tfr_a,3) % Loop over channels
                outputSpectrogram(:,:,c) = repmat( noisespect_est(c,:)' , [1, duration*srtfr]); % Estimated noise spectrum
            end
        end
        outputSpectrogram( outputSpectrogram<0 ) = 0; % Ensure all values are non-negative
        randomangle = repmat( 2*pi*rand( size(outputSpectrogram,1), size(outputSpectrogram,2) )-pi , [1, 1, size(outputSpectrogram,3)]); % Create random angles mathcing the size of the output spectrogram
        outputWave = name_invtfr(outputSpectrogram, randomangle, win, srtfr, srtime, frametime); % Apply inverse TFR
        disp(['Adding white noise dither at ',num2str(20*log10(1/2^(bitdepth-1))), ' dB to mask quantization errors in audio sampled at ',num2str(bitdepth),' bit depth.'])
        outputWave = outputWave + 1/2^(bitdepth-1)*repmat(2*rand(size(outputWave,1),1)-1, [1, size(outputWave,2)]); % Add white noise dither at relevant bit depth amplitude limit
        if strcmpi(filetype,'.wav') || strcmpi(filetype,'.flac')
            audiowrite( fullfile([outpath, filesep, filename, suffix, filetype]) , outputWave , srtime, 'BitsPerSample', bitdepth)
        else % Else if is lossy format that does not support definition of BitsPerSample
            audiowrite( fullfile([outpath, filesep, filename, suffix, filetype]) , outputWave , srtime)
        end
    end
    
end

noise_audmed = 20*log10(noise_amplitude); % Output the audio medium noise amplitude estimates in dB (1 = Brownian (red), 2 = Pink, 3 = Blue, 4 = Violet, 5 = White)

end