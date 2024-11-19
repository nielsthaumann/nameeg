function name_onsed

% Sound onset editor. 
% 
% Use this editor to combine the speed of automatic sound onset detection 
% with the reliability of manual editing.
% 
% Load an audio file (.wav, .flac, .mp3, .m4a, .mp4, .ogg) for auditory validation
% and visual inspection of energy increases in a time-frequency representation. 
% 
% Load and save a .mat file with the variable 'onsets' containing a numerical vector 
% with the onsets stored in seconds. 
% 
% Changes are saved automatically (along with the undo/redo history). 
% 
% The default keyboard and mouse configuration is: 
% 
% Select a time point ------- Mouse click on time-frequency representation
% 
% Playback/Stop ------------- Space key ------- Mouse click on Play/Stop button
% Seek backwards on/off ----- Left arrow key -- Mouse click on << button
% Seek forwards on/off ------ Right arrow key - Mouse click on >> button
% Playback speed decrease --- Page Down key --- Mouse click on x1/? button
% Playback speed increase --- Page Up key ----- Mouse click on x1/? button
% Onsets feedback on/off ---- o key ----------- Mouse click on Onsets button
% 
% Previous onset ------------ Home key -------- Mouse click on Prev. button
% Next onset ---------------- End key --------- Mouse click on Next button
% Add onset ----------------- Insert key ------ Mouse click on Add button
% Delete onset -------------- Delete key ------ Mouse click on Delete button
% Delete visible onsets ----- d key ----------- Mouse click on Del.vis. button
% 
% Undo ---------------------- z key ----------- Mouse click on Undo button
% Redo ---------------------- x key ----------- Mouse click on Redo button
% 
% Time zoom out ------------- - key ----------- Mouse wheel scroll down or click on ? s. button
% Time zoom in -------------- + key ----------- Mouse wheel scroll up or click on ? s. button
% 
% (The keyboard configuration is stored in 'name_onsed_keys.mat' and can be modified if relevant.)
% 
% Beta version 20230530. Developed at Center for Music in the Brain by Niels Trusbak Haumann. https://musicinthebrain.au.dk/ 
% 
% The sound onset editor is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 


%% Load/save onsets and audio

% Load the keyboard configuration
name_onsed_keys = struct; 
load('name_onsed_keys.mat')

% Load saved or create new onsets
selection = questdlg('Edit onsets', '', 'Create new', 'Load saved', 'Load saved');
if strcmp(selection, 'Load saved')
    
    [onsetsfile, onsetspath] = uigetfile('*.mat', 'Load saved onsets...');
    if onsetsfile==0
        error('Please select an onset file.')
    end
    disp(['Loading the onsets in ''',onsetsfile,'''.'])
    load( fullfile(onsetspath, onsetsfile) );
    if ~exist('onsets','var')
        error('Please select an onset file.')
    end
    
    if exist('onsetsedit', 'var') % Resume previous editing if 'onsetsedit' variable was loaded
        audiopath = onsetsedit.audiopath;
        audiofile = onsetsedit.audiofile;
        
        if exist( fullfile(audiopath, audiofile) , 'file' )==2 % Ensure audio file exists in previous location
            
            disp(['Welome back! Resuming editing of ''',onsetsfile,'''.'])
            
        else
            
            warning(['Could not find ',fullfile(audiopath, audiofile)])
            
            [audiofile, audiopath] = uigetfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Load audio file...');
            if audiofile==0
                error('Please select an audio file.')
            end
            
        end
        
    else
        
        [audiofile, audiopath] = uigetfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Load audio file...');
        if audiofile==0
            error('Please select an audio file.')
        end
        
    end
    
elseif strcmp(selection, 'Create new')
    
    onsets = []; 
    
    [audiofile, audiopath] = uigetfile({'*.wav';'*.flac';'*.mp3';'*.m4a';'*.mp4';'*.ogg'}, 'Load audio file...');
    if audiofile==0
        error('Please select an audio file.')
    end
    
end

% Load the audio file
disp(['Loading the audio in ''',audiofile,'''.'])
[inputWave, fs] = audioread( fullfile(audiopath, audiofile) );

% Select the output file name for saving the edited onsets
if exist('onsetsedit', 'var')
    outputpath = onsetsedit.outputpath; 
    outputfile = onsetsedit.outputfile; 
else
    [outputfile, outputpath] = uiputfile('*.mat', 'Save onsets as...');
    if outputfile==0
       error('Please select an output file.')
    end
end
output = fullfile(outputpath, outputfile); 
disp(['Auto-saving the edited onsets to ''',output,'''.'])


%%  Time-frequency analysis

spectrumfile = fullfile([outputpath,audiofile,'tfr']);

if exist(spectrumfile,'file')~=2 % If previously created TFR data is not found in the output path
    
    disp('Preparing the time-frequency analysis with vibrato suppression.')
    ons = struct;
    ons.fftdur  = 100; % FFT window duration (milliseconds) (recommended = 100 ms)
    ons.tres    = 100; % Desired output time resolution or sampling rate (Hz) (recommended = 100 Hz)
    ons.freqrange = [10 22050]; % Frequencies in Hz are within the range [min max] (recommended is [10 22050] Hz)
    ons.fftSize = 2*round(ons.fftdur/1000*fs/2); % FFT size
    ons.shiftSize = round(fs/ons.tres); % Time frame shift size
    ons.tres = 1/(ons.shiftSize/fs); % Time resolution
    if license('test','Signal_Toolbox')
        bhwindow = blackmanharris(ons.fftSize); % Blackman-Harris window (NB: Requires the Signal Processing Toolbox)
    else
        n = (0:ons.fftSize-1); 
        bhwindow = ( 0.35875 - 0.48829*cos(2*pi*n/(ons.fftSize-1)) + 0.14128*cos(4*pi*n/(ons.fftSize-1)) - 0.01168*cos(6*pi*n/(ons.fftSize-1)) )'; % Blackman-Harris window
    end
    ons.frequencies = (fs*(0:floor(ons.fftSize/2))/ons.fftSize)'; % Frequencies in Hz
    ons.fs_spectrum = fs/ons.shiftSize; % Time-frequency representation sampling rate (Hz)
    
    % Prepare matrix indexing for parallel Fourier transform across time windows
    frameaddsamples = (0:ons.fftSize-1)'; % Index offset for FFT time window
    frameoffsets = (1:ons.shiftSize:size(inputWave,1))'; % Index offset for each time frame
    index = repmat(frameaddsamples, [1, length(frameoffsets)]) + repmat(frameoffsets', [length(frameaddsamples), 1]); % Matrix indexing (FFT window time samples, time frame)
    
    % Prepare zero-padded input wave with centered frame (NB: Timing differs from onset detection TFR)
    inputWavezp = vertcat(zeros(ons.fftSize/2-1, size(inputWave,2)), inputWave, zeros(ons.fftSize/2-1, size(inputWave,2)));
    
    % Prepare the TFR matrix
    inputSpectrogram = zeros(floor(ons.fftSize/2)+1, length(frameoffsets), size(inputWave,2)); % (frequency, time, channel)
    for c=1:size(inputWave,2) % Loop over channels
        
        disp(['Performing fast Fourier transform of ',num2str(length(ons.frequencies)),' frequency bands in ',num2str(length(frameoffsets)),' time frames in channel ',num2str(c), ' of ',num2str(size(inputWave,2)),'.'])
        
        % Parallel FFT calcuation across time windows
        specttemp = fft( repmat(bhwindow(:), [1, length(frameoffsets)]) .* reshape( inputWavezp(index, c) , [length(frameaddsamples), length(frameoffsets)]) );
        
        % Obtain the single-sided amplitude spectra
        specttemp = abs(specttemp/ons.fftSize);
        inputSpectrogram(:,:,c) = specttemp(1:floor(ons.fftSize/2)+1, :);
        inputSpectrogram(2:end-1,:,c) = 2*inputSpectrogram(2:end-1,:,c);
        
    end
    clear('inputWavezp'); % Cleanup memory
    
    % Convert the TFR to logarithmic scale mono version
    disp('Converting the time-frequency representation to a logarithmic scale mono version.')
    frequencybins = find( ons.frequencies >= ons.freqrange(1) & ons.frequencies <= ons.freqrange(2) ); % Find the frequency bins within the relevant range
    logscale = round(logspace( log10(frequencybins(1)), log10(frequencybins(end)) , length(frequencybins) )); % Convert the frequency bin scale from linear to logarithmic spacing
    logscale( logscale < 1 ) = 1; % Ensure the first bin is number 1 (not rounded to 0)
    logscale( logscale > length(frequencybins) ) = length(frequencybins); % Ensure the last log-spaced bin number corresponds to the last linear-spaced bin number
    logscaleinterp = logspace( log10(frequencybins(1)), log10(frequencybins(end)) , length(frequencybins) ); % Convert the frequency bin scale from linear to logarithmic spacing with interpolation for Gaussian kernel frequency smoothing (includes non-integers)
    
    % Create a gain matrix for converting from linear to log frequency scale
    lintologfreq = zeros(length(frequencybins), length(frequencybins)); % Gain matrix (output frequency bin, input frequency bin)
    for f=1:length(frequencybins) % Loop over the output frequency bins
        if f==1 % If is first frequency bin
            bins = logscale(f+1) - logscale(f); % Input frequency bins are the linearly spaced bins between the logarithmic spaced bins 1 and the middle between bins 1 and 2
              lintologfreq(f, logscale(1) : logscale(2) - ceil(1/2*((logscale(2)-logscale(1))))  ) = 1; % The output gain is the sum over the input frequency bins
              
          elseif f==length(frequencybins) % If is last frequency bin
              bins = logscale(f) - logscale(f-1); % Input frequency bins are the linearly spaced bins between the logarithmic spaced middle between the next last and last bins and the last bin
              lintologfreq(f, logscale(end) - floor(1/2*((logscale(end)-logscale(end-1)))) : logscale(end)  ) = 1; % The output gain is the sum over the input frequency bins
              
        else
              bins = logscale(f+1) - logscale(f-1); % Input frequency bins are the linearly spaced bins between the logarithmic spaced middle between the previous and current frequency bin and the current and next frequency bin
              lintologfreq(f, logscale(f) - floor(1/2*((logscale(f)-logscale(f-1)))) : logscale(f+1) - ceil(1/2*((logscale(f+1)-logscale(f))))  ) = 1; % The output gain is the sum over the input frequency bins
              
        end
        if bins <= 1 % Use a Gaussian kernel to interpolate the low frequency bands
            lintologfreq(f, :) = exp(-(( (1:length(frequencybins)) - logscaleinterp(f) )/(0.5)).^2); % Extend a Gaussian kernel centred on the interpolated input frequency bin (includes non-integers)
        end
    end
    lintologfreq = lintologfreq ./ repmat( sum(lintologfreq, 2), [1, length(frequencybins)] ); % Calculate the mean frequency bin weighting based on the summed frequency bin weighting
    
    % Multiply the mean TFR across audio channels by the linear to log frequency transfer matrix, and convert the intensity into dB scale
    inputSpectrogram_scaled = 20*log10( mean(inputSpectrogram(frequencybins,:,:),3)' * lintologfreq' )'; 
    clear('InputSpectrogram'); % Cleanup memory
    inputSpectrogram_scaled_max = max(max(inputSpectrogram_scaled)); % Store the global maximum intensity for adjusting visualization of limited viewing range
    
    disp(['Saving time-frequency analysis with vibrato suppression for later use to ''',audiofile,'tfr','''.'])
    save(spectrumfile, 'ons', 'inputSpectrogram_scaled', 'inputSpectrogram_scaled_max', 'logscale', '-v7.3')
    
else
    
    disp(['Loading a previous time-frequency analysis from ''',audiofile,'tfr','''.'])
    load(spectrumfile, '-mat')
    
end


%% GUI

disp('Starting the onset editor.')

% Initialize GUI
guifig = figure('Pointer', 'ibeam', 'MenuBar','none', 'name',['Sound onset editor [ ',outputfile,' <~> ',audiofile,' ]'], 'NumberTitle','off'); 

% Prepare GUI data
datagui = struct;
datagui.name_onsed_keys = name_onsed_keys; % Keyboard key assignments
datagui.timezoom = 5; % Time zoom in seconds
if exist('onsetsedit', 'var') % If resuming previous editing
    datagui.timesample = onsetsedit.timesample; % Last edited time sample
    datagui.undo = onsetsedit.undo; % Actions for undoing
    datagui.redo = onsetsedit.redo; % Actions for redoing
else % Prepare new editing
    datagui.timesample = 1; % Start time sample
    datagui.undo = []; % Actions for undoing
    datagui.redo = []; % Actions for redoing
end
datagui.onsets = onsets(:); % Onset time points in seconds
datagui.backwards = false; % Fast backwards state
datagui.forwards = false; % Fast forwards state
datagui.outputpath = outputpath; % Path to output onsets file
datagui.outputfile = outputfile; % Name of output onsets file
datagui.output = output; % Full path and name of output onsets file
datagui.close = false; % Closing state
datagui.audiopath = audiopath; % Path to audio file
datagui.audiofile = audiofile; % Name of audio file
guidata(guifig, datagui); % Save the GUI data

% Save the initial onsets and edits to the output file
if exist('onsetsedit', 'var')
    save(datagui.output, 'onsets', 'onsetsedit')
else
    onsetsedit = struct;
    onsetsedit.audiopath = audiopath;
    onsetsedit.audiofile = audiofile;
    datagui.outputpath = outputpath;
    datagui.outputfile = outputfile;
    onsetsedit.timesample = 1; 
    onsetsedit.undo = []; 
    onsetsedit.redo = []; 
    save(datagui.output, 'onsets', 'onsetsedit')
end

% Initialize audioplayer object
datagui.gain = 1; % Gain applied to avoid clipping when adding onsets audio feedback to audio
datagui.inputWave = inputWave; % Store audio waveform as GUI data
clear('inputWave'); % Cleanup memory
datagui.feedback = zeros(size(datagui.inputWave)); % Initialize onsets audio feedback toggled off
datagui.fs = fs; % Store audio sampling rate
datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.fs); % Create the audioplayer object
guidata(guifig, datagui); % Store the GUI data

% Create onsets feedback audio
timesamples = (0:round(0.050*datagui.fs)-1); % 50 ms feedback time samples
sinesamples = repmat( sin( 1000*2*pi * timesamples/fs)', [1, size(datagui.inputWave,2)]); % 1000 Hz sine wave feedback
for i=1:length(datagui.onsets) % Loop over onsets
    amplitude = max( abs( datagui.inputWave(round(datagui.onsets(i)*fs)+1 + timesamples) )); % Adjust amplitude of onset marker to sound intensity of input
    amplitude = 10^(6/20)*amplitude; % Increase onset marker by 6 dB
    datagui.feedback( round(datagui.onsets(i)*datagui.fs)+1 + timesamples , : ) = amplitude*sinesamples; % Add the onset audio
end
datagui.feedback = datagui.feedback(1:size(datagui.inputWave,1),:); % Trim to duration of input audio
if max(max(abs( datagui.inputWave + datagui.feedback ))) <= 1 % Check that amplitude is within maximum tolerated
    datagui.gain = 1; % Amplitude gain parameter
else % If audio will be clipped, lower gain to maximum tolerated
    datagui.gain = 1/max(max(abs( datagui.inputWave + datagui.feedback ))); % Amplitude gain parameter
end
datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.aobj.SampleRate); % Update the audioplayer object
datagui.feedbackstatus = true; % Set audio feedback state to true
guidata(guifig, datagui); % Store the GUI data

% Show the TFR
screenupdatefun(guifig,[], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
try
    set(guifig, 'SizeChangedFcn',@(src,event)screenupdatefun(src, event, ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)); 
catch
    set(guifig, 'ResizeFcn',@(src,event)screenupdatefun(src, event, ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)); 
end

% Initialize backwards button object
datagui.backbutton = uicontrol('Style','Togglebutton', 'String','<<', 'Position', [20, 20, 50, 30], 'Callback', @(src,event)backbuttonfun(src,event)); 

% Initialize playback button object
datagui.playbutton = uicontrol('Style','Pushbutton', 'String','Play', 'Position', [70, 20, 50, 30], 'Callback', @(src,event)playbuttonfun(src,event)); 

% Initialize forwards button object
datagui.forwardbutton = uicontrol('Style','Togglebutton', 'String','>>', 'Position', [120, 20, 50, 30], 'Callback', @(src,event)forwardbuttonfun(src,event)); 

% Initialize speed button object
datagui.speedbutton = uicontrol('Style','Popupmenu', 'String',{'x1/1','x1/2','x1/4','x1/8'}, 'Position', [170, 20, 50, 30], 'Callback', @(src,event)speedbuttonfun(src,event)); 

% Initialize audio feedback button object
datagui.feedbackbutton = uicontrol('Style','Togglebutton', 'String','Onsets', 'Position', [220, 20, 50, 30], 'value', 1, 'Callback', @(src,event)feedbackbuttonfun(src,event)); 

% Initialize previous button object
datagui.previousbutton = uicontrol('Style','Pushbutton', 'String','Prev.', 'Position', [300, 20, 50, 30], 'Callback', @(src,event)previousbuttonfun(src,event)); 

% Initialize next button object
datagui.nextbutton = uicontrol('Style','Pushbutton', 'String','Next', 'Position', [350, 20, 50, 30], 'Callback', @(src,event)nextbuttonfun(src,event)); 

% Initialize add button object
datagui.addbutton = uicontrol('Style','Pushbutton', 'String','Add', 'Position', [400, 20, 50, 30], 'Callback', @(src,event)addbuttonfun(src,event)); 

% Initialize delete button object
datagui.deletebutton = uicontrol('Style','Pushbutton', 'String','Delete', 'Position', [450, 20, 50, 30], 'Callback', @(src,event)deletebuttonfun(src,event)); 

% Initialize delete onsets in visible range button object
datagui.deletevisbutton = uicontrol('Style','Pushbutton', 'String','Del.vis.', 'Position', [530, 20, 50, 30], 'Callback', @(src,event)deletevisbuttonfun(src,event)); 

% Initialize time scale button object
datagui.timscalebutton = uicontrol('Style','Popupmenu', 'String',{'0.5 s.', '1 s.', '2 s.', '5 s.', '10 s.', '30 s.', '60 s.'}, 'value', 4, 'Position', [610, 20, 60, 30], 'Callback', @(src,event)timscalefun(src,event)); 

% Initialize undo button object
datagui.undobutton = uicontrol('Style','Pushbutton', 'String','Undo', 'Position', [20, 110, 50, 30], 'Callback', @(src,event)undobuttonfun(src,event)); 

% Initialize redo button object
datagui.redobutton = uicontrol('Style','Pushbutton', 'String','Redo', 'Position', [20, 80, 50, 30], 'Callback', @(src,event)redobuttonfun(src,event)); 

% Initialize keypress callback function
set(guifig, 'WindowKeyPressFcn',@(src,event)keypressfun(src,event))

% Initialize mouse scroll wheel callback function
set(guifig, 'WindowScrollWheelFcn',@(src,event)mousescrollfun(src,event))

% Initialize closerequest callback function
set(guifig, 'CloseRequestFcn', @(src,event)closefun(src, event))

guidata(guifig, datagui); % Save the GUI data


%% Screen update function

function screenupdatefun(src, ~, ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)

datagui = guidata(src); % Retrieve the GUI data
timesample = datagui.timesample; % Get the current audio time sample
timezoom = datagui.timezoom; % Get the current time zoom
onsets = datagui.onsets; % Get the onsets

% Show the TFR
curserframe = round((timesample)/datagui.fs*ons.fs_spectrum); % Find the TFR time frame for the audio time sample
if curserframe < 1, curserframe = 1; end % Constrain time frame to at least the first
if curserframe > size(inputSpectrogram_scaled,2), curserframe = size(inputSpectrogram_scaled,2); end % Constrain to maximum the last time frame
timeframes = round((-1/2*timezoom*ons.fs_spectrum:+1/2*timezoom*ons.fs_spectrum) + curserframe-1); % Time frames to show within visible range
timeframes = timeframes( timeframes>=1 & timeframes<=size(inputSpectrogram_scaled,2) ); % Constrain to at least the first and maximum the last time frame
imageplot = imagesc(inputSpectrogram_scaled(:,timeframes)); % Show the TFR image for visible time frames
set(gca,'ydir','normal')
set(gca,'clim',([-90 inputSpectrogram_scaled_max]))
ylim([1 size(inputSpectrogram_scaled,1)])
set(gca,'YTickLabel', ons.frequencies(logscale(get(gca,'YTick'))))
xlabel('Time (s.)')
ylabel('Frequency (Hz)')

% Show onsets within the visible time frame
visons = onsets( onsets >= (timeframes(1)-1)/ons.fs_spectrum & onsets <= (timeframes(end)-1)/ons.fs_spectrum ); % Visible onsets
if ~isempty(visons)
    lineplot = line( repmat( round(visons*ons.fs_spectrum) +1 - timeframes(1) , [1,2])', [0;size(inputSpectrogram_scaled,1)], 'linewidth', 1.5, 'color',[1 1 1], 'LineStyle',':'); 
end

% Show centered curser
cursersample = find( curserframe == timeframes ); % Courser time sample in relation to the visible time frames
xlim( cursersample + timezoom*1/2*ons.fs_spectrum*[-1 +1])
set(gca,'XTickLabel', (get(gca,'XTick') + timeframes(1) -1) / ons.fs_spectrum )
line( cursersample*[1,1], [0,size(inputSpectrogram_scaled,1)], 'linewidth', 1.5, 'color',[.7 .7 .7], 'LineStyle',':')

drawnow

try
    % Initialize mouse left-click callback function for TFR content
    set(imageplot, 'ButtonDownFcn',@(src,event)mouseleftclickfun(src,event))
    
    % Initialize mouse left-click callback function for onset line content
    set(lineplot, 'ButtonDownFcn',@(src,event)mouseleftclickfun(src,event))
    
catch % Fixes an issue with call from fast forwards/backwards state
end

end


%% Backwards button function

function backbuttonfun(src, ~)

datagui = guidata(src); % Retrieve the GUI data

timesample = datagui.timesample; % Current audio time sample

if strcmp(datagui.aobj.Running, 'off') && datagui.forwards == false % If audio is not playing and fast forwards is not toggled on
    
    timeoffsetprev = 0; % Reset time offset from previous clock time in seconds
    tic; % Reset timer
    
    while get(src, 'value') == 1 % While toggled on
        
        datagui = guidata(src);
        datagui.backwards = true; % Set fast forwards state to on
        timesample = datagui.timesample; % Current audio time sample
        timesamplediff = toc*fs - timeoffsetprev*fs; % Time difference based on clock time
        timeoffsetprev = toc; % Time offset from previous clock time in seconds
        
        if toc < 2 % Toggled on for less than 2 seconds
            
            timesample = round( timesample - 2*timesamplediff ); % Fast backwards at 2x speed
            
        elseif toc >= 2 && toc < 4 % Toggled on for at least 2 seconds and less than 4 seconds
            
            timesample = round( timesample - 4*timesamplediff ); % Fast backwards at 4x speed
            
        elseif toc >= 4 && toc < 6 % Toggled on for at least 4 seconds and less than 6 seconds
            
            timesample = round( timesample - 8*timesamplediff ); % Fast backwards at 8x speed
            
        elseif toc >= 6 % Toggled on for at least 6 seconds
            
            timesample = round( timesample - 16*timesamplediff ); % Fast backwards at 16x speed
            
        end
        
        if timesample < 1, timesample = 1; end % Stop at the beggining of the audio
        datagui.timesample = timesample; % Update the audio time sample
        guidata(src, datagui); % Save the data to the GUI
        
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        
    end
    
    datagui.timesample = timesample; % Update the audio time sample
    datagui.backwards = false; % Set fast forwards state to off
    guidata(src, datagui); % Save the data to the GUI
        
end

end


%% Play button function

function playbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if datagui.timesample < datagui.aobj.TotalSamples && datagui.backwards == false && datagui.forwards == false % If not reached the end of the audio or in fast backwards/forwards state
    
    % Replace play button with stop button
    stopbutton = uicontrol('Style','Pushbutton', 'String','Stop', 'Position', [70, 20, 50, 30], 'Callback', @(src,event)stopbuttonfun(src,event));
    
    play(datagui.aobj, datagui.timesample); % Playback from time sample
    
    % Update screen while playing
    while strcmp(datagui.aobj.Running, 'on')
        
        % Get the current audio time sample
        datagui = guidata(src); % Load the GUI data
        timesample = get(datagui.aobj,'CurrentSample');
        if timesample >= datagui.timesample % As long as time sample increases
            datagui.timesample = timesample;
        else % Avoid returning to the first time sample after playing to the end of the audio
            timesample = datagui.aobj.TotalSamples;
            datagui.timesample = timesample; 
        end
        guidata(src, datagui); % Save the GUI data
        
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        
    end
    
    % Replace stop button with play button
    datagui.playbutton = uicontrol('Style','Pushbutton', 'String','Play', 'Position', [70, 20, 50, 30], 'Callback', @(src,event)playbuttonfun(src,event));
    drawnow
    
end

end


%% Stop button function

function stopbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

stop(datagui.aobj); % Stop playing

end


%% Forwards button function

function forwardbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

timesample = datagui.timesample; % Get the current audio time sample

if strcmp(datagui.aobj.Running, 'off') && datagui.backwards == false % If audio is not playing and fast backwards is not toggled on
    
    timeoffsetprev = 0; % Reset time offset from previous clock time in seconds
    tic; % Reset timer
    
    while get(src, 'value') == 1 % While toggled on
        
        datagui = guidata(src);
        datagui.forwards = true; % Set fast forwards state to true
        timesample = datagui.timesample; % Current audio time sample
        timesamplediff = toc*fs - timeoffsetprev*fs; % Time difference from previous time
        timeoffsetprev = toc; % Time offset from previous clock time in seconds
        
        if toc < 2 % Toggled on for less than 2 seconds
            
            timesample = round( timesample + 2*timesamplediff ); % Fast forwards at 2x speed
            
        elseif toc >= 2 && toc < 4 % Toggled on for at least 2 seconds and less than 4 seconds
            
            timesample = round( timesample + 4*timesamplediff ); % Fast forwards at 4x speed
            
        elseif toc >= 4 && toc < 6 % Toggled on for at least 4 seconds and less than 6 seconds
            
            timesample = round( timesample + 8*timesamplediff ); % Fast forwards at 8x speed
            
        elseif toc >= 6 % Toggled on for at least 6 seconds
            
            timesample = round( timesample + 16*timesamplediff ); % Fast forwards at 16x speed
            
        end
        
        if timesample > datagui.aobj.TotalSamples, timesample = datagui.aobj.TotalSamples; end % Stop at the end of the audio
        datagui.timesample = timesample; % Update the audio time sample
        guidata(src, datagui); % Save data to the GUI
        
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        
    end
    
    datagui.timesample = timesample; % Audio time sample
    datagui.forwards = false;  % Set fast forwards state to false
    guidata(src, datagui); % Save data to the GUI
        
end

end


%% Playback speed button function

function speedbuttonfun(src, ~)

datagui = guidata(src);

if get(src,'value') == 1
    datagui.aobj.SampleRate = datagui.fs; % x1
elseif get(src,'value') == 2
    datagui.aobj.SampleRate = datagui.fs/2; % x1/2
elseif get(src,'value') == 3
    datagui.aobj.SampleRate = datagui.fs/4; % x1/4
elseif get(src,'value') == 4
    datagui.aobj.SampleRate = datagui.fs/8; % x1/8
end

guidata(src, datagui);

end


%% Feedback button function

function feedbackbuttonfun(src, ~)

datagui = guidata(src); % Load the data from the GUI

% Ensure audio is stopped while processing the onsets
if strcmp(datagui.aobj.Running, 'on')
    stop(datagui.aobj); 
    resume = true;
else
    resume = false;
end

% Toggle feedback on
if get(datagui.feedbackbutton,'value')==1
    
    % Create/update onsets feedback audio
    timesamples = (0:round(0.050*datagui.fs)-1); % 50 ms feedback time samples
    sinesamples = repmat( sin( 1000*2*pi * timesamples/fs)', [1, size(datagui.inputWave,2)]); % 1000 Hz sine wave feedback
    datagui.feedback = zeros(size(datagui.inputWave)); % Reset feedback audio
    for i=1:length(datagui.onsets) % Loop over onsets
        amplitude = max( abs( datagui.inputWave(round(datagui.onsets(i)*fs)+1 + timesamples) )); % Adjust amplitude of onset marker to sound intensity of input
        amplitude = 10^(6/20)*amplitude; % Increase onset marker by 6 dB
        datagui.feedback( round(datagui.onsets(i)*datagui.fs)+1 + timesamples , : ) = amplitude*sinesamples; % Add the onset audio
    end
    datagui.feedback = datagui.feedback(1:size(datagui.inputWave,1),:); % Trim to duration of input audio
    if max(max(abs( datagui.inputWave + datagui.feedback ))) <= 1 % Check that amplitude is within maximum tolerated
        datagui.gain = 1; % Amplitude gain parameter
    else % If audio will be clipped, lower gain to maximum tolerated
        datagui.gain = 1/max(max(abs( datagui.inputWave + datagui.feedback ))); % Amplitude gain parameter
    end
    datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.aobj.SampleRate); % Update the audioplayer object
    datagui.feedbackstatus = true; % Set feedback to true
    guidata(src, datagui); % Store the data to the GUI
    
% Toggle feedback off
elseif get(datagui.feedbackbutton,'value')==0
    
    datagui.gain = 1; % Amplitude gain parameter
    datagui.feedback = zeros(size(datagui.inputWave)); % Reset feedback audio
    datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.aobj.SampleRate); 
    datagui.feedbackstatus = false; % Set feedback to false
    guidata(src, datagui); % Store the data to the GUI
    
end

% Resume playback if it was stopped while processing the onsets
if resume == true
    playbuttonfun(src, event);
end

end


%% Previous button function

function previousbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not being played
    
    timesample = datagui.timesample; % Get the current audio time sample
    onsets = datagui.onsets; % Get the onsets
    
    % Find the last onset preceding the current time sample
    preconset = onsets(find( round(onsets*fs)+1 < timesample, 1, 'last' ));
    if ~isempty(preconset)
        datagui.timesample = round(preconset*fs)+1; % Audio time sample
        guidata(src, datagui); % Save the GUI data
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
    end
    
end

end


%% Next button function

function nextbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not being played
    
    timesample = datagui.timesample; % Get the current audio time sample
    onsets = datagui.onsets; % Get the onsets
    
    % Find the first onset following after the current time sample
    folonset = onsets(find( round(onsets*fs)+1 > timesample, 1, 'first' ));
    if ~isempty(folonset)
        datagui.timesample = round(folonset*fs)+1; % Audio time sample
        guidata(src, datagui); % Save data to the GUI
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
    end
    
end

end

%% Add button function

function addbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not being played

    % Add onset at the current audio time point if not already existing
    if sum(ismember(datagui.onsets, (datagui.timesample-1)/fs))==0

        disp(['Adding onset at ',num2str((datagui.timesample-1)/fs),' s.'])
        datagui.onsets = unique(vertcat(datagui.onsets, (datagui.timesample-1)/fs));
        datagui.undo = vertcat(datagui.undo, [(datagui.timesample-1)/fs, true]); % Add actions for undoing
        datagui.redo = []; % Reset actions for redoing
        guidata(src, datagui); % Save the GUI data

        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        feedbackbuttonfun(src, []) % Update the audio feedback
        autosavefun(src, [])

    end
end

end


%% Delete button function

function deletebuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not being played

    % Delete any existing onset at the current audio time point
    if sum(ismember( round(datagui.onsets*fs)+1 , datagui.timesample))>0

        disp(['Deleting onset at ',num2str((datagui.timesample-1)/fs),' s.'])
        datagui.onsets = datagui.onsets(~ismember(round(datagui.onsets*fs)+1, datagui.timesample));
        datagui.undo = vertcat(datagui.undo, [(datagui.timesample-1)/fs, false]); % Add actions for undoing
        datagui.redo = []; % Reset actions for redoing
        guidata(src, datagui); % Save the GUI data

        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        feedbackbuttonfun(src, []) % Update the audio feedback
        autosavefun(src, [])

    end
end

end

%% Delete all onsets in the visible range button function

function deletevisbuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not being played

    if strcmp( questdlg('Delete all onsets in visible range?', '', 'Delete', 'Cancel', 'Cancel') , 'Delete') % Request user to verify batch deletion before proceeding

        % Find any existing onsets in the visible time range
        timesec = (datagui.timesample-1)/fs; % Selected time point in seconds
        onsetvis = datagui.onsets( datagui.onsets > timesec - datagui.timezoom/2 & datagui.onsets < timesec + datagui.timezoom/2 ); % Find onsets in visible time range

        % Delete all onsets in the visible time range
        for i=1:length(onsetvis)

            disp(['Deleting onset at ',num2str(onsetvis(i)),' s.'])
            datagui.onsets = datagui.onsets(~ismember(datagui.onsets, onsetvis(i)));
            datagui.undo = vertcat(datagui.undo, [onsetvis(i), false]); % Add actions for undoing
            datagui.redo = []; % Reset actions for redoing
            guidata(src, datagui); % Save the GUI data

        end

        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        feedbackbuttonfun(src, []) % Update the audio feedback
        autosavefun(src, [])

    end
end

end


%% Time scale button function

function  timscalefun(src, ~)

datagui = guidata(src);
if get(src,'value') == 1
    datagui.timezoom = 0.5; % 0.5 s.
elseif get(src,'value') == 2
    datagui.timezoom = 1; % 1 s.
elseif get(src,'value') == 3
    datagui.timezoom = 2; % 2 s.
elseif get(src,'value') == 4
    datagui.timezoom = 5; % 5 s.
elseif get(src,'value') == 5
    datagui.timezoom = 10; % 10 s.
elseif get(src,'value') == 6
    datagui.timezoom = 30; % 30 s.
elseif get(src,'value') == 7
    datagui.timezoom = 60; % 60 s.
end
guidata(src, datagui);

screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)

end


%% Undo button function

function undobuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not playing
    
    if ~isempty(datagui.undo) % If there are actions to undo
        
        if datagui.undo(end,2)==false % If the last edit was to delete an onset
            
            % Add the onset again
            disp(['Undoing deletion of onset at ',num2str(datagui.undo(end,1)),' s.'])
            datagui.onsets = unique(vertcat(datagui.onsets, datagui.undo(end,1)));
            datagui.redo = vertcat(datagui.redo, datagui.undo(end,:)); % Add last edit for redoing
            datagui.undo(end,:) = []; % and remove it from undoing
            
        elseif datagui.undo(end,2)==true % If the last edit was to add an onset
            
            % Delete the added onset
            disp(['Undoing adding of onset at ',num2str(datagui.undo(end,1)),' s.'])
            datagui.onsets = datagui.onsets(~ismember(datagui.onsets, datagui.undo(end,1)));
            datagui.redo = vertcat(datagui.redo, datagui.undo(end,:)); % Add last edit for redoing
            datagui.undo(end,:) = []; % and remove it from undoing
            
        end
        
        datagui.timesample = round(datagui.redo(end,1)*fs)+1; % Move to the affected time sample
        if datagui.timesample < 1, datagui.timesample = 1; end % Ensure the time sample is at least the first
        if datagui.timesample > datagui.aobj.TotalSamples, datagui.timesample = datagui.aobj.TotalSamples; end % Ensure the time sample is maximum the last
        
        guidata(src, datagui); % Save the GUI data
        
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        feedbackbuttonfun(src, []) % Update the audio feedback
        autosavefun(src, [])
        
    else
        
        disp('(All previously performed edits are undone.)')
        
    end
    
end

end


%% Redo button function

function redobuttonfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not playing
    
    if ~isempty(datagui.redo) % If there are any undone actions to redo
        
        if datagui.redo(end,2)==false % If the last undone edit was to add a deleted onset
            
            % Delete the onset again
            disp(['Redoing deletion of onset at ',num2str(datagui.redo(end,1)),' s.'])
            datagui.onsets = datagui.onsets(~ismember(datagui.onsets, datagui.redo(end,1)));
            datagui.undo = vertcat(datagui.undo, datagui.redo(end,:)); % Add last edit for undoing
            datagui.redo(end,:) = []; % and remove it from redoing
            
        elseif datagui.redo(end,2)==true % If the last undone edit was to delete an added onset
            
            % Add the onset again
            disp(['Redoing adding of onset at ',num2str(datagui.redo(end,1)),' s.'])
            datagui.onsets = unique(vertcat(datagui.onsets, datagui.redo(end,1)));
            datagui.undo = vertcat(datagui.undo, datagui.redo(end,:)); % Add last edit for undoing
            datagui.redo(end,:) = []; % and remove it from redoing
            
        end
        
        datagui.timesample = round(datagui.undo(end,1)*fs)+1; % Move to the affected time sample
        if datagui.timesample < 1, datagui.timesample = 1; end % Ensure the time sample is at least the first
        if datagui.timesample > datagui.aobj.TotalSamples, datagui.timesample = datagui.aobj.TotalSamples; end % Ensure the time sample is maximum the last
        
        guidata(src, datagui); % Save the GUI data
        
        screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
        feedbackbuttonfun(src, []) % Update the audio feedback
        autosavefun(src, [])
        
    else
        
        disp('(All previously undone edits are redone.)')
        
    end
    
end

end


%% Key press function

function keypressfun(src, event)
    
datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not playing
    
    % Start playing when pressing play ket (default is space)
    if strcmp(event.Key, datagui.name_onsed_keys.play)
        
        playbuttonfun(src, event);
        
        
    % Fast backwards when pressing left (default key) and not already running
    elseif strcmp(event.Key, datagui.name_onsed_keys.backwards) && datagui.backwards == false && datagui.forwards == false
        
        datagui.backwards = true; % Set fast backwards state to true
        set(datagui.backbutton, 'value', 1) % Set backwards button to toggled on
        guidata(src, datagui); % Save the GUI data
        
        timeoffsetprev = 0; % Reset time offset from previous time
        tic; % Reset timer
        
        while datagui.backwards == true && datagui.timesample > 1 % While toggled on and the beginning of the audio has not been reached
            
            datagui = guidata(src); % Load the GUI data
            timesample = datagui.timesample; % Audio time sample
            timesamplediff = toc*fs - timeoffsetprev*fs; % Time difference from previous
            timeoffsetprev = toc; % Time offset from previous clock time
            
            if toc < 2 % Toggled on for less than 2 seconds
                
                timesample = round( timesample - 2*timesamplediff ); % Fast backwards at 2x speed
                
            elseif toc >= 2 && toc < 4 % Toggled on for at least 2 seconds and less than 4 seconds
                
                timesample = round( timesample - 4*timesamplediff ); % Fast backwards at 4x speed
                
            elseif toc >= 4 && toc < 6 % Toggled on for at least 4 seconds and less than 6 seconds
                
                timesample = round( timesample - 8*timesamplediff ); % Fast backwards at 8x speed
                
            elseif toc >= 6 % Toggled on for at least 6 seconds
                
                timesample = round( timesample - 16*timesamplediff ); % Fast backwards at 16x speed
                
            end
            
            if timesample < 1, timesample = 1; end % Ensure the time sample is at least the first
            datagui.timesample = timesample;
            guidata(src, datagui); % Save the GUI data
            
            screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
            
        end
        
        datagui.timesample = timesample;
        datagui.backwards = false; % Set fast backwards to false
        set(datagui.backbutton, 'value', 0) % Set backwards button to toggled off
        guidata(src, datagui); % Save the GUI data
        
        
    % Fast forwards when pressing right and not already running
    elseif strcmp(event.Key, datagui.name_onsed_keys.forwards) && datagui.backwards == false && datagui.forwards == false
        
        datagui.forwards = true; % Set fast forwards to true
        set(datagui.forwardbutton, 'value', 1) % Set forwards button to toggled on
        guidata(src, datagui); % Save the GUI data
        
        timeoffsetprev = 0; % Reset time offset from previous clock time
        tic; % Reset timer
        
        while datagui.forwards == true && datagui.timesample < datagui.aobj.TotalSamples % While toggled on and not reached the end of the audio
            
            datagui = guidata(src); % Load the GUI data
            timesample = datagui.timesample;
            timesamplediff = toc*fs - timeoffsetprev*fs; % Difference in time
            timeoffsetprev = toc; % Time offset from previous clock time
            
            if toc < 2 % Toggled on for less than 2 seconds
                
                timesample = round( timesample + 2*timesamplediff ); % Fast forwards at 2x speed
                
            elseif toc >= 2 && toc < 4 % Toggled on for at least 2 seconds and less than 4 seconds
                
                timesample = round( timesample + 4*timesamplediff ); % Fast forwards at 4x speed
                
            elseif toc >= 4 && toc < 6 % Toggled on for at least 4 seconds and less than 6 seconds
                
                timesample = round( timesample + 8*timesamplediff ); % Fast forwards at 8x speed
                
            elseif toc >= 6 % Toggled on for at least 6 seconds
                
                timesample = round( timesample + 16*timesamplediff ); % Fast forwards at 16x speed
                
            end
            
            if timesample > datagui.aobj.TotalSamples, timesample = datagui.aobj.TotalSamples; end % Ensure the time sample is maximum the last
            datagui.timesample = timesample;
            guidata(src, datagui); % Save the GUI data
            
            screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
            
        end
        
        datagui.timesample = timesample;
        datagui.forwards = false; % Set fast forwards state to false
        set(datagui.forwardbutton, 'value', 0) % Set forwards button to toggled off
        guidata(src, datagui); % Save the GUI data
        
        
    % Previous onset when pressing home (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.previous)
        
        previousbuttonfun(src, [])
        
        
    % Next onset when pressing end (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.next)
        
        nextbuttonfun(src, [])
        
        
    % Add onset when pressing insert (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.add)
        
        addbuttonfun(src, [])
        autosavefun(src, [])
        
        
    % Delete onset when pressing delete (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.delete)
        
        deletebuttonfun(src, [])
        autosavefun(src, [])
        
        
    % Delete all onsets in visible range when pressing d (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.deletevis)
        
        deletevisbuttonfun(src, [])
        autosavefun(src, [])
        
        
    % Undo last action when pressing z (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.undo)
        
        undobuttonfun(src, [])
        autosavefun(src, [])
        
        
    % Redo last action when pressing x (default)
    elseif strcmp(event.Key, datagui.name_onsed_keys.redo)
        
        redobuttonfun(src, [])
        autosavefun(src, [])
        
    end
    
    
    % Stop any running fast backwards/forwards when pressing any key
    if datagui.backwards == true || datagui.forwards == true
        
        datagui.backwards = false;
        datagui.forwards = false;
        guidata(src, datagui);
        
    end
    
    
elseif strcmp(datagui.aobj.Running, 'on')
    
    % Stop playing when pressing space (default) and playing
    if strcmp(event.Key, datagui.name_onsed_keys.play)
        
        stop(datagui.aobj);
        
    end
end


% Decrease/increase playback speed when pressing pagedown/pageup (default)
if strcmp(event.Key, datagui.name_onsed_keys.speed_decrease)
    
    % Increase speed in speedbutton object
    if get(datagui.speedbutton,'value') + 1 <= 4
        set(datagui.speedbutton,'value', get(datagui.speedbutton,'value') + 1);
    end
    
elseif strcmp(event.Key, datagui.name_onsed_keys.speed_increase)
    
    % Decrease speed in speedbutton object
    if get(datagui.speedbutton,'value') - 1 >= 1
        set(datagui.speedbutton,'value', get(datagui.speedbutton,'value') - 1);
    end
    
end
if strcmp(event.Key, datagui.name_onsed_keys.speed_increase) || strcmp(event.Key, datagui.name_onsed_keys.speed_decrease)
    
    % Update the playback sampling rate in the audioplayer object
    if get(datagui.speedbutton,'value') == 1
        datagui.aobj.SampleRate = datagui.fs; % x1
    elseif get(datagui.speedbutton,'value') == 2 
        datagui.aobj.SampleRate = datagui.fs/2; % x1/2
    elseif get(datagui.speedbutton,'value') == 3
        datagui.aobj.SampleRate = datagui.fs/4; % x1/4
    elseif get(datagui.speedbutton,'value') == 4
        datagui.aobj.SampleRate = datagui.fs/8; % x1/8
    end
    guidata(src, datagui);
    
end

% Toggle onsets audio feedback on/off when pressing o (default)
if strcmp(event.Key, datagui.name_onsed_keys.onsetsfeedback)
    
    % Ensure audio is stopped while processing the onsets
    if strcmp(datagui.aobj.Running, 'on')
        stop(datagui.aobj);
        resume = true;
    else
        resume = false;
    end
    
    if datagui.feedbackstatus == false % If audio feedback is off, toggle it on
        
        set(datagui.feedbackbutton,'value', 1) % Set feedback button to toggled on
        
        % Create/update onsets feedback audio
        timesamples = (0:round(0.050*datagui.fs)-1); % 50 ms feedback time samples
        sinesamples = repmat( sin( 1000*2*pi * timesamples/fs)', [1, size(datagui.inputWave,2)]); % 1000 Hz sine wave feedback
        datagui.feedback = zeros(size(datagui.inputWave)); % Reset feedback audio
        for i=1:length(datagui.onsets) % Loop over onsets
            amplitude = max( abs( datagui.inputWave(round(datagui.onsets(i)*fs)+1 + timesamples) )); % Adjust amplitude of onset marker to sound intensity of input
            amplitude = 10^(6/20)*amplitude; % Increase onset marker by 6 dB
            datagui.feedback( round(datagui.onsets(i)*datagui.fs)+1 + timesamples , : ) = amplitude*sinesamples; % Add the onset audio
        end
        datagui.feedback = datagui.feedback(1:size(datagui.inputWave,1),:); % Trim to duration of input audio
        if max(max(abs( datagui.inputWave + datagui.feedback ))) <= 1 % Check that amplitude is within maximum tolerated
            datagui.gain = 1; % Amplitude gain parameter
        else % If audio will be clipped, lower gain to maximum tolerated
            datagui.gain = 1/max(max(abs( datagui.inputWave + datagui.feedback ))); % Amplitude gain parameter
        end
        datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.aobj.SampleRate);
        datagui.feedbackstatus = true; % Set feedback to true
        guidata(src, datagui); % Save the GUI data
        
    elseif datagui.feedbackstatus == true % If audio feedback is on, toggle it off
        
        set(datagui.feedbackbutton,'value', 0) % Set feedback button to toggled off
        
        datagui.gain = 1; % Amplitude gain parameter
        datagui.feedback = zeros(size(datagui.inputWave)); % Reset feedback audio
        datagui.aobj = audioplayer( datagui.gain*( datagui.inputWave + datagui.feedback ) , datagui.aobj.SampleRate);
        datagui.feedbackstatus = false; % Set feedback to false
        guidata(src, datagui); % Save the GUI data
        
    end
    
    % Resume playback if it was stopped while processing the onsets
    if resume == true
        playbuttonfun(src, event);
    end
    
end


% Decrease/increase time scale when pressing -/+ (default)
if strcmp(event.Key, datagui.name_onsed_keys.zoomout)
    
    % Increase time scale in the time scale button object
    if get(datagui.timscalebutton,'value') + 1 <= 7
        set(datagui.timscalebutton,'value', get(datagui.timscalebutton,'value') + 1);
    end
    
elseif strcmp(event.Key, datagui.name_onsed_keys.zoomin)
    
    % Decrease time scale in the time scale button object
    if get(datagui.timscalebutton,'value') - 1 >= 1
        set(datagui.timscalebutton,'value', get(datagui.timscalebutton,'value') - 1);
    end
    
end
if strcmp(event.Key, datagui.name_onsed_keys.zoomout) || strcmp(event.Key, datagui.name_onsed_keys.zoomin)
    
    % Increase/decrease time scale in the visible TFR range
    if get(datagui.timscalebutton,'value') == 1
        datagui.timezoom = 0.5; % 0.5 s.
    elseif get(datagui.timscalebutton,'value') == 2
        datagui.timezoom = 1; % 1 s.
    elseif get(datagui.timscalebutton,'value') == 3
        datagui.timezoom = 2; % 2 s.
    elseif get(datagui.timscalebutton,'value') == 4
        datagui.timezoom = 5; % 5 s.
    elseif get(datagui.timscalebutton,'value') == 5
        datagui.timezoom = 10; % 10 s.
    elseif get(datagui.timscalebutton,'value') == 6
        datagui.timezoom = 30; % 30 s.
    elseif get(datagui.timscalebutton,'value') == 7
        datagui.timezoom = 60; % 60 s.
    end
    guidata(src, datagui); % Save the GUI data
    
    screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
    
end

end


%% Mouse left-click function

function mouseleftclickfun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'off') % If audio is not playing
    
    % Center the time frame on left mouse click
    timesample = datagui.timesample; % Audio time sample
    timezoom = datagui.timezoom; % Time zoom in s.
    curserframe = round((timesample-1)/fs*ons.fs_spectrum)+1; % Current time frame
    if curserframe < 1, curserframe = 1; end % Ensure current time frame is at least the first
    if curserframe > size(inputSpectrogram_scaled,2), curserframe = size(inputSpectrogram_scaled,2); end % Ensure current time frame is at maximum the last
    timeframes = (-1/2*timezoom*ons.fs_spectrum:+1/2*timezoom*ons.fs_spectrum) + curserframe-1; % Visible time frames
    timeframes = timeframes(timeframes>0); % Ensure the visible time frames are at least the first
    axesHandle  = get(src,'Parent'); % Get the mouse left click coordinates...
    coordinates = get(axesHandle,'CurrentPoint');
    coordinates = ceil(coordinates(1, 1));
    curserframe = timeframes(coordinates); % Set the current time frame to the left mouse click x-coordinate
    if curserframe < 1, curserframe = 1; end % Ensure current time frame is at least the first
    if curserframe > size(inputSpectrogram_scaled,2), curserframe = size(inputSpectrogram_scaled,2); end % Ensure current time frame is at maximum the last
    timesample = floor(fs/ons.tres*(curserframe-1))+1; % Convert the current time frame in the TFR to the corresponding audio time sample
    if timesample < 1, timesample = 1; end % Ensure the audio time sample is at least the first
    if timesample > datagui.aobj.TotalSamples, timesample = datagui.aobj.TotalSamples; end % Ensure the audio time sample is maximum the last
    datagui.timesample = timesample; % Store the audio time sample
    guidata(src, datagui); % Save the GUI data
    
    screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)
    
end

end


%% Mouse scroll wheel function

function mousescrollfun(src, event)

datagui = guidata(src); % Load the GUI data

if event.VerticalScrollCount > 0 % If scrolling backwards
    
    % Zoom out in the time scale button object
    if get(datagui.timscalebutton,'value') + 1 <= 7
        set(datagui.timscalebutton,'value', get(datagui.timscalebutton,'value') + 1);
    end

elseif event.VerticalScrollCount < 0 % If scrolling forwards
    
    % Zoom in in the time scale button object
    if get(datagui.timscalebutton,'value') - 1 >= 1
        set(datagui.timscalebutton,'value', get(datagui.timscalebutton,'value') - 1);
    end

end
if event.VerticalScrollCount > 0 || event.VerticalScrollCount < 0
    
    % Zoom in/out in the TFR time frames
    if get(datagui.timscalebutton,'value') == 1
        datagui.timezoom = 0.5; % 0.5 s.
    elseif get(datagui.timscalebutton,'value') == 2
        datagui.timezoom = 1; % 1 s.
    elseif get(datagui.timscalebutton,'value') == 3
        datagui.timezoom = 2; % 2 s.
    elseif get(datagui.timscalebutton,'value') == 4
        datagui.timezoom = 5; % 5 s.
    elseif get(datagui.timscalebutton,'value') == 5
        datagui.timezoom = 10; % 10 s.
    elseif get(datagui.timscalebutton,'value') == 6
        datagui.timezoom = 30; % 30 s.
    elseif get(datagui.timscalebutton,'value') == 7
        datagui.timezoom = 60; % 60 s.
    end
    guidata(src, datagui); % Save the GUI data

    screenupdatefun(src, [], ons, inputSpectrogram_scaled, inputSpectrogram_scaled_max, logscale)

end

end


%% Auto-save function

function autosavefun(src, ~)

datagui = guidata(src); % Load the GUI data

onsetsedit = struct; % Reset the onsets editing data
onsetsedit.audiopath = datagui.audiopath; % Store the path to the audio file
onsetsedit.audiofile = datagui.audiofile; % Store name of the audio file
onsetsedit.outputpath = datagui.outputpath; % Store the path to the onsets edit output file
onsetsedit.outputfile = datagui.outputfile; % Store the name of the onsets edit output file
onsetsedit.timesample = datagui.timesample; % Store the current audio time sample
onsetsedit.undo = datagui.undo; % Store actions to undo
onsetsedit.redo = datagui.redo; % Store actions to redo
onsets = datagui.onsets; % Store the latest edited onsets in the 'onsets' vector variable
save(datagui.output, 'onsets', 'onsetsedit') % Save the edited onsets to the output file

end


%% Close function

function closefun(src, ~)

datagui = guidata(src); % Load the GUI data

if strcmp(datagui.aobj.Running, 'on') % If audio is playing
    
    stop(datagui.aobj); % Ensure audio is stopped when closing the editor
    
else
    
    datagui.close = true; % Set close state to true
    guidata(src, datagui); % Save the GUI data
    
    disp('Goodbye! See you later.')
    onsetsedit = struct;
    onsetsedit.audiopath = datagui.audiopath; % Store the path to the audio file
    onsetsedit.audiofile = datagui.audiofile; % Store the name of the audio file
    onsetsedit.outputpath = datagui.outputpath; % Store the path to the onsets edit output file
    onsetsedit.outputfile = datagui.outputfile; % Store the name of the onsets edit output file
    onsetsedit.timesample = datagui.timesample; % Store the audio time sample
    onsetsedit.undo = datagui.undo; % Store actions to undo
    onsetsedit.redo = datagui.redo; % Store actions to redo
    onsets = datagui.onsets; % Store the latest edited onsets in the 'onsets' vector variable
    save(datagui.output, 'onsets', 'onsetsedit') % Save the edited onsets to the output file
    
    delete(gcf) % Close the GUI
    
end
end


end
