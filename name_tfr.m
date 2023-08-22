function [tfr_a, tfr_p, win, srtfr, freq, time] = name_tfr(timeseries, srtime, varargin)

% Calculate the time-frequency representation (TFR) for a time series
% 
% Use as: 
% 
% [tfr_a, tfr_p, win, srtfr, freq, time] = name_tfr(timeseries, srtime)
% 
% 
% Required inputs are: 
% 
% timeseries    numeric time by channel array or matrix with time samples along the first
%               dimension and channels along the second dimension
%
% srtime        sampling rate of the time series in Hz
% 
%  
% Optional input arguments are: 
% 
% [...] = name_tfr(..., 'srtfr',srtfr)           srtfr is the sampling rate of the TFR in Hz (default is 100)
% 
% [...] = name_tfr(..., 'framesize',framesize)   frame size for FFT in number of time samples
%                                                (default is 100 ms rounded to an even number of time samples)
% 
% [...] = name_tfr(..., 'frametime',frametime)   frame time reference can be 'center' or 'end' (default is 'center')
% 
% 
% Outputs are: 
% 
% tfr_a         amplitude spectra over time in a numeric frequency by time by channel matrix 
%               with the frequnecy bins along the first dimension, time frames
%               along the second dimension, and channels along the third dimension
%    
% tfr_p         phase spectra over time in radians (between -pi and +pi) in a numeric frequency 
%               by time by channel matrix with the frequency bins along the first dimension, 
%               time frames along the second dimension, and channels along the third dimension
% 
% win           Blackman-Harris window
% 
% srtfr         actual sampling rate of the TFR in Hz
% 
% freq          frequencies in Hz
%
% time          time in seconds
%
% 
% Beta version 20230607. 
% 
% name_tfr is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 


% Parse the input arguments
p = inputParser; 
addOptional(p, 'srtfr', 100) % (default sampling rate of the TFR is 100 Hz)
addOptional(p, 'framesize', 2*round(0.100*srtime/2) ) % (default frame size for FFT is 100 ms rounded to an even number of time samples)
addOptional(p, 'frametime', 'center') % (default frame time reference is 'center')
parse(p, varargin{:})
srtfr = p.Results.srtfr; % Desired sampling rate of output
framesize = p.Results.framesize; % Frame size for FFT
frametime = p.Results.frametime; % Frame time reference

% Prepare the input and output variables
shiftsize = round(srtime/srtfr); % Time frame shift size
srtfr = 1/(shiftsize/srtime); % Actual sampling rate of output
if license('test','Signal_Toolbox')
    win = blackmanharris(framesize); % Blackman-Harris window (NB: Requires the Signal Processing Toolbox)
else
    n = (0:framesize-1);
    win = ( 0.35875 - 0.48829*cos(2*pi*n/(framesize-1)) + 0.14128*cos(4*pi*n/(framesize-1)) - 0.01168*cos(6*pi*n/(framesize-1)) )'; % Blackman-Harris window
end
freq = (srtime*(0:floor(framesize/2))/framesize)'; % Frequencies in Hz

% Prepare matrix indexing for parallel Fourier transform across time frames
frameaddsamples = (0:framesize-1)'; % Time samples in the FFT time frame
frameoffsets = (1:shiftsize:size(timeseries,1))'; % Offset in time samples for each FFT time frame
index = repmat(frameaddsamples, [1, length(frameoffsets)]) + repmat(frameoffsets', [length(frameaddsamples), 1]); % Matrix indexing ( FFT frame time samples , time frame number )

if strcmp(frametime, 'center')
    
    % Prepare zero-padded input time series centered at the center of the time frame
    disp('Input time series is centered on the center time sample in the analysis time frame.')
    timeserieszp = vertcat(zeros(framesize/2-1, size(timeseries,2)), timeseries, zeros(framesize/2-1, size(timeseries,2))); % Zero-padded input wave
    
elseif strcmp(frametime, 'end')
    
    % Prepare zero-padded input time series centered at the end of the time frame
    disp('Input time series is centered on the last time sample in the analysis time frame.')
    timeserieszp = vertcat(zeros(framesize-1, size(timeseries,2)), timeseries); % Zero-padded input wave
    
end

% Calculate the TFR
tfr_a = zeros(floor(framesize/2)+1, length(frameoffsets), size(timeseries,2));        % Prepare TFR amplitude (frequency, time, channel) matrix
tfr_p = zeros(floor(framesize/2)+1, length(frameoffsets), size(timeseries,2)); % Prepare TFR phase in radians (frequency, time, channel) matrix
for c=1:size(timeseries,2) % Loop over channels
    
    disp(['Performing fast Fourier transform of ',num2str(length(freq)),' frequency bands in ',num2str(length(frameoffsets)),' time frames in channel ',num2str(c), ' of ',num2str(size(timeseries,2)),'.'])
    
    % Parallel FFT calcuation across time windows
    specttemp = fft( repmat(win, [1, length(frameoffsets)]) .* reshape( timeserieszp(index, c) , [length(frameaddsamples), length(frameoffsets)]) );
    
    % Save the spectra angles in radians for later possibility of inverse Fourier transform
    tfr_p(:,:,c) = angle(specttemp(1:floor(framesize/2)+1, :));
    
    % Obtain the single-sided amplitude spectra
    specttemp = abs(specttemp/framesize);
    tfr_a(:,:,c) = specttemp(1:floor(framesize/2)+1, :);
    tfr_a(2:end-1,:,c) = 2*tfr_a(2:end-1,:,c);
    
end

time = 0:1/srtfr:(size(tfr_a,2)-1)/srtfr; % Store the TFR time in seconds

end