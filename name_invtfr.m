function timeseries = name_invtfr(tfr_a, tfr_p, win, srtfr, srtime, frametime)

% Calculate the inverse time-frequency representation (TFR) for a time series
% 
% Use as: 
% 
% timeseries = name_invtfr(tfr_a, tfr_p, win, srtfr, srtime, frametime)
% 
% 
% Required inputs are: 
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
% srtfr         sampling rate of the TFR in Hz
% 
% srtime        sampling rate of the time series in Hz
% 
% frametime     frame time reference, either 'center' or 'end'
% 
%  
% Outputs is: 
% 
% timeseries    numeric time by channel array or matrix with time samples along the first
%               dimension and channels along the second dimension
%
% 
% Beta version 20230607. 
% 
% name_invtfr is part of the Naturalistic Auditory MEG/EEG (NAME) package. https://github.com/nielsthaumann/nameeg
% 


% Calculate the inverse TFR
framesize = 2*(size(tfr_a,1)-1); % FFT frame size
shiftsize = round(srtime/srtfr); % FFT frame shift size
invwindow = win/sum(win(1:shiftsize:end).^2); % Inversion window function
tfr_a(2:end-1,:) = 1/2*tfr_a(2:end-1,:); % Revert the amplitude weighting...
tfr_a = tfr_a*framesize; 
timeseries = zeros(size(tfr_a,2)*srtime/srtfr+framesize, size(tfr_a,3)); % Prepare the output time series
for c=1:size(tfr_a,3) % Loop over channels
    
    disp(['Performing inverse fast Fourier transform of ',num2str(size(tfr_a,1)),' frequency bands in ',num2str(size(tfr_a,2)),' time frames in channel ',num2str(c), ' of ',num2str(size(tfr_a,3)),'.'])
    
    % Calculate the inverse FFT weighted by the inversion window function
    timetemp = 2 * repmat(invwindow , [1, size(tfr_a,2)] ).* real(ifft( tfr_a(:,:,c) .* exp(1i*tfr_p(:,:,c)) , framesize));
    
    % Sum the inverse FFT frames over time to reconstruct the time series
    for j=1:size(tfr_a,2) % Loop over time frames
        timesamples = (j-1)*shiftsize +1 : (j-1)*shiftsize +framesize; % Time samples
        timeseries(timesamples, c) = timeseries(timesamples, c) + timetemp(:,j); % Add the inverse FFT frame to the time series
    end
end
clear('timetemp'); % Cleanup memory

if strcmp(frametime, 'center')
    
    % Remove zero-padding from input time series centered at the center of the time frame
    disp('Removing zero-padding from time series centered on the center time sample in the analysis time frame.')
    timeseries = timeseries(framesize/2:(size(tfr_a,2))*shiftsize+framesize/2, :,:);
    
elseif strcmp(frametime, 'end')
    
    % Remove zero-padding from input time series centered at the end of the time frame
    disp('Removing zero-padding from time series centered on the last time sample in the analysis time frame.')
    timeseries = timeseries(framesize:(size(tfr_a,2))*shiftsize+framesize, :,:);
    
end

end