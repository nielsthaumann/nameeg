# Naturalistic auditory MEG/EEG (NAME) package

Under development at Center for Music in the Brain ( https://musicinthebrain.au.dk/ )

Please cite this site ( https://github.com/nielsthaumann/nameeg ) if you are publishing results using functions from this package. 

----------------------------------------
## Automatic sound onset detection with noise suppression

Use as: 

onsets = name_ons(audio)  

to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'

The onsets are given as an array of time points in seconds and saved to an output file 
(by default with the same path and name as the input file added the ending ..._onsets.mat)
