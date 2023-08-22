# Naturalistic auditory MEG/EEG (NAME) package

Under development at Center for Music in the Brain ( https://musicinthebrain.au.dk/ )

Please cite this site ( https://github.com/nielsthaumann/nameeg ) if you are publishing results using functions from this package. 

----------------------------------------
## Automatic sound onset detection with noise suppression

The algorithm includes suppression of noise in the audio medium and vibrato/tremolo </br>
and the ability to detect slow attacks. 

Use as: 

onsets = name_ons(audio)  

to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'

The onsets are given as an array of time points in seconds and saved to an output file </br>
(by default with the same path and name as the input file added the ending ..._onsets.mat).

----------------------------------------
## Sound onset editor

Use this editor to combine the speed of automatic sound onset detection 
with the reliability of manual editing.

Call the editor from the Matlab Command Window and hit Enter: </br>
name_onsed

Load an audio file (.wav, .flac, .mp3, .m4a, .mp4, .ogg) for auditory validation
and visual inspection of energy increases in a time-frequency representation. 

Load and save a .mat file with the variable 'onsets' containing a numerical vector 
with the onsets stored in seconds. 

Changes are saved automatically (along with the undo/redo history). 
