# Naturalistic auditory MEG/EEG (NAME) package

Under development at Center for Music in the Brain ( https://musicinthebrain.au.dk/ )

Please cite this site ( https://github.com/nielsthaumann/nameeg ) if you are publishing results using functions from this package. 

----------------------------------------
## Automatic sound onset detection with noise suppression

Automatic detection of sound onsets in naturalistic music
presents a methodologically complicated challenge, as the
audio signal typically contains multiple overlapping voices
with complex spectral-temporal patterns and background
noise (e.g., [Smith & Fraser, 2004](https://doi.org/10.1109/TNN.2004.832831); [Thoshkahna &
Ramakrishnan, 2008](http://dx.doi.org/10.1109/ICOSP.2008.4697399); [Alías et al., 2016](https://doi.org/10.3390/app6050143)). 
An additional challenge in measuring evoked responses 
to sound onsets in recorded music is the brevity of the 
auditory cortical evoked responses, which last only a few tens of 
milliseconds reversing polarity, thereby requiring high temporal accuracy
([Haumann et al., 2021](https://doi.org/10.1016/j.brainres.2020.147248)). 

(https://github.com/nielsthaumann/nameeg/blob/main/01%20mironsets%20accuracy%20improvement.png)

Brownian noise affecting low frequency bands in the audio
can be suppressed with the non-negative least squares
function (lsqnonneg) to fit a 1/f curve to the average
amplitude spectrum. 

(https://github.com/nielsthaumann/nameeg/blob/main/02%20simulated%20noise.png)

This fitted curve serves as a noise floor
threshold, below which all audio can be removed. 

(https://github.com/nielsthaumann/nameeg/blob/main/02%20simulated%20noise.png) (https://github.com/nielsthaumann/nameeg/blob/main/03%20noise%20suppression.png)


The algorithm includes suppression of noise in the audio medium. 

Use as: 

onsets = name_ons(audio)  

to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'

The onsets are given as an array of time points in seconds.

Please make sure that a recent version of the MIRtoolbox is installed. 
(This function was tested with MIRtoolbox 1.8.1.)
E.g, visit: https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mirtoolbox

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
