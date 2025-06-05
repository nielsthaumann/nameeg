# Naturalistic auditory MEG/EEG (NAME) package

Under development at Center for Music in the Brain ( https://musicinthebrain.au.dk/ )

Please cite this site ( https://github.com/nielsthaumann/nameeg ) if you are publishing results using functions from this package. 


----------------------------------------
## Automatic sound onset detection with noise suppression

Use the onset detection function with the following command: 

onsets = name_ons(audio)  

to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'

The onsets are given as an array of time points in seconds.
<br></br>
Automatic detection of sound onsets in naturalistic music
presents a methodologically complicated challenge, 
as the audio signal typically contains multiple overlapping voices 
with complex spectral-temporal patterns and background noise 
(e.g., [Smith & Fraser, 2004](https://doi.org/10.1109/TNN.2004.832831); [Thoshkahna &
Ramakrishnan, 2008](http://dx.doi.org/10.1109/ICOSP.2008.4697399); [Alías et al., 2016](https://doi.org/10.3390/app6050143)). <br>

An additional challenge in measuring evoked responses to sound onsets
in recorded music is the brevity of the auditory cortical evoked responses,
which last only a few tens of milliseconds reversing polarity, 
thereby requiring high temporal accuracy
([Haumann et al., 2021](https://doi.org/10.1016/j.brainres.2020.147248)). 

<p align="center">
  <image width="567" height="567" src="https://ars.els-cdn.com/content/image/1-s2.0-S0006899320306065-ga1_lrg.jpg">
</p>
<p align="center">
  <a href=https://doi.org/10.1016/j.brainres.2020.147248> (Haumann et al., 2021, "Extracting human cortical responses to sound onsets <br> and acoustic feature changes in real music, and their relation to event rate",<i> Brain Research</i>)</br>
  </a>
</br></br>

The automatic sound onset detection algorithm includes suppression of noise in the audio medium. 

<p align="center">
  <image width="700" height="526" src="https://github.com/nielsthaumann/nameeg/blob/main/name_ons_01.png">
</p>

Brownian noise affecting low frequency bands in the audio
is suppressed with the non-negative least squares
function (lsqnonneg) to fit a 1/f curve to the average
amplitude spectrum. 

<p align="center">
  <image width="700" height="526" src="https://github.com/nielsthaumann/nameeg/blob/main/name_ons_02.png">
</p>

The fitted curve serves as a noise floor
threshold, below which all audio is removed. 

<p align="center">
  <image width="700" height="526" src="https://github.com/nielsthaumann/nameeg/blob/main/name_ons_03.png">
</p>

The automatic onset detection is performed using finetuned functions
from the Music Information Retrieval
(MIR) toolbox for Matlab ([Lartillot &
Toiviainen, 2007](https://ismir2007.ismir.net/proceedings/ISMIR2007_p127_lartillot.pdf)). <br></br>
The mirspectrum function is used with
a 100 ms frame duration and 10 ms hop size, applying the
'Blackmann-Harris’ window to suppress scalloping loss in
the Fourier transform. <br></br>
The ‘Terhardt’ filter ([Terhardt, 1979](https://doi.org/10.1016/0378-5955(79)90025-x);
[Pampalk et al., 2004](https://doi.org/10.1162/014892604323112248)) 
is applied to simulate the effects of 
the human outer ear on the perceived spectra. <br></br>
Spectro-temporal amplitude changes are computed using
the mirflux function, incorporating the ‘Emerge’ filter
([Lartillot et al., 2013](https://www.academia.edu/download/121050832/olivier_lartillot_-_estimating_tempo_and_metrical_features_by_tracking_the_whole_metrical_hierarchy.pdf)) 
to minimize the influence of tremolo and vibrato. <br></br>
Sound onsets are then identified at peak values in the spectral 
flux using the mirevents function with the default settings.

Please make sure that a recent version of the MIRtoolbox is installed. 
(This function was tested with MIRtoolbox 1.8.1.)<br></br>
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

<p align="center">
  <image width="960" height="513" src="https://github.com/nielsthaumann/nameeg/blob/main/name_onsed_01.png">
</p>

