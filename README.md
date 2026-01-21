# Naturalistic auditory MEG/EEG (NAME) package

Auditory neuroscience using MEG/EEG methods has primarily been conducted in controlled laboratory settings, utilizing simplified stimuli, such as pure sinusoidal tones or brief sound sequences. However, recent research questions to which extent the human brain responds to simplified stimuli the same way as to sounds heard in the real world in everyday situations. For example, recent research shows that different neural populations are involved in the processing of naturalistic than simplified sound sequences ([Haumann et al., 2023](https://doi.org/10.1016/j.biopsycho.2023.108566); [Khalighinejad et al., 2021](http://dx.doi.org/10.1016/j.neuroimage.2021.118003)). 

<p align="center">
  <image width="750" height="358" src=https://ars.els-cdn.com/content/image/1-s2.0-S1053811921002809-gr1_lrg.jpg>
</p>
<p align="center">
  <a href=https://doi.org/10.1016/j.neuroimage.2021.118003> (Khalighinejad et al., 2021. "Functional characterization of human Heschl's gyrus in response to natural speech." <i> Neuroimage</i>, 235, 118003.)
  </a>
</p>
</br>

A popular assumption is that a one-to-one copy of the rich audio details in the acoustic environment is encoded in the human brain's cortex, which can be recorded using MEG/EEG and directly converted back into audio waveforms. However, the auditory cortex rather encodes interpretations of relevant sound events, not the rich audio details in the acoustic environment for making the interpretations of the sound events ([Haumann et al., 2021](https://doi.org/10.1016/j.brainres.2020.147248); [Haumann et al., 2018](https://doi.org/10.3390/app8050716); [Poikonen et al. 2016](https://doi.org/10.1016/j.neuroscience.2015.10.061)).

<p align="center">
  <image width="700" height="700" src="https://github.com/user-attachments/assets/a90f1238-6cd4-43d3-b54a-a65faf2e12fb">
</p>
<p align="center">
  <a href=https://doi.org/10.1016/j.brainres.2020.147248> (Haumann et al., 2021, "Extracting human cortical responses to sound onsets <br> and acoustic feature changes in real music, and their relation to event rate",<i> Brain Research</i>)</br>
  </a>
</br>

An analogue to information encoding in the cortex is the cave metaphor by Plato, which illustrates how our interpretation of what happens in our sorroundings is a simplified, but often useful, understanding of what actually happens in reality. </br>
E.g., we see a horse, not the actions and light reflections that made us draw the interpretation that there is a horse:
<p align="center">
  <image width="850" height="300" src=https://upload.wikimedia.org/wikipedia/commons/thumb/8/8d/An_Illustration_of_The_Allegory_of_the_Cave%2C_from_Plato%E2%80%99s_Republic.jpg/1920px-An_Illustration_of_The_Allegory_of_the_Cave%2C_from_Plato%E2%80%99s_Republic.jpg>
</p>
<p align="center">
  <a href=https://commons.wikimedia.org/wiki/File:An_Illustration_of_The_Allegory_of_the_Cave,_from_Plato%E2%80%99s_Republic.jpg> (4edges, Wikimedia Commons)
  </a>
</p>
</br>
Similarly, we hear a sound, not the continuous acoustical changes that made us draw the interpretation that there is a sound: 
<p align="center">
  <image width="1127" height="766" src=https://github.com/nielsthaumann/nameeg/blob/main/event_related_vs_TRF.png>
</p>
<p align="center"><i> The above test results show that cortical brain responses are predicted by sound onset events, but not from continuous changes in spectral flux or sound intensity when excluding the rapid increases correlating with the sound onsets. </i>
</p>
    
<br>
</br>

Please cite this site ( https://github.com/nielsthaumann/nameeg ) if you are publishing results using functions from this package. 

<br>
</br>

----------------------------------------
## Automatic sound onset detection with noise suppression

Use the onset detection function with the following command: 

onsets = [name_ons](https://github.com/nielsthaumann/nameeg/blob/main/name_ons.m)(audio)  

to estimate sound onsets in an audio file, where, e.g., audio = 'C:\folder\audio.wav'

The onsets are given as an array of time points in seconds.

The automatic sound onset detection algorithm includes suppression of noise in the audio medium (see [name_ns](https://github.com/nielsthaumann/nameeg/blob/main/name_ns.m)). 
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
  <image width="700" height="526" src="https://github.com/nielsthaumann/nameeg/blob/main/name_ons_01.png">
</p>

Brownian noise affecting low frequency bands in the audio
is suppressed with the non-negative least squares
function (lsqnonneg) to fit a 1/f curve to the average
amplitude spectrum. 

<p align="center">
  <image width="240" height="240" src="https://upload.wikimedia.org/wikipedia/commons/c/c2/Brownian_motion_large.gif">
</p>
<p align="center">
  <a href=https://da.m.wikipedia.org/wiki/Fil:Brownian_motion_large.gif> (Lookang, Brownian motion in air, Wikimedia Commons)
  </a>
</p>
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
the human outer ear on the perceived spectra. 
<br></br>
Spectro-temporal amplitude changes are computed using
the mirflux function, incorporating the ‘Emerge’ filter
([Lartillot et al., 2013](https://www.academia.edu/download/121050832/olivier_lartillot_-_estimating_tempo_and_metrical_features_by_tracking_the_whole_metrical_hierarchy.pdf)) 
to minimize the influence of tremolo and vibrato. <br></br>
Sound onsets are then identified at peak values in the spectral 
flux using the mirevents function with the default settings.

Please make sure that a recent version of the MIRtoolbox is installed. 
(This function was tested with MIRtoolbox 1.8.1.)<br></br>
E.g, visit: [https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mirtoolbox](https://github.com/olivierlar/mirtoolbox)

----------------------------------------
## Sound onset editor

Use this editor to combine the speed of automatic sound onset detection 
with the reliability of manual editing.

Call the editor from the Matlab Command Window and hit Enter: </br>
[name_onsed](https://github.com/nielsthaumann/nameeg/blob/main/name_onsed.m)

Load an audio file (.wav, .flac, .mp3, .m4a, .mp4, .ogg) for auditory validation
and visual inspection of energy increases in a time-frequency representation. 

Load and save a .mat file with the variable 'onsets' containing a numerical vector 
with the onsets stored in seconds. 

Changes are saved automatically (along with the undo/redo history). 

<p align="center">
  <image width="960" height="513" src="https://github.com/nielsthaumann/nameeg/blob/main/name_onsed_01.png">
</p>



## For further information, please refer to: 

Haumann, N. T., Petersen, B., Seeberg, A. B., Vuust, P., & Brattico, E. (2025). <b> "It takes experience to tango: Experienced cochlear implant users show cortical evoked potentials to naturalistic music."</b> bioRxiv, 2025.2006.2004.657805. https://doi.org/10.1101/2025.06.04.657805

Haumann, N. T., Petersen, B., Vuust, P., & Brattico, E. (2023).  <b> "Age differences in central auditory system responses to naturalistic music."</b> Biological Psychology, 179, 108566. https://doi.org/10.1016/j.biopsycho.2023.108566

Haumann, N. T., Lumaca, M., Kliuchko, M., Santacruz, J. L., Vuust, P., & Brattico, E. (2021). <b> "Extracting human cortical responses to sound onsets and acoustic feature changes in real music, and their relation to event rate". </b> Brain Research, 1754, 147248. https://doi.org/10.1016/j.brainres.2020.147248

Haumann, N. T., Kliuchko, M., Vuust, P., & Brattico, E. (2018).  <b> "Applying acoustical and musicological analysis to detect brain responses to realistic music: A case study." </b> Applied sciences, 8(5), 716. https://doi.org/10.3390/app8050716

