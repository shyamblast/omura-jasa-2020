# omura-jasa-2020
Implementation of algorithms described in

> Shyam Madhusudhana, Anita Murray, and Christine Erbe. "Automatic detectors for low-frequency vocalizations of Omura’s whales, Balaenoptera omurai: A performance comparison." _The Journal of the Acoustical Society of America_ (2020)

for the detection of Omura's whale vocalizations.

Of the three methods described in the article, implementations of the following two are provided here -
* Blob Detector (implemented in Python)
* Entropy Detector (implemented in Matlab)

Instructions for using them are provided below. If you use these (in entirety or as parts), we request that you please cite the aforementioned article.

---

## Blob Detector
The Python3 implementation provided in **omura_detector.py** requires the [DetectSound](https://github.com/shyamblast/DetectSound) module (which contains the general-purpose blob detection functionality), [librosa](https://librosa.github.io/librosa/) (for loading audio files) and [matplotlib](https://matplotlib.org/) (used for displaying detection results). These can be installed from your command prompt as -
```shell
$> pip3 install librosa matplotlib
$> pip3 install git+https://github.com/shyamblast/DetectSound.git
```
The detector can then be run as
```shell
$> python3 omura_detector.py /path/to/audio_file.wav
```
to display the detection results on screen, or as
```shell
$> python3 omura_detector.py /path/to/audio_file.wav /save/to/output_file.txt
```
to save the detection results to a file in Raven Selection Table format.

---

## Entropy Detector
The Matlab implementation provided in **omura_detector.m** requires the following Matlab toolboxes:
* signalprocessing
* any others?

It can be run as ...
```matlab
> ...
```
