# omura-jasa-2020 [![DOI](https://zenodo.org/badge/249245872.svg)](https://zenodo.org/badge/latestdoi/249245872)
Implementation of algorithms described in

> Shyam Madhusudhana, Anita Murray, and Christine Erbe. "Automatic detectors for low-frequency vocalizations of Omura’s whales, _Balaenoptera omurai_: A performance comparison." _The Journal of the Acoustical Society of America_ (2020)

for the detection of Omura's whale vocalizations.

Of the three methods described in the article, implementations of the following two are provided here -
* Blob Detector (implemented in Python)
* Entropy Detector (implemented in Matlab)

Instructions for using them are provided below. If you use these (in entirety or as parts), we request that you please cite the aforementioned article.

---

## Blob Detector
The Python3 implementation provided in [blob_detector/omura_detector.py](blob_detector/omura_detector.py) requires the Python packages [DetectSound](https://github.com/shyamblast/DetectSound) (which contains the general-purpose blob detection functionality), [librosa](https://librosa.github.io/librosa/) (for loading audio files) and [matplotlib](https://matplotlib.org/) (used for displaying detection results). These can be installed from command prompt as -
```shell
$> pip3 install librosa matplotlib
$> pip3 install git+https://github.com/shyamblast/DetectSound.git
```
After downloading **omura_detector.py** on your local machine, the detector can then be run as
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
The Matlab implementation of the entropy detector is provided as a collection of 3 files in [entropy_detector](entropy_detector). After downloading the entire folder, set directory paths appropriately on lines 33 & 34 in [OmuraEntropy2XSTD.m](entropy_detector/OmuraEntropy2XSTD.m) and then execute the script within Matlab.
