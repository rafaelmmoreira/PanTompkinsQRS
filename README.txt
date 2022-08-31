An ANSI-C implementation of the Pan-Tompkins Real-Time QRS Detection Algorithm
With Embedded Version
Author: Rafael de Moura Moreira <rafaelmmoreira@gmail.com>
        https://github.com/rafaelmmoreira/PanTompkinsQRS
		
		Thorold Tronrud <ttronrud@starfishmedical.com>
		https://github.com/ttronrud/PanTompkinsQRS
License: MIT License
Copyright (C) 2018-2022

USING THE ALGORITHM
Just import the .c and .h to your project, or paste them in the same folder and include "panTompkins.h".
To use the algorithm "as is", you must first call the init() function passing 2 arguments: the name of
your input file (which must be a list of integers in ASCII) and the name of your output file (be careful,
it's an existing file, it will be overwritten!).
It will output a list of 0's and 1's, where 0 means a given sample didn't trigger a R-peak detection,
while an 1 means it did.

EMBEDDED ALGORITHM
The embedded version of the algorithm has various components shifted around for use in more memory-limited systems,
where you don't want an extensive buffer of input data. The filters have been moved globally, so the learning
is retained in the state of the system, instead of being re-acquired each iteration. This is demonstrated with the
two "init" functions, which swap the signal used as input. The two arrays are separate segments of the same recording
and once the R peaks have been detected in one, the other is loaded. The system provides accurate peak locations
immediately.

MODIFYING THE CODE
The code was designed to be easy to change and port: your input source (file, serial comms etc), the
input format (signed or unsigned int, float, double etc), sampling frequency, fine-tunings to the algorithm, 
and so on.
The .c file is well documented, with every meaningful line or chunk of code explained in the comments,
besides a long description which suggests all the pertinent changes to make it work on different applications
and systems.

TESTING
The code, "as is", should work on Windows and Linux for x86. A test input file (the lead A for patient 100 from
the MIT-BIH database converted to ASCII) and the output for this signal are included in the examples folder.
There's also a plot of the first 10000 samples (~27.8 seconds) plus the output. One can note on this plot two
limitations from the algorithm: it takes about 2 R-R intervals to learn (before that, its thresholds are still
adjusting and there are a couple of false positives). After the first 2 seconds, the algorithm stabilizes. For 
patients with anomalous ECG signals, chances of false positives or missed detections increase. However, this 
algorithm is known for a very high precision. 
Also, the output is delayed by a few milisseconds due to the filtering stages. A fix has been added by ignoring
the first few samples so that the input and output signals' peaks match one another. The down side is missing a
few samples.
I've also added (April 2019) another plot showing a few heartbeats from the original signal plus the output
from the bandpass filter, the derivative, the squared derivative, the integral and the Pan-Tompkins classification.

For further information about the data used on the test:
https://www.physionet.org/physiobank/database/mitdb/  
Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. 
IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)

For further information about the algorithm, its details and its limitations:
Pan, J., & Tompkins, W. J. (1985). A real-time QRS detection algorithm. 
IEEE transactions on biomedical engineering, (3), 230-236.





