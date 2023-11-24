The following program takes data from roughness measurements using a confocal microscopy and transforms it into 2 dimensional Power Spectral Distribution plots.

To use this program the following Python 3.11 libraries are necessary: numpy, scipy and matplotlib. The use of Anaconda and VSCode is recommended, and makes the process much faster.

## Installation
1. Clone this directory into a desired folder.
2. Install the necessary libraries or Anaconda.

## How to setup before running

For this script to run smoothly it is important that the data associated with the files is well separated into individual folders with the name of the samples. For example:
![[Pasted image 20231124161945.png]]
**Where each folder has inside one and only one confocal data of the sample at a fixed magnification.** 

Now open the program 'Analysis_2D.py' and navigate to the code snippet between lines 173 and 178:
![[Pasted image 20231124162301.png]]
In here you should place the path and name of the results folder, the folder where the output of the program will be and the path to where you have your data. The magnification should also be explicit and **constant through all samples**. Add the names of the samples to the 'sample_list' and make sure that the folders where you store the confocal data have the same names as the samples. Let the variable 'PLACE' as it is.
## How to run
1. From a terminal window identify your python env and run using `& C:/Users/YOURUSER/AppData/Local/anaconda3/python.exe "PATH_TO_SCRIPT/Confocal2PSD/Analysis_2D.py"`
2. This will go over each one of your folders, run and script and save the data.


Any questions contact tomas.sousa@unibas.ch
