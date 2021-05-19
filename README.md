# SmarTrace
Chain tracing program developed to analyze AFM images of filaments such as collagen and DNA. Can be used on any type of image with bright chain-like structures on a darker background. 

Software was developed by Naghmeh Rezaei while a PhD student at Simon Fraser University.  The program is described in her PhD thesis (2016, Simon Fraser University).

More information about SmarTrace can found in the following publication: N. Rezaei, A. Lyons, & N.R. Forde, Environmentally controlled curvature of single collagen proteins. Biophys. J. 115, 1457–1469 (2018).

Graphical User Interface is adapted from the EasyWorm GUI: G. Lamour, J.B. Kirkegaard, H. Li, T.P.J. Knowles & J. Gsponer, Source Code for Biology and Medicine 9,16 (2014).

**BEFORE USE**

Use of SmarTrace requires MATLAB.

Images to be analyzed should be in the same directory.

To use SmarTrace, you need to download DIPlib, which can be found on GitHub. 

**Installation**

Move the SmarTrace folder into a directory of your choice, then save the directory to the MATLAB path.

ex: |savepath('C:\Program Files\matlab\SmarTrace');
    |savepath('C:\Project\Chain Tracing\SmarTrace');


**Using SmarTrace**

1- To use SmarTrace, you will need to initialise DIPimage first:

|addpath('C:\...\DIPimage\common\dipimage')
|dip_initialise

2- Change into the directory that contains your images, then hit Run.

3- The SmarTrace GUI window will appear, where you can "Select Height Map" in the top right-hand corner. 

4- Once you load the image of interest, enter the image resolution: nm per px. 

5- Click on "Select Chain" and use your cursor to define a few inputs on the chain you want to trace. 

6- You can press "Z" to zoom into the chain of interest to make it easier to trace, and "D" to delete the already selected points and re-select them.

7- Once you have a few user-defined inputs on the chain, click "enter" and wait for SmarTrace to define the backbone. After the analysis a red curve will be plotted on the chain. You will also retrieve a contour length of your chain in the Matlab Command Window. 

8- If you are happy with the traced chain, click "Add Chain" to add the data of the chain and click "+" to increase the chain number by 1 unit. This will help you to keep track of the number of chains traced in your data set. 

9- You may feel free to save the traced images using "Save Image", and save the data containing traced chains using "Save Data". The saved ".mat" file is ready for analysis.  

10- If you wish to stop tracing, and want to continue later on, you may use "Reload Data" to open the data set and continue adding chains to the file. 
