#SmarTrace Manual

Chain tracing program developed to analyze AFM images of filaments such as collagen and DNA. Can be used on any type of image with bright chain-like structures on a darker background. 

Software was developed by Naghmeh Rezaei while a PhD student at Simon Fraser University.  The program is described in her PhD thesis (2016, Simon Fraser University).

More information about what SmarTrace can be used for is described in the following publication: Rezaei, N., Lyons, A. & Forde, N. R. Environmentally controlled curvature of single collagen proteins. Biophys. J. 115, 1457–1469 (2018).

Graphical User Interface is adapted from the EasyWorm GUI: G. Lamour, J.B. Kirkegaard, H. Li, T.P.J. Knowles & J. Gsponer, Source Code for Biology and Medicine 9,16 (2014).

**BEFORE USE**

Use of SmarTrace requires a MATLAB distribution of 2018a or later for the Deep Learning Toolbox.

Images to be analyzed should be in the same directory.


**Installation**

Move the SmarTrace folder into a directory of your choice, then save the directory to the MATLAB path.

ex: |savepath('C:\Program Files\matlab\SmarTrace');
    |savepath('C:\Project\Chain Tracing\SmarTrace');

The following MATLAB toolboxes are required:
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox


**Using SmarTrace**

1- Change into the directory that contains your images, then hit Run.

2- The SmarTrace GUI window will appear, where you can "Select Height Map" in the top right-hand corner. 

3- Once you load the image of interest, enter the image resolution: nm per px. 

4- Click on "Select Chain" and use your cursor to define a few inputs on the chain you want to trace. 

5- You can press "Z" to zoom into the chain of interest, 'A' to zoom out, and "D" to delete the unwanted selected points.

6- Once you have a few user-defined inputs on the chain, click "enter" and wait for SmarTrace to define the backbone. After the analysis a red curve will be plotted on the chain. You will also retrieve a contour length of your chain in the Matlab Command Window. 

7- If you are happy with the traced chain, click "Add Chain" to add the data of the chain and click "+" to increase the chain number by 1 unit. This will help you to keep track of the number of chains traced in your data set. 

8- You may feel free to save the traced images using "Save Image", and save the data containing traced chains using "Save Data". The saved ".mat" file is ready for analysis.  

9- If you wish to stop tracing, and want to continue later on, you may use "Reload Data" to open the data set and continue adding chains to the file. 
