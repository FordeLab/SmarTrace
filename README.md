#SmarTrace Manual

Chain tracing program developed to analyze and extract persistence lengths from AFM images of filaments such as collagen and DNA. Can be used on any type of image with bright chain-like structures on a darker background. 

Software was developed by Naghmeh Rezaei while a PhD student at Simon Fraser University.  The program is described in her PhD thesis (2016, Simon Fraser University).

More information about what SmarTrace can be used for is described in the following publication: Rezaei, N., Lyons, A. & Forde, N. R. Environmentally controlled curvature of single collagen proteins. Biophys. J. 115, 1457–1469 (2018).

Graphical User Interface is adapted from the EasyWorm GUI: G. Lamour, J.B. Kirkegaard, H. Li, T.P.J. Knowles & J. Gsponer, Source Code for Biology and Medicine 9,16 (2014).

**Before Use**

Use of SmarTrace requires a MATLAB distribution of 2018a or later.

The following Add-Ons are required to use MatLab: 
	- Deep Learning Toolbox
	- Curve Fitting Toolbox 
	- Computer Vision Toolbox
	- Signal Processing Toolbox
	- Statistics and Machine Learning Toolbox 
	- Image processing Toolbox

These toolboxes can be added during the process of installing MatLab as well as afterwards by going to the 'Home' tab, clicking on 'Add-Ons' and searching for the appropriate add-ons.


**Installation**

Move the SmarTrace folder into a directory of your choice, then save the directory to the MATLAB path.

ex: |savepath('C:\Program Files\matlab\SmarTrace');
    |savepath('C:\Project\Chain Tracing\SmarTrace');


**Using SmarTrace**
Description updated by Alyssa Abake Oke, July 2021

1- Open the file 'SmarTrace.m' from the folder 'SmarTrace-main\SmarTrace\SmarTrace\ChainTracing'.

2- Change the 'Current Folder' on the left hand side to the folder that contains your images. 

3- Click 'Run' under the 'Editor' tab. 

4- The SmarTrace GUI window will appear, where you can 'Select Height Map' in the top right-hand corner. 

5- Once you load the image of interest, enter the image resolution: nm per px. 

6- Click on 'Select Chain' and use your cursor to define a few input points on the chain you want to trace. 

7- You can press 'Z' to zoom into the chain of interest, 'A' to zoom out, and 'D' to delete the unwanted selected points.

8- Once you have a few user-defined inputs on the chain, click 'enter' and wait for SmarTrace to define the backbone. After the analysis a red curve will be plotted on the chain. You can also retrieve a contour length of your chain in the Matlab Command Window. 

9- If you are happy with the traced chain: 
	a. Change the 'set_sample_name', 
	b. Click 'Add Chain' to add the data of the chain 
	c. and click '+' to increase the chain number by 1 unit. This will help you to keep track of the number of chains traced in your data set. 

10- Repeat Steps 6-9 for the following chains. You may also select a new height map at this point before hitting 'Select Chain' if you have several images with chains you want to include in your analysis.

11- You may wish to save the traced images using 'Save Image', and save the data containing traced chains using "Save Data". The saved ".mat" file is ready for analysis.  

12- If you wish to stop tracing, and want to continue later on, you may use "Reload Data" to open the data set and continue adding chains to the file. 


**Data Analysis**

1- When you are done tracing and your last chain is saved, exit out of the SmarTrace window. There will be a saved data file in the images' folder and type the following in the command window:

	analyseMOD('set_sample_name_SmTr.mat', 'output_name_of_choice',10,0,0) --> For the WLC model analysis

	OR

	analyseMOD('set_sample_name_SmTr.mat', 'output_name_of_choice',10,0,1) --> For the curved WLC model analysis

2- Press 'enter'

3- View the output folder with your output name of choice in the images' directory.
