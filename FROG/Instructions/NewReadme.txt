On MATLAB versions: we have tested this code with MATLAB 2012a and 2011b on Windows(64) and 2010b on Linux(64) and 2009b on Mac(64).
2008a probably works too.  Other versions might work, but we haven't checked.

***********SETUP*****************************************

1. Download 'FROG_v1.2.0.7z'. This contains all of the code, both the MATLAB and C++ source code.

2. Download 'Window_64bit.7z' or, if you're not using 64 bit Windows, the file that matches your operating system.  This contains precompiled C++ code for the library functions.  
If there isn't a version of the precompiled code that matches your operating system, you will have to compile the code yourself.  Skip down to the section on 'Compiling the Library' for instructions on how to do this.

3. Unzip them - you will need a program like 7Zip or Winrar to do this.

4. Make sure both folders you just downloaded (and all their subfolders) are added to your MATLAB path.

5. Try running FROG_v1.2.0/Test_QuickFROG.m .  If this works, everything should be successfully set up!
   If you have issues, check your MATLAB version and your operating system version, and double-check your MATLAB path information.

We provide a GUI version and a non-GUI version of the code.  Check out the instructions for the version that most interests you!  If you're curious about the inner workings or have lots and lots of data, we suggest the non-GUI code.  If you already feel a bit out of your element or don't have much data, try the GUI code.




***********GUI CODE INSTRUCTIONS**************************

Overview:
 Calibration
 binner
 frogger

1. Calibration 

The binner program expects a certain file format, and the main purpose of this step is to create the correct file format.  This file format is:

[width of trace (pixels)] \t [height of trace (pixels)] \t [temporal calibration (fs/pixel)] \t [spectral calib. (nm/pixel)] \t [center wavelength] \n
[raw data matrix]

IMPORTANT: The width should be the time axis, and the height should be the wavelength axis.
Save the resulting trace with the format given above as a .frg file.

2. binner

Run binner in MATLAB by typing >> binner;  A GUI should come out, open the .frg file you prepared. A pop-up will appear and ask if "Delay" or "Wave" is the first axis. If you prepared the file according to step 1, you should choose "Wave". Make sure the FROG trace displayed looks reasonable. Extraction, Background Subtraction, Centering and Filtering can be performed using different tabs. 

Binning tab
-Array Size: Depending on the complexity of the pulse, an appropiate size should be chosen. For pulses with short temporal range and narrow spectral bandwidth, a smaller size can be used. Since the Fourier transform relation has to be satisfied, the temporal range is proportional to the spectral resolution and the spectral range is proportional to the temporal resolution. 
-Binned Width (%): The program allow users to pad zeros or remove empty space of the traces to increase or decrease of the resolution and range. For example, binned width of 50% means fill the binned trace with 50% of the raw trace and pad the rest with zeros. 
-Axis Fit: Use along with Binned Width, fit delay means using the temporal domain to consider the binned width. 
-Spec. Eff. Corr. and Freq. Marginal. Corr.: If a wide range of spectrum is presented in the trace, the response of the grating, camera and lens could play a part to the intensity. These two functions allow users to correct for the response of the optics component. 
-Bin Data: Once everything is set, click this button to bin the trace. The binned trace will be displayed and if anything goes wrong, click the 'Undo' button to restore. A good binned trace should have energy in most of the area. 

After binning the trace, save the binned trace. The program should automatically change the file to xxx.bin.frg. 

3. frogger

Make sure you have a binned trace. You can use the 'binner' GUI or command version of binner, 'binner_cmd' (discussed in detail in the non-GUI code section)

Run in MATLAB >> frogger;   A GUI should come out, open the binned trace you prepared (should have extension .bin.frg if prepared by binner). A pop-up will appear and ask if "Delay" or "Wave" is the first axis, if you prepare the file using GUI 'binner', you should choose "Delay". 

Choose the 'Nonlinearity' according to your experimental setup. The current version of frogger provides 'SHG', 'PG' and 'THG'.

Click 'Run' button'. The program should start running, stop when the desired error is reached. You can save the retrieved E-field in temporal and spectral domain and also the measured and retrieved trace using 'save' function. The saved files have the following structure: 

=== file structure of Ek.dat and Speck.dat ===
x=1\t Intensity\t Phase\t Real()\t Imaginary()\t\n
:
:
:
x=N\t Intensity\t Phase\t Real()\t Imaginary()\t\n
=== EOF ===
where x:= time in Ek.dat and x:= lambda in Speck.dat

To use XFROG, switch the type to XFROG. Click 'Extra Information' to choose the gate pulse. You can choose from a computer generated function like sech and Gaussin or read the gate from a file. Check 'use file' inside 'extra information' to use a measured gate pulse. If you use frogger to retrieve the gate pulse, choose the file contains the temporal E-field, which by default is Ek.dat. If you obtain the gate pulse from other sources, construct a file where the first column is the x-axis and the 4th and 5th columns are the real and imaginary part of the temporal E-field. 





***********NON-GUI CODE INSTRUCTIONS**********************

We try to provide everything you need to go from a raw camera FROG trace to a retrieved pulse.  If you are dealing with lots of data, you may want to try to automate some of the processing.  An example of how to do this is PG_XFROG_auto.m .

Overview:
1. Calibration / creating a .frg file
2. Binning
3. Retrieval

1. Calibration 

The binner program expects a certain file format, and the main purpose of this step is to create the correct file format.  This file format is:

[width of trace (pixels)] \t [height of trace (pixels)] \t [temporal calibration (fs/pixel)] \t [spectral calib. (nm/pixel)] \t [center wavelength] \n
[raw data matrix]

IMPORTANT: The width should be the time axis, and the height should be the wavelength axis.

This isn't always necessary, but if you're doing any background subtraction or correction for spectral response or other calibration-type tasks, this is probably a good place to put that.  Save the resulting trace with the format given above as a .frg file.

Note for more advanced users: you can look at calibrate.m (in the demo folder) for an example of how to automate this.  Most of the code in there is for manipulating file names.  Given a raw trace file name, that program automatically locates a background measurement (this code was used for PG FROG and needs to subtract the polarizer leakage) associated with the given trace and subtracts the background.


2. Binning

The FROG code uses a lot of fast Fourier transforms, and consequently needs the time and frequency/wavelength axes to be the same size and also a power of 2.  Also, the FROG algorithm gets slower the larger the trace is, so it's a good idea to use a smaller trace if and only if your FROG trace is simple.

binner_cmd.m accomplishes two things: converting the data from wavelength to frequency, and putting it on the specified grid size using 2D interpolation.  When you call binner_cmd, give it the name of your .frg file.  You can also specify the number of pixels to use by passing it ([filename], 'sz', [size]).  The default is 128x128.

By default, the full delay range in your data (calculated from the delay calibration times the number of pixels in time dimension) will be put into half the delay range of the binned trace.  To change this, you can pass a different percent of the time data to put into half the binned range to binner_cmd by adding the arguments 'width', [width] to the function call.  Remember, though, that the trace really needs to go to zero at the edges. 

The frequency range and resolution is then calculated based on the time range and resolution, using the Fourier transform relations.  If your data doesn't fit on the resulting axis, you will have to either a) change the width of the time data, or b) fit the frequency dimension first.  To fit frequency first, add the arguments 'method', 2 to the binner_cmd call.

binner_cmd will pass back the file name of the binned data file.

If your trace does not go to zero at the edges, the retrieval algorithm will not work well.  


3. Retrieval

The type of nonlinearity used to generate the data is part of the FROG algorithm, so you will have to call the correct version of the FROG retrieval for the nonlinearity in your experiment.  We currently provide code for SHG, PG, and THG (as well as XFROG using any nonlinearity).  FROG_v1.2.0/lib/frog/QuickFrog_tT.m is set up by default to do SHG FROG.  To change this, 

-for PG FROG:
change line 61: Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,abs(Et).^2);
change line 81: Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,abs(Et).^2);
change line 122:Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,abs(Et).^2);

change line 72: dZ = -dZdE_shg(Esig,Et); ->  dZ = -dZdE_pg(Esig,Et);
change line 118:dZ = -dZdE_shg(Esig,Et); ->  dZ = -dZdE_pg(Esig,Et);

change line 74: [Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);  ->  [Et, Z(k)] = MinZerr_pg(Esig, Et, dZ);
change line 120:[Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);  ->  [Et, Z(k)] = MinZerr_pg(Esig, Et, dZ);


-similarly, for THG FROG:
change line 61: Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,Et.^2);
change line 81: Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,Et.^2);
change line 122:Esig = CalcEsig(Et,Et);  ->  Esig = CalcEsig(Et,Et.^2);

change line 72: dZ = -dZdE_shg(Esig,Et); ->  dZ = -dZdE_thg(Esig,Et);
change line 118:dZ = -dZdE_shg(Esig,Et); ->  dZ = -dZdE_thg(Esig,Et);

change line 74: [Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);  ->  [Et, Z(k)] = MinZerr_thg(Esig, Et, dZ);
change line 120:[Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);  ->  [Et, Z(k)] = MinZerr_thg(Esig, Et, dZ);


-BUT for ANY version of XFROG, use FROG_v1.2.0/demo/qFROG_TX.m .  Be sure to pass to the method the correct form of the gate pulse depending on the nonlinearity (i.e. abs(E_gate).^2 for PG).




***********COMPILING THE LIBRARY**************************
Setup MEX: run mex -setup, follow the instruction for a compiler. If a compiler is already installed, Matlab should be able to locate it and use it. Otherwise, follow the instructions from Matlab to download and install a compatible compiler. Details can be found here http://www.mathworks.com/support/compilers/R2011a/win64.html or similar link according to Matlab and OS version. 

Compile: once the MEX is setup correctly, compile the MEX file. For example, compile CalcEsig.cpp by using the command >> mex CalcEsig.cpp. Please be aware that different OS and Matlab version gives different file extension, and MEX file compile on a Mac is not compatible with the one compiled on Linux or Windows.  Copy the compiled mex file to a folder which is under the Matlab path. Make sure the path contain all the mex files is at the highest priority since a helper file with the same filename but with extension .m will also be present. 

Batch compile: batch_mex.m is include in the source directory (/lib/Mexlib/src). Simply run it to compile all the cpp file to generate all the MEX file needed.  

* An issue has been found in MagRepl.cpp with standard gcc library in Linux and Mac. If you are using some other C compilers, please be aware that similar kinds of issues may also occur. 