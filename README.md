MC Metalicity
====================

Monte Carlo method to calculate metalicity uncertainty from flux data.

====================
Usage:
====================
Place the  _meas and _err  files in the "input" directory (example files are provded in the directory).


From the commandline simply use as:
```
fedMCZ_err.py <filename> nsample --path PATH --clobber --delog --verbose
```
-\<name\>: the SN name which should be the common name of the _min and _max files (e.g. testdata13 for testdata13_meas.txt and testdata13_err.txt)

-nsample: the number of MC samples desired 

--path: the directory in which subdirectory "input" is located. If not provided, will default to environmental variable MCMetdata that should be set to point to that directory. _err.txt _meas.txt must live in <path>/input

--unpickle: if it exists, a pickle output files for this SN and this number of samples is read, instead of recalculating the metalicity

--binmode: how to choose the number of bins for plotting the histogram:'d' is based on Doane's formula (wilipedia's version), 's' is the sqrt of number of data, 't' is on 2*n**1/3 (default), 'kb' uses Knuth's block rule, 'bb' uses bayesian blocks (must have astroML installed or it defaults to 't')

--clobber: If set to true, will overwrite existing output files. Default False.

--delog: If set to true, result will be in natural space instead of log space. Default False.

--verbose: verbose mode. Default False.


====================
Input file format
====================
each flux data should be stored in the directory sn_data that exists in the directory provided by --path arg. 

with common filename \<fi\>:

\<fi\>_err.txt

\<fi\>_meas.txt 


where max = mesured+err, min=mesured-err. The _mes file is optional.

The format for each of the txt files should be as follows:


;# galnum,[OII]3727,Hb,[OIII]4959,[OIII]5007,[OI]6300,Ha,[NII]6584,[SII]6717,[SII]6731,[SIII]9069,[SIII]9532,E(B-V),dE(B-V)
       1     0.0     0.0     0.0     0.0     0.0   5.117   0.998     0.0     0.0     0.0     0.0
       2     0.0     0.0     0.0     0.0     0.0   5.031   1.012     0.0     0.0     0.0     0.0
       
       
galnum is the index used for matching the flux with the corresponding radius for calculating the gradient.

The first row is optional - if can contain more keys than column, as long as the columns are listed in the same order as specified above, up to whichever line is of interest. Missing data should be filled in with 0s.

the data should be in the above order, separated by any number of white spaces.


====================
Output
====================
All the results will be saved in the directory "\<fi\>" inside "outputs" which will be creates in the directory provided by the --path arg.

For a given \<fi\> and nSample, the following output files are generated:

"\<fi\>_n"nSample"_X.csv": the metallicity and its uncertainty that was calculated, where X will correspond to the number of rows the input data contains

"\<fi\>_n"nSample"_X.pkl": the metallicity and its uncertainty that was calculated stored in binary (pickled) format in a python dictionary

"hist" folder, containing all the histograms generated


 
