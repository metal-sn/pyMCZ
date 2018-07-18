[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.17880.svg)](http://dx.doi.org/10.5281/zenodo.17880) 
[![Code Health](https://landscape.io/github/fedhere/pyMCZ/master/landscape.svg?style=flat)](https://landscape.io/github/fedhere/pyMCZ/master)
[![ASCL](ascl-1505.025-blue.png)](http://ascl.net/1505.025)
MC Metallicity
====================

This code allows the user to calculate metallicity according to a number of <i> strong line metallicity diagnostics </i> from spectroscopy line measurements, and obtain uncertainties from the line flux errors in a Monte Carlo framework. If you use this code, <b> please cite our [paper](http://www.sciencedirect.com/science/article/pii/S2213133716300178) which is published in [Astronomy and Computing](http://www.journals.elsevier.com/astronomy-and-computing/) and also available on the [arxiv](http://arxiv.org/abs/1505.06213).</b>

    @ARTICLE{2016A&C....16...54B,
       author = {{Bianco}, F.~B. and {Modjaz}, M. and {Oh}, S.~M. and {Fierroz}, D. and 
	    {Liu}, Y.~Q. and {Kewley}, L. and {Graur}, O.},
        title = "{Monte Carlo method for calculating oxygen abundances and their uncertainties from strong-line flux measurements}",
      journal = {Astronomy and Computing},
    archivePrefix = "arXiv",
      eprint = {1505.06213},
    primaryClass = "astro-ph.IM",
    keywords = {Galaxy, Abundances-ISM, HII regions-supernovae, General},
         year = 2016,
        month = jul,
      volume = 16,
        pages = {54-66},
          doi = {10.1016/j.ascom.2016.03.002},
      adsurl = {http://adsabs.harvard.edu/abs/2016A%26C....16...54B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }



This code is released under MIT licence: see LICENSE.txt

<b>Most recent stable release: v1.3</b>


Usage:
====================
This code requires ascii input files, their format is described below and examples are provided in the \<input\> directory in this package.  We provide an example dataset names 'exampledata', i.e. we provide input files 'exampledata_meas.txt' and 'exampledata_err.txt' are the input files. 

From the commandline simply use as:
```
python mcz.py <filename> nsample --path PATH 
```
<b>\<filename\></b>: \<filename\> is the root name for the input as well as for all the output files: it could be for example the name of the supernova at the location of which HII regions have been measured (\<filename\> is 'exampledata' in our examples)

<b>nsample</b>: the number of MC samples desired (for exampledata 2000 is an appropriate sample size. 
                       nsample must be >100 or it can be 0 to produce a single value, and no confidence region 

additional command line arguments

 <b> --path PATH  </b>         the directory in which subdirectory "input" is located. If not provided, will default
                        environmental variable MCMetdata that should be set to point to that directory. 
                        In our example the input directory is in pyMCZ (level 1).
                        _err.txt _meas.txt must live in \<path\>/input
                        
 <b> --md MD      </b>         metallicity scales to calculate. default is 'all',
                        options are: D02, Z94, M91, M08, P05,P10, PP04, M13, D13, KD02,
                        KD02comb, DP00 (deprecated), P01 (deprecated), C01 (deprecated)
                        
<b>  --unpickle    </b>        if it exists, a pickle file generated in a previous run for this \<filename\> and this 
                        number of samples is read in instead of recalculating the metalicity

<b>  --binmode  BM  </b>         how to choose the number of bins for plotting the histogram:
                            'd' is based on Doane's formula (wikipedia's version),  
                            's' is the sqrt of number of data,        
                            't' is on 2*n**1/3 , 
                            'kb' uses Knuth's block rule (default), 
                            'bb' uses bayesian blocks (must have astroML installed or it defaults to 'kb')
                            'kd' is the kernel density, which requires sklearn installed and additioinally plots the                             histogram with 'kb' mode

<b>  --clobber   </b>          If set to true, will overwrite existing output files without asking. Default False.

<b>  --verbose   </b>          verbose mode. Default False.

<b>  --nodust    </b>          don't do dust corrections (default is to do it)

<b>  --noplot    </b>          don't plot individual distributions (default is to
                        plot all distributions)

<b>  --asciiout   </b>         write distribution medians and 68% inclusion regions in an ascii output (default is not
                        to)
                        
<b>  --asciidistrib  </b>       write the entire distribution for every scale in an ascii file output (default is not to)
                        
                        
<b>  --multiproc (--multi)  </b>         multiprocess, with number of threads nps=max(available cores-1, MAXPROCESSES)

<b>  --log LOGFILE  </b>       outputting messages to a log file instead of standard output. disabled when multiprocessing



Input file format
====================
Flux data for each object should be stored in the directory \<input\> that exists in the directory provided by the --path arg or by the environmental variable MCMetdata. Their ascii file names are expected to end by "\_meas.txt" for the file containing the measurements, and "\_err.txt" for the error file.

with common filename \<filename\>:

\<filename\>_err.txt

\<filename\>_meas.txt 

The format for each of the txt files should be as follows, e.g.,


;# galnum,[OII]3727,Hb,[OIII]4959,[OIII]5007,[OI]6300,Ha,[NII]6584,[SII]6717,[SII]6731,[SIII]9069,[SIII]9532
       1     0.0     0.0     0.0     0.0     0.0   5.117   0.998     0.0     0.0     0.0     0.0
       2     0.0     0.0     0.0     0.0     0.0   5.031   1.012     0.0     0.0     0.0     0.0
       
       
galnum is a sequential index


The first row is optional - it can contain more keys than columns, as long as the columns are listed in the same order as specified above, up to whichever line is of interest. Missing data should be filled in with 'NaN's.

The data should be in the above order, separated by any number of white spaces.


Output
====================
All the results will be saved in the directory "output/\<filename\>" which will be creates in the directory provided by the --path arg.

Read our [paper](http://www.sciencedirect.com/science/article/pii/S2213133716300178), or the README.md file in the output directory in this package https://github.com/nyusngroup/pyMCZ/blob/master/output/README.md for more details on the output products.


Tests implemented
====================

A few tests are implemented to make sure your inputs are valid, and your outputs are robust. 
In order to assess if the number of samples requested is sufficiently large, we provide the module testcompleteness.py. 

Use as, for example: 
	
	
	import testcompleteness as tc

	tc.fitdistrib(\<path to pickle file\>)
	
This reads in the pickle file generated by the code, and compares the cumulative distributions for 4 parameters: E(B-V), D02, Z94 and KD02comb_updated (the user can choose to use whichever scale of course!) (choosing subsets with minimal invalid outpus for each subset).  

The cumulative distribution for 1/10, 1/4, 1/2, 3/4 and the full sample generated by the code are compared with a KS tests and plotted. 

If the number of samples is sufficient, you should expect the probability of the 3/4 sample and full sample to come from the same parent distribution to be close to 1, and generally the KD p-value to increase with the sample fraction (but not always it turns out...). More importantly though this can be used visually: the cumulative distribution should be nearly identical (and smooth). If all distributions overlap, then you are sure that you are well above the critical sample size that achieves smoothness (by a factor 10).

The figures below show an undersampled realization and a well- (possibly) over-sampled realization.




![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/exampledata_n200_testcomplete.png)


![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/exampledata_n2000_testcomplete.png)


![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/exampledata_n20000_testcomplete.png)


 
Known Issues and TODO
====================


Plot formatting is designed for a platform using latex and with Times New Roman serif font available to matplotlib. The module pylabsetup loaded early on assures that the matplotlib rcparameters are set up appropriately, including font choice, and in case of missing fonts an error message is streamed (but not paused upon). If your plots don't look good, change the necessary parameters in pylabsetup.py. This is an unfortunately common occurrence when not running on a Mac.


Note that while a PDF version of our paper is included in the repository, the PDF does not render properly in github. To see the paper please go to the arxiv link (http://arxiv.org/abs/1505.06213) or compile it with pdflatex in the repo "paper" directory on your local machine.

The input file is a but awkward, for historical reasons: we decided to make minimal modification tot he input files required by the original IDL code. An alternative, more pythonic input, may be made available in the future, when we have time.

It would be great to create sphyinx documentation.


Packages required
====================
required packages:

numpy,pylab,matplotlib,scipy

desirable packages:

pickle,multiprocessing,itertools,csv,cProfile,pyqz


Some figures from the paper as examples of code products
===================


![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/pyMCZ_KD02distrib.png)


![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/pyMCZ_KD02distrib_KDE.png)


![alt tag](https://github.com/nyusngroup/pyMCZ/blob/master/figures/pyMCZ_boxplot.png)

