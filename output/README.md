MC Metalicity Output
====================

Monte Carlo method to calculate metalicity uncertainty from flux data.

====================
Output
====================
All the results will be saved in this direcoty or that  which will be created in the directory provided by the --path arg.

For a given \<fi\> and nSample, the following output files are generated:

"\<fi\>_n"nSample"_X.csv": the metallicity and its uncertainty that was calculated, where X will correspond to the number of rows the input data contains

"\<fi\>_n"nSample"_X.pkl": the metallicity and its uncertainty that was calculated stored in binary (pickled) format in a python dictionary

"hist" folder, containing all the histograms generated


 
