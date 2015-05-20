MC Metalicity Output
====================

Monte Carlo sampling to calculate metalicity uncertainty from strong line flux data.

====================
Output
====================
All the results will be saved in this direcoty by the name of the input file (e.g. exampledata for input exampledata_meas.txt and exampledata_err.txt), which will be created upon running the code in the directory provided by the --path argument.

For a given file <fname> and nsample, the following output files are generated:

\<fname\>\_n\<nsample\>.pkl": (e.g. exampledata_n2000.pkl) the metallicities calculated and its uncertainties stored in binary (pickled) format in a python dictionary

\<fname\>\_boxplot\_n\<nsample\>\_\<measurement\>.pdf": (e.g. exampledata\_boxplot\_n2000\_1.pdf, exampledata_boxplot2000_m2.pdf) box and whiskers plot showing all metallicities and their uncertainties in a single plot, one plot (and one PDF file) for each measurement in input 

"hist" folder, containing all the histograms generated. For each scale calculated, scale <scale_name>, and each input set of lines, or measurement, the distribution histogram is saved in a PDF file \<fname\>\_\<nsample\>\_\<scale\>\_\<measurement\>.pdf (e.g. exampledata\_n2000\_KD02\_N2O2\_1.pdf for scale KD02_N2O2, measurement 1)

if the keyword --asciiout is given a text file \<fname\>\_n\<nsample\>\_\<measurement\>.txt (e.g. exampledata\_n2000\_1.txt) stores the distribution 50th, 16th and 84th percentiles in the format: 
E(B-V)	 0.121000	 0.014000	 0.019000
logR23	 0.428000	 0.003000	 0.004000
...

if the keyword --asciidistrib is given a csv file \<fname\>\_n\<nsample\>\_\<scale\>\_\<measurement\>.txt (e.g. exampledata\_n2000\_KD02\_N2O2\_1.csv) is created for each scale and each measurement, containing the full distribution of metallicities. 

if the keyword --binmode is set to 'kd' (kernel density) a the kernel density is saved in a binary (pickle) file \<fname\>\_n\<nsample\>\_\<scale\>\_KDE\_\<measurement\>.pkl (e.g. exampledata\_n2000\_KD02\_N2O2\_1\_KDE.pkl) is created for each scale and each measurement. The kernel density is saved as a python sklearn KernelDensity object)



====================
Interpret the Output Values
====================
 
'-1' means: the emissions lines don't satisfy some sort of minimum criteria to derive metallicity  

0.000 (no distribution) means no distribution was generated for some reason: you may have selected nsample=0, or have errors = 0.0, and that is the appropriate behavor in this case. otherwise the distribution may be really screwed up and take a single value, or a shape such that the 16th and 84th percentiles cannot be calculated (e.g. two values), and you should definitely don't trust the output  
