MC Metalicity
====================

Monte Carlo method to calculate metalicity uncertainty from flux data.

====================
Usage:
====================
Place the _min, _max, and _med(optional) files in the sn_data folder.


From the commandline simply use as:
```
MCZ_err.py <filename> nsample
```
-\<filename\>: the common filename of the _min and _max files

-nsample: the number of iterations desired


====================
Input file format
====================
each data should be in the following format, with common filename \<f\i>_max.txt, _min.txt, or _med.txt where max = med+err, min=med-err .

'''
;# galnum,[OII]3727,Hb,[OIII]4959,[OIII]5007,[OI]6300,Ha,[NII]6584,[SII]6717,[SII]6731,[SIII]9069,[SIII]9532,E(B-V),dE(B-V)
       1     0.0     0.0     0.0     0.0     0.0   5.117   0.998     0.0     0.0     0.0     0.0
       2     0.0     0.0     0.0     0.0     0.0   5.031   1.012     0.0     0.0     0.0     0.0
'''
(example)

galnum is just the index used for matching with the radius for calculating the gradient with.

The first row is optional - any row with more strings than numbers will be ignored

the data should be in the above order, separated by any number of white spaces.

