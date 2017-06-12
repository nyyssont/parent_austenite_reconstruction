# parent_austenite_reconstruction

A package for reconstructing the austenite orientation map from martensitic EBSD data with Matlab and MTEX 4.5.0. MTEX should be up and running on Matlab before the scripts can be used. MTEX is available here:
http://mtex-toolbox.github.io/

Scripts pag1, pag2 and pag3 should be run sequentially. The script pag4 will produce a text file that can be imported to Channel 5 analysis software.

There is an example *.cpr datafile in the attached *.zip package that should work with the script.

The script package is meant for the analysis of Oxford Instruments EBSD data containing Iron bcc (old) and Iron fcc data points.

In the script "pag1_import_data.m" you will need to  edit line 32 with the filename of the *.cpr file located in the same folder to be analyzed. The script can then be run, no manual data import necessary.

In the script "pag2_fcc2bcc_detection.m" you can try changing the cutoff value on line 23 to a higher value if you have problems.

In the script "pag3_merge_mcl.m" you can try changing the exponent value p to a higher value to produce smaller clusters and a more rigorous result.

The orientation relationship determination algorithm in script 2 is explained in more detail here:

https://link.springer.com/article/10.1007/s11661-016-3462-2

The algorithm for reconstructing the parent austenite orientation map is explained here (download the document):

http://urn.fi/URN:ISBN:978-952-15-3896-4

If you have problems, please send feedback to tuomoknyyssonen@gmail.com.
