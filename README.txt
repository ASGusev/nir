This program is for getting some data about MS Deconv and Thermo Xtract deconvolution output.

The command line arguements are:
-f to add file name
-l to load file names from a text file
-e to specify good EValue border
-a to specify mass accuracy
-t to add task

The following tasks are possible:
1 Check for scans with zero mass.
2 Check for scans without peaks.
3 Check for scans with peaks found by MS Deconv only.
4 Check for scans that have the same mass found by both programs.
5 Check how much scans by peptide length have correct mass found by each program.
6 Find scans distribution by ratio between Thermo Xtract mass and theoretical mass.
7 Check for scans with certain ratio between mass found by Thermo Xtract and theoretic mass.
