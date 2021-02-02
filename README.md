# MIRT
A general program for item-response analysis is described which uses the stabilized Newton-Raphson algorithm. This program is written to be compliant with Fortran 2003 standards and is sufficiently general to handle independent variables, multidimensional ability parameters, and matrix sampling. The ability variables may be either polytomous or multivariate normal. Items may be dichotomous or polytomous.

Please see irtprogram.pdf for details of the program. An earlier version of this document was published as Haberman (2013, ETS Research Report RR-13-32). 

## Running the program

Binaries can be found here: https://github.com/EducationalTestingService/MIRT/releases/tag/v1.0.0. The latest release is v1.0.0. 

Input for the program consists of a data file which contains the observations and a control file which follows Fortran 2003 rules for namelist input. For details on how to create contol files, see also the MIRT wiki page: https://github.com/EducationalTestingService/MIRT/wiki. The program is run from a command line. Within a Unix/Linux environment (including MacOS), the command line is opened by opening a terminal. In Windows, command prompt is used to obtain a command line. The user is expected to be able to change directories and perform other basic tasks associated with the basic commands of an operating system. Due to weaknesses in memory management in Windows, program performance is typically much better in a Linux/Unix environment. The program name is mirt, and piping is normally employed on the command line to specify a control file. If control.txt is the name of the control file and if the executable file is in the path used to find commands, then the program is invoked with the following command:

mirt < control.txt

The control file uses namelist input records to specify the data, the model parameters, the data files, other input files, and output files. Output generally is designed to produce files with comma-separated values (csv files) readily treated by standard software for
spreadsheets. 

It is also possible to invoke the program from other programs. Many professional text editors (e.g., Emacs (ESS), UltraEdit, TextPad) can be useful for this purpose, because input files can be edited and run directly from within the editor (e.g., using short cut keys).
If the control file is in the working directory, it is generally possible to invoke the program from R with the following command:

shell("mirt < control.txt")

This approach can be particularly useful for jackknife and bootstrap procedures and simulations. R also has many convenient graphical output options. Furthermore, invocation of the program can also be done via other statistical software (e.g., SAS).

## Compiling from source

The source code can be compiled using gfortran. On Windows, Mingw-w64 also needs to be downloaded and installed. On Windows, the code can be compiled with the following command:

gfortran *.f95 -O3 -ffree-line-length-none -o mirt.exe

Fresh compilations need to be done twice due to module files.
