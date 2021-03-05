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

## Installing MIRT on Windows

Since the MIRT binary for Windows is created with gfortran and Mingw-w64, the latter also needs to be installed. The following steps are needed:
1. Make sure you have administrator rights on your computer.
2. Copy the MIRT binary to a directory (e.g., C:\Program Files\MIRT).
3. Download mingw-w64-install.exe from http://mingw-w64.org/.
4. Run the installer and make sure to specify the correct architecture. This will be "i686" for a 32-bit system, or "x86_64" for a 64-bit system.
5. Once installation is completed, you need to set the path for Mingw-w64: Right click on the "My Computer" or "This PC" icon on your desktop, select "Properties", click "Advanced system settings", and click "Environment variables". In the "Environment Variables" window, select "Path" and then click on "Edit". Add the directory of the relevant Mingw‐w64 libraries. If you installed the 32-bit version, these will be in: C:\Program Files (x86)\mingw‐w64\i686‐8.1.0‐posix‐dwarf‐rt_v6‐rev0\mingw32\binb. If you installed the 64-bit version, these will be in: C:\Program Files\mingw‐w64\x86_64‐8.1.0‐posix‐seh‐rt_v6‐rev0\mingw64\bin.
6. Set the path for MIRT: Add the directory of the MIRT binary in the same way as above.
7. Verify that the path variable has been set correctly by going to the command prompt (Search Windows > cmd > Command Prompt), and type "path". Check that the directories appear correctly.
8. You are now ready to run MIRT from any directory. 






