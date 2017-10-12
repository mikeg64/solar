Routines are used to creat VAC/SAC format input files

Creating atmospheric profiles

The data from the valIIIc and McWhirter model is contained in atmos.xls
we also use the hdf inspired data structures

Main routines to use here are
createmodel - uses valiiic data to interpolate between photosphere and 6Mm
createextmodel  - uses valiiic data to interpolate between photosphere and 6Mm. Use data fitting to extrapolate data
                  to corona upto 12Mm

Both of these routines can call generatefield
generatefield  - build magnetic flux tube uses hydrostatic pressure balance to compute correct pressure and energy for the model

testgeneratefield
testsac3dread
testsac3dwrite