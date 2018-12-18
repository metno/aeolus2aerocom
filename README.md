# aeolus2netcdf
aeolus2netcdf is a converter for L2A data of the ESA satellite 
[ADM Aeolus](http://www.esa.int/Our_Activities/Observing_the_Earth/Aeolus) from the ESA binary format to netcdf files.

The reading is done using ESA's reading library [coda](http://stcorp.nl/coda/) 
(provided by [Science [&] Technology Corporation](http://stcorp.nl)). Because the coda library is not packaged 
in conda or for Ubuntu, it has to be compiled for every update of either coda itsself or the libraries used by it. 

The purpose of aeolus2netcdf is therefore, to just convert the binary file to a netcdf file that can be easier read by 
other programs. The contents of the resulting netcdf files are focused on (Met Norway's)[https://www.met.no/] cal / val effort
for the ADM Aeolus mission at this point and might therefore not contain all the information in the binary file.
