# aeolus2netcdf
aeolus2netcdf is a converter for L2A data of the ESA satellite 
[ADM Aeolus](http://www.esa.int/Our_Activities/Observing_the_Earth/Aeolus) from the ESA binary format to netcdf files.

The reading is done using ESA's reading library [coda](http://stcorp.nl/coda/) 
(provided by [Science [&] Technology Corporation](http://stcorp.nl)). coda is now also packaged for conda [link](https://anaconda.org/stcorp/coda)

The purpose of aeolus2netcdf is therefore, to just convert the binary file to a netcdf file that can be easier read by 
other programs. The contents of the resulting netcdf files are focused on [Met Norway's](https://www.met.no/) cal / val effort
for the ADM Aeolus mission at this point and might therefore not contain all the information in the binary file.

# prerequests
The converter needs a working python interface to [coda](http://stcorp.nl/coda/). You can either install it if you use anaconda as
your python package manager, or install it yourself using the information provided at stcorp's [website for coda](http://stcorp.nl/coda/).
You also need the coda definition file for the Aeolus L2A product that you can download [here](https://github.com/stcorp/codadef-aeolus/releases/download/20170913/AEOLUS-20170913.codadef).
You also have to tell coda where to find this file by either putting it at the default location 
(/usr/local/share/coda/definitions) or set the environment variable CODA_DEFINITION to the folder where the dinition file is located. Please look at [here](http://stcorp.nl/coda/registration/download/) for details.
