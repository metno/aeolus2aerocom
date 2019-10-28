# aeolus2aerocom
aeolus2aerocom is a converter for L2A data of the ESA satellite 
[Aeolus](http://www.esa.int/Our_Activities/Observing_the_Earth/Aeolus) from the ESA binary format to netcdf files used by Met Norway's general evaluation tool [pyaerocom](https://pyaerocom.met.no/). It is intended to be included into pyaerocom at a later stage in the development.

The reading is done using ESA's reading library [coda](http://stcorp.nl/coda/) 
(provided by [Science [&] Technology Corporation](http://stcorp.nl)). coda is now also [packaged for conda](https://anaconda.org/stcorp/coda).

The purpose of aeolus2aerocom is to convert the binary file to a netcdf file that can be read without coda by 
other programs. The contents of the resulting netcdf files are focused on [Met Norway's](https://www.met.no/) cal/val effort
for the Aeolus mission at this point and might therefore not contain all the information of the binary file.

It can also plot the location of satellite overpaths ![satellite overpaths](https://github.com/metno/aeolus2aerocom/blob/master/example_files/20180912_map_colored.png) 

and colocated curtain plots

![curtain plots](https://github.com/metno/aeolus2aerocom/blob/master/example_files/20180912_curtain.png) 

# prerequests
The converter needs a working python interface to [coda](http://stcorp.nl/coda/). You can either install it into your conda environment, or compile it yourself using the information provided at stcorp's [website for coda](http://stcorp.nl/coda/).

You need the coda definition file for the Aeolus L2A product that you can download [here](https://github.com/stcorp/codadef-aeolus/releases/download/20170913/AEOLUS-20170913.codadef) in case you did not use the conda package.
In addition, [coda](http://stcorp.nl/coda/) needs to know where to find this file. 
This can be done by either putting it at the default location 
(/usr/local/share/coda/definitions) or set the environment variable CODA_DEFINITION to the folder where 
the definition file is located. Please look [here](http://stcorp.nl/coda/registration/download/) for details.
