#!/usr/bin/env python3
################################################################
# read_aeolus_l2a_data.py
#
# read binary ESA L2A files of the ADM Aeolus mission
#
# this file is part of the pyaerocom package
#
#################################################################
# Created 20190104 by Jan Griesfeller for Met Norway
#
# Last changed: See git log
#################################################################

# Copyright (C) 2019 met.no
# Contact information:
# Norwegian Meteorological Institute
# Box 43 Blindern
# 0313 OSLO
# NORWAY
# E-mail: jan.griesfeller@met.no
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA

"""
object to read colocation netcdf files



"""
import os
import glob
import numpy as np

import logging
import time
import geopy.distance
#import coda



class ReadCoLocationData:

    _FILEMASK = '*.nc'
    __version__ = "0.01"
    DATASET_NAME = 'co-location data'
    DATASET_PATH = '/lustre/storeB/project/fou/kl/admaeolus/data.rev.TD01/download/'
    # Flag if the dataset contains all years or not
    DATASET_IS_YEARLY = False

    FILE_MASK = '*AE_OPER_ALD_U_N_2A*'

    _TIMEINDEX = 0
    _LATINDEX = 1
    _LONINDEX = 2
    _ALTITUDEINDEX = 3
    _EC355INDEX = 4
    _BS355INDEX = 5
    _SRINDEX = 6
    _LODINDEX = 7
    # for distance calculations we need the location in radians
    # so store these for speed in self.data
    # the following indexes indicate the column where that is stored
    _RADLATINDEX = 8
    _RADLONINDEX = 9
    _DISTINDEX = 10

    _COLNO = 11
    _ROWNO = 100000
    _CHUNKSIZE = 10000
    _HEIGHTSTEPNO = 24

    # variable names
    # dimension data
    _LATITUDENAME = 'latitude'
    _LONGITUDENAME = 'longitude'
    _ALTITUDENAME = 'altitude'

    _EC355NAME = 'ec355aer'
    _BS355NAME = 'bs355aer'
    _LODNAME = 'lod'
    _SRNAME = 'sr'
    _CASENAME = 'case'
    _TIME_NAME = 'time'

    # create a dict with the aerocom variable name as key and the index number in the
    # resulting numpy array as value.
    INDEX_DICT = {}
    INDEX_DICT.update({_LATITUDENAME: _LATINDEX})
    INDEX_DICT.update({_LONGITUDENAME: _LONINDEX})
    INDEX_DICT.update({_ALTITUDENAME: _ALTITUDEINDEX})
    INDEX_DICT.update({_TIME_NAME: _TIMEINDEX})
    INDEX_DICT.update({_EC355NAME: _EC355INDEX})
    INDEX_DICT.update({_BS355NAME: _BS355INDEX})
    INDEX_DICT.update({_LODNAME: _LODINDEX})
    INDEX_DICT.update({_SRNAME: _SRINDEX})

    # NaN values are variable specific
    NAN_DICT = {}
    NAN_DICT.update({_LATITUDENAME: -1.E-6})
    NAN_DICT.update({_LONGITUDENAME: -1.E-6})
    NAN_DICT.update({_ALTITUDENAME: -1.})
    NAN_DICT.update({_EC355NAME: -1.E6})
    NAN_DICT.update({_BS355NAME: -1.E6})
    NAN_DICT.update({_LODNAME: -1.})
    NAN_DICT.update({_SRNAME: -1.})

    TEX_UNITS = {}
    TEX_UNITS['ec355aer'] = r'$10^{-6} \cdot m^{-1}$'
    TEX_UNITS['bs355aer'] = ''

    PROVIDES_VARIABLES = list(INDEX_DICT.keys())


    def __init__(self, index_pointer=0, loglevel=logging.INFO, verbose=False):
        self.verbose = verbose
        self.metadata = {}
        self.data = []
        self.index = len(self.metadata)
        self.files = []
        self.index_pointer = index_pointer
        # that's the flag to indicate if the location of a data point in self.data has been
        # stored in rads in self.data already
        # trades RAM for speed
        self.rads_in_array_flag = False

        if loglevel is not None:
            self.logger = logging.getLogger(__name__)
            if self.logger.hasHandlers():
                # Logger is already configured, remove all handlers
                self.logger.handlers = []
            # self.logger = logging.getLogger('pyaerocom')
            default_formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(default_formatter)
            self.logger.addHandler(console_handler)
            self.logger.setLevel(loglevel)
            self.logger.debug('init')

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.metadata[float(self.index)]

    def __str__(self):
        stat_names = []
        for key in self.metadata:
            stat_names.append(self.metadata[key]['station name'])

        return ','.join(stat_names)

    ###################################################################################
    def ndarr2data(self, file_data):
        """small helper routine to put the data read by the read_file method into
        the ndarray of self.data"""

        # start_read = time.perf_counter()
        # return all data points
        num_points = len(file_data)
        if self.index_pointer == 0:
            self.data = file_data
            self._ROWNO = num_points
            self.index_pointer = num_points

        else:
            # append to self.data
            # add another array chunk to self.data
            self.data = np.append(self.data, np.zeros([num_points, self._COLNO], dtype=np.float_),
                                  axis=0)
            self._ROWNO = num_points
            # copy the data
            self.data[self.index_pointer:, :] = file_data
            self.index_pointer = self.index_pointer + num_points

            # end_time = time.perf_counter()
            # elapsed_sec = end_time - start_read
            # temp = 'time for single file read seconds: {:.3f}'.format(elapsed_sec)
            # self.logger.warning(temp)

    ###################################################################################

    def read_file(self, filename, vars_to_read=['ec355aer'], engine='xarray'):
        """method to read one file

        >>> import read_colocation_files
        >>> obj=read_colocation_files.ReadCoLocationData()
        >>> aeolus_file = '/lustre/storeB/project/fou/kl/admaeolus/data.rev.TD01/netcdf_emep_domain/AE_TD01_ALD_U_N_2A_20181130T032226039_005423993_001574_0001.DBL.nc'
        >>> model_file = '/lustre/storeB/project/fou/kl/admaeolus/EMEPmodel.colocated/AE_TD01_ALD_U_N_2A_20181130T032226039_005423993_001574_0001.DBL.colocated.nc'
        >>> aeolus_data = obj.read_file(aeolus_file)
        >>> model_data = obj.read_file(model_file)

        >>> import xarray as xr
        >>> aeolus_data = xr.open_dataset(aeolus_file)
        >>> model_data = xr.open_dataset(model_file)

        """

        vars_to_read_arr = [self._TIME_NAME, self._LATITUDENAME, self._LONGITUDENAME, self._ALTITUDENAME]
        vars_to_read_arr.extend(vars_to_read)

        if engine == 'xarray':
            # read data using xarray
            import xarray as xr
            file_data = xr.open_dataset(filename)
            file_data.close()
            num_points = len(file_data[self._TIME_NAME])
            ret_data = np.empty([num_points, self._COLNO], dtype=np.float_)
            ret_data[:] = np.nan
            for var in vars_to_read_arr:
                if var == self._TIME_NAME:
                    # convert time to datetime64
                    # xarray represents the time as datetime64[ns], but we use datetime64[s]
                    # internally
                    ret_data[:,self.INDEX_DICT[var]] = file_data[var].data.astype('datetime64[s]')
                else:
                    ret_data[:,self.INDEX_DICT[var]] = file_data[var].data

            return ret_data

        elif engine == 'netcdf4':
            # read data using netcdf 4 python
            pass
        else:
            # print error message and return
            return []
            pass

    ###################################################################################

    def plot_profile(self, data_dict, plotfilename, vars_to_plot = ['ec355aer'], title=None,
                     retrieval_name=None,
                     plot_range=(0., 2000.),
                     plot_nbins=20.,
                     linear_time=False):
        """plot sample profile plot

        >>> import read_colocation_files
        >>> obj=read_colocation_files.ReadCoLocationData()
        >>> aeolus_file = '/lustre/storeB/project/fou/kl/admaeolus/data.rev.TD01/netcdf_emep_domain/AE_TD01_ALD_U_N_2A_20181130T032226039_005423993_001574_0001.DBL.nc'
        >>> model_file = '/lustre/storeB/project/fou/kl/admaeolus/EMEPmodel.colocated/AE_TD01_ALD_U_N_2A_20181130T032226039_005423993_001574_0001.DBL.colocated.nc'
        >>> data_dict = {}
        >>> data_dict['aeolus'] = obj.read_file(aeolus_file)
        >>> data_dict['emep'] = obj.read_file(model_file)
        >>> obj.plot_profile(data_dict, './test.png')

        """
        import matplotlib.pyplot as plt
        from scipy import interpolate
        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator


        height_step_no = self._HEIGHTSTEPNO
        target_height_no = 2001
        target_heights = np.arange(0, target_height_no) * 10
        # target_heights = np.flip(target_heights)
        plot_row_no = len(data_dict)
        # enable TeX
        # plt.rc('text', usetex=True)
        # plt.rc('font', family='serif')
        fig, _axs = plt.subplots(nrows=plot_row_no, ncols=1, constrained_layout=True )
        # fig.subplots_adjust(hspace=0.3)
        try:
            axs = _axs.flatten()
        except:
            axs = [_axs]

        plot_handle = []
        levels = []
        cmap = []
        norm = []
        yticks = []
        times = {}
        times_no = {}
        unique_times = {}
        unique_indexes = {}
        unique_height_step_no = {}
        time_step_no = {}

        vars_to_plot_arr = ['altitude']
        vars_to_plot_arr.extend(vars_to_plot)
        for plot_index, data_name in enumerate(data_dict):
            # read returning a ndarray
            data = data_dict[data_name]

            times[data_name] = data[:,self._TIMEINDEX]
            times_no[data_name] = len(times)
            plot_data = {}
            plot_data_masks = {}
            unique_times[data_name], unique_indexes[data_name], unique_height_step_no[data_name] = \
                np.unique(times[data_name], return_index=True, return_counts=True)
            time_step_no[data_name] = len(unique_times[data_name])

            target_x = np.arange(0,time_step_no[data_name])

            for data_var in vars_to_plot_arr:
                # plot_data[data_var] = \
                #     self.data[:, self.INDEX_DICT[data_var]]
                plot_data_masks[data_var] = np.isnan(data[:, self.INDEX_DICT[data_var]])

            # in case of a cut out area, there might not be all the height steps
            # in self.data (since the Aeolus line of sight is tilted 35 degrees)
            # or due to the fact the the slection removes points where longitude or
            # latitude are NaN
            # unfortunately the number of height steps per time code is not necessarily equal
            # to self._HEIGHTSTEPNO anymore
            # e.g. due to an area based selection or due to NaNs in the profile
            # we therefore have to go through the times and look for changes

            # idx_time = times[0]
            # time_cut_start_index = 0
            # time_cut_end_index = 0
            time_index_dict = {}
            for idx, time in enumerate(unique_times[data_name]):
                time_index_dict[time] = np.arange(unique_indexes[data_name][idx],
                                                  unique_indexes[data_name][idx]+unique_height_step_no[data_name][idx])
            #     if time == idx_time:
            #         time_cut_end_index = idx
            #     else:
            #         time_cut_end_index = idx
            #         time_index_dict[idx_time] = np.arange(time_cut_start_index, time_cut_end_index )
            #         time_cut_start_index = idx
            #         idx_time = time
            # time_index_dict[idx_time] = np.arange(time_cut_start_index, time_cut_end_index + 1)
            # time_index_dict[idx_time] = np.arange(unique_indexes[data_name]time_cut_start_index, time_cut_end_index + 1)

            for var in vars_to_plot:
                # this loop has not been optimised for several variables
                out_arr = np.zeros([time_step_no[data_name], target_height_no])
                out_arr[:] = np.nan
                for time_step_idx, unique_time in enumerate(unique_times[data_name]):
                    var_data = data[time_index_dict[unique_time],self.INDEX_DICT[var]]
                    var_data = data_dict[data_name][time_index_dict[unique_time],self.INDEX_DICT[var]]
                    # scipy.interpolate cannot cope with nans in the data
                    # work only on profiles with a nansum > 0

                    nansum = np.nansum(var_data)
                    if nansum > 0:
                        height_data = data_dict[data_name][time_index_dict[unique_time],self.INDEX_DICT['altitude']]
                        if np.isnan(np.sum(var_data)):
                            height_data = height_data[~plot_data_masks[var][time_index_dict[unique_time]]]
                            var_data = var_data[~plot_data_masks[var][time_index_dict[unique_time]]]

                        f = interpolate.interp1d(height_data, var_data, kind='nearest', bounds_error=False, fill_value=np.nan)
                        interpolated = f(target_heights)
                        out_arr[time_step_idx,:] = interpolated
                    elif nansum == 0:
                        # set all heights of the plotted profile to 0 since nothing was detected
                        out_arr[time_step_idx,:] = 0.

                # levels = MaxNLocator(nbins=15).tick_values(np.nanmin(out_arr), np.nanmax(out_arr))
                # levels = MaxNLocator(nbins=20).tick_values(0., 2000.)
                levels = MaxNLocator(nbins=plot_nbins).tick_values(plot_range[0],plot_range[1])
                # cmap = plt.get_cmap('PiYG')
                cmap = plt.get_cmap('autumn_r')
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

                plot_handle.append(axs[plot_index].pcolormesh(out_arr.transpose(), cmap=cmap, norm=norm))
                yticks.append(plot_handle[plot_index].axes.get_yticks())
                yticklabels = plot_handle[plot_index].axes.set_yticklabels(yticks[-1]/100.)
                plot_handle[plot_index].axes.set_xlabel('time step number')
                plot_handle[plot_index].axes.set_ylabel('height [km]')
                if title:
                    plot_handle[plot_index].axes.set_title(title, fontsize='small')
                else:
                    plot_handle[plot_index].axes.set_title(data_name)
                #plot_simple2.axes.set_aspect(0.05)
                # plt.show()
                plt.annotate('start={}'.format(unique_times[data_name][0].astype('datetime64[s]')),
                             (-0.12, -0.24),
                             xycoords=plot_handle[plot_index].axes.transAxes,
                             fontsize=8,
                             horizontalalignment='left')
                plt.annotate('end={}'.format(unique_times[data_name][-1].astype('datetime64[s]')),
                             (1.02, -0.24),
                             xycoords=plot_handle[plot_index].axes.transAxes,
                             fontsize=8,
                             horizontalalignment='right')

        clb = plt.colorbar(plot_handle[0], ax=axs, orientation='vertical', fraction=0.05,
                           aspect=30)
        if retrieval_name:
            clb.set_label('{} [{}] {} retrieval'.format(var, self.TEX_UNITS[var], retrieval_name),
                             )
        else:
            clb.ax.set_title('{} [{}]'.format(var, self.TEX_UNITS[var]), fontsize='small', orientation='vertical')
        plt.savefig(plotfilename, dpi=300)
        plt.close()
            # print('test')

    ###################################################################################

    ###################################################################################


if __name__ == '__main__':
    import logging

    default_plot_file = './out.png'
    default_topo_file = '/lustre/storeB/project/fou/kl/admaeolus/EMEP.topo/MACC14_topo_v1.nc'
    default_aeolus_dir = '/lustre/storeB/project/fou/kl/admaeolus/data.rev.TD01/netcdf_emep_domain_mca'
    default_model_dir = '/lustre/storeB/project/fou/kl/admaeolus/EMEPmodel.colocated'

    import argparse
    options = {}
    parser = argparse.ArgumentParser(
        description='command line interface to aeolus2netcdf.py\n\n\n')
    parser.add_argument("--aeolusfile", help="aeolus file to read")
    parser.add_argument("--modelfile", help="model file to read")
    parser.add_argument("-v", "--verbose", help="switch on verbosity",
                        action='store_true')
    parser.add_argument("--starttime",
                        help="startdate as time string understood by numpy's datetime64; defaults to None; example: 2018-12-11T00:00:00",
                        default=None)
    parser.add_argument("--endtime",
                        help="enddate as time string; defaults to None.",
                        default=None)
    # parser.add_argument("--outdir", help="output directory; the filename will be extended with the string '.nc'")
    parser.add_argument("--logfile", help="logfile; defaults to /home/jang/tmp/aeolus2netcdf.log",
                        default="/home/jang/tmp/plot_colocation_files.log")

    # parser.add_argument("-O", "--overwrite", help="overwrite output file", action='store_true')
    # parser.add_argument("--emep", help="flag to limit the read data to the cal/val model domain", action='store_true')
    # parser.add_argument("--latmin", help="min latitude to return", default=np.float_(30.))
    # parser.add_argument("--latmax", help="max latitude to return", default=np.float_(76.))
    # parser.add_argument("--lonmin", help="min longitude to return", default=np.float_(-30.))
    # parser.add_argument("--lonmax", help="max longitude to return", default=np.float_(45.))
    parser.add_argument("--aeolusdir",
                        help="base directory of the aeolus files; defaults to {}/".format(default_aeolus_dir),
                        default=default_aeolus_dir)
    parser.add_argument("--modeldir", help="base directory of the model files; defaults to {}/".format(default_model_dir),
                        default=default_model_dir)
    # parser.add_argument("--filemask", help="file mask to find data files",
    #                     default='*AE_OPER_ALD_U_N_2A_*')
    # parser.add_argument("--tempdir", help="directory for temporary files",
    #                     default=os.path.join(os.environ['HOME'], 'tmp'))
    # parser.add_argument("--plotmap", help="flag to plot a map of the data points; files will be put in outdir",
    #                     action='store_true')
    # parser.add_argument("--plotprofile", help="flag to plot the profiles; files will be put in outdir",
    #                     action='store_true')
    parser.add_argument("--plotfile", help="plot file name (abs. path); defaults to {}.".format(default_plot_file),
                        action='store_true',default=default_plot_file)
    parser.add_argument("--variables", help="comma separated list of variables to write; default: ec355aer,bs355aer",
                        default='ec355aer,bs355aer')

    args = parser.parse_args()

    if args.plotfile:
        options['plotfile'] = args.plotfile

    if args.aeolusfile:
        options['aeolusfile'] = args.aeolusfile
        
    if args.modelfile:
        options['modelfile'] = args.modelfile
        
    if args.logfile:
        options['logfile'] = args.logfile
        logging.basicConfig(filename=options['logfile'], level=logging.INFO)

    if args.starttime:
        options['starttime'] = args.starttime

    if args.endtime:
        options['endtime'] = args.endtime

    # if args.outdir:
    #     options['outdir'] = args.outdir
    #
    #
    # if args.tempdir:
    #     options['tempdir'] = args.tempdir
    #
    # if args.emep:
    #     options['emepflag'] = args.emep
    #     options['latmin'] = np.float(30.)
    #     options['latmax'] = np.float(76.)
    #     options['lonmin'] = np.float(-30.)
    #     options['lonmax'] = np.float(45.)
    # else:
    #     options['emepflag'] = False
    #
    # if args.latmin:
    #     options['latmin'] = np.float_(args.latmin)
    #
    # if args.latmax:
    #     options['latmax'] = np.float_(args.latmax)
    #
    # if args.lonmin:
    #     options['lonmin'] = np.float_(args.lonmin)
    #
    # if args.lonmax:
    #     options['lonmax'] = np.float_(args.lonmax)
    #
    # if args.readpaths:
    #     options['readpaths'] = args.readpaths.split(',')

    if args.variables:
        options['variables'] = args.variables.split(',')

    if args.verbose:
        options['verbose'] = True
    else:
        options['verbose'] = False

    # import read_data_fieldaeolus_l2a_data
    import os
    import sys
    import glob
    import pathlib
    import xarray as xr



    aeolus_name = 'aeolus TD01'
    model_name = 'EMEP_MSC-W'
    obj = ReadCoLocationData()
    data_dict = {}
    obj.logger.info('reading aeolus file: {}'.format(options['aeolusfile']))

    data_dict[aeolus_name] = obj.read_file(options['aeolusfile'])
    obj.logger.info('reading model file: {}'.format(options['modelfile']))
    data_dict[model_name] = obj.read_file(options['modelfile'])
    # adjust extinction value
    data_dict[model_name][:,obj.INDEX_DICT['ec355aer']] = data_dict[model_name][:,obj.INDEX_DICT['ec355aer']] *1E6

    obj.logger.info('plotted file: {}'.format(options['plotfile']))
    obj.plot_profile(data_dict, options['plotfile'],retrieval_name='mca',
                                    plot_range=(0.,200.))

