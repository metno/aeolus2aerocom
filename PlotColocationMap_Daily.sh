#!/usr/bin/env bash

retrieval=${1}
if [[ $# -lt 1 ]]
    then echo "usage: ${0} <retrieval>"
    echo "retrieval can be one of sca, mca or ica"
    exit
fi

set -x
basedir='/lustre/storeB/project/fou/kl/admaeolus/'
datadir="${basedir}EMEPmodel.colocated.${retrieval}/"
plotdir="${datadir}plots/"
#downloaddir="${datadir}download/"
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181120T144402034_005448000_001423_0001/"
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181[1-2]*/"
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181*/"
#netcdfdir="${datadir}netcdf/"
#netcdfdir="${datadir}netcdf_emep_domain_${retrieval}/"
#modeloutdir="${basedir}EMEPmodel.colocated.${retrieval}/"

dates=`find ${datadir} -name '*.nc' | cut -d_ -f7 | cut -dT -f1 | uniq`

jobfile="./TD01.run.${retrieval}.emep.plot.daily.txt"
rm -f "${jobfile}"

mkdir -p "${plotdir}"
set -x
for date in `find ${datadir} -name '*.nc' | cut -d_ -f7 | cut -dT -f1 | sort | uniq`
    do echo ${date}
    startdate=`date '+%C%y-%m-%d' -d ${date}`
    enddate=`date '+%C%y-%m-%d' "-d ${date} +1 day"`
    namedate=`echo ${startdate} | sed -e 's/\-//g'`
    plotfile="${plotdir}${namedate}_profile.png"
    mapplotfile="${plotdir}${namedate}_map.png"
    cmd="./read_colocation_files.py --retrieval ${retrieval} --starttime ${startdate} --endtime ${enddate} --plotfile ${plotfile} --plotmapfile ${mapplotfile}"
#    cmd="./read_colocation_files.py --retrieval ${retrieval} --starttime ${startdate} --endtime ${enddate} --plotmapfile ${mapplotfile}"
    echo ${cmd}
    echo ${cmd} >> "${jobfile}"
done

# start using gnu parallel
/usr/bin/parallel -vk -j 5 -a "${jobfile}"
