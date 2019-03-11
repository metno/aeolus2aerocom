#!/usr/bin/env bash

retrieval=${1}
if [[ $# -lt 1 ]]
    then echo "usage: ${0} <retrieval>"
    echo "retrieval can be one of sca, mca or ica"
    exit
fi

set -x
basedir='/lustre/storeB/project/fou/kl/admaeolus/'
datadir="${basedir}data.rev.TD01/"
#retrieval='mca'
#downloaddir="${datadir}download/"
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181120T144402034_005448000_001423_0001/"
downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_201*/"
netcdfdir="${datadir}netcdf/"
netcdfdir="${datadir}netcdf_${retrieval}/"
#netcdfdir="${datadir}netcdf_emep_domain/"

jobfile="./TD01.run.${retrieval}.txt"
rm -f "${jobfile}"

mkdir -p "${netcdfdir}"
for file in `find ${downloaddir} -name '*.DBL' | grep 201810 | sort`
    do echo ${file}
#    cmd="./read_aeolus_l2a_data.py --outdir ${netcdfdir} --plotmap --plotprofile --tempdir /tmp/ --file ${file}"
    cmd="./read_aeolus_l2a_data.py --retrieval ${retrieval} --outdir ${netcdfdir} --plotmap --plotprofile --tempdir /tmp/ --file ${file}"
    echo ${cmd}
    echo ${cmd} >> "${jobfile}"
done

# command to start using gnu parallel
/usr/bin/parallel -vk -j 4 -a "${jobfile}"
