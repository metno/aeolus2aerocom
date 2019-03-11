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
#downloaddir="${datadir}download/"
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181120T144402034_005448000_001423_0001/"
downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181[1-2]*/"
downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_201811*/"
#netcdfdir="${datadir}netcdf/"
netcdfdir="${datadir}netcdf_emep_domain_${retrieval}/"
modeloutdir="${basedir}EMEPmodel.colocated.${retrieval}/"

jobfile="./TD01.run.${retrieval}.emep.colocation.txt"
rm -f "${jobfile}"

mkdir -p "${netcdfdir}"
for file in `find ${downloaddir} -name '*.DBL' | sort`
    do echo ${file}
    cmd="./read_aeolus_l2a_data.py --emep --retrieval ${retrieval} --netcdfcolocate --modeloutdir ${modeloutdir} --outdir ${netcdfdir} --plotmap --plotprofile --tempdir /tmp/ --file ${file}"
    echo ${cmd}
    echo ${cmd} >> "${jobfile}"
done

# start using gnu parallel
/usr/bin/parallel -vk -j 5 -a "${jobfile}"
