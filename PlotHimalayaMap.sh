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
#downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_20181[1-2]*/"
downloaddir="${datadir}download/AE_TD01_ALD_U_N_2A_201809*/"
#netcdfdir="${datadir}netcdf/"
netcdfdir="${datadir}netcdf_himalaya_domain_${retrieval}/"
#modeloutdir="${basedir}EMEPmodel.colocated.${retrieval}/"

jobfile="./TD01.run.${retrieval}.himalaya.colocation.txt"
rm -f "${jobfile}"

mkdir -p "${netcdfdir}"
for file in `find ${downloaddir} -name '*.DBL' | sort`
    do echo ${file}
    # read_aeolus_l2a_data.py --himalaya --plotmap --outdir ./ --file /lustre/storeB/project/fou/kl/admaeolus/data.rev.TD01/download/AE_TD01_ALD_U_N_2A_20180908T120926033_005387992_000264_0002/AE_TD01_ALD_U_N_2A_20180908T120926033_005387992_000264_0002.DBL
#    cmd="./read_aeolus_l2a_data.py --himalaya --retrieval ${retrieval}  --outdir ${netcdfdir} --plotmap --plotprofile --tempdir /tmp/ --file ${file}"
    cmd="./read_aeolus_l2a_data.py --himalaya --retrieval ${retrieval}  --outdir ${netcdfdir} --plotmap --tempdir /tmp/ --file ${file}"
    echo ${cmd}
    echo ${cmd} >> "${jobfile}"
done

# start using gnu parallel
/usr/bin/parallel -vk -j 5 -a "${jobfile}"
