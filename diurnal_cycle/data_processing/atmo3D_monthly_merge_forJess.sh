#!/bin/bash

module load cdo

indir="/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/"
outdir="/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/"
tmpdir="/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/tmp/"
codedir="/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/data_processing/"

cdoout=${codedir}cdoout_forJess
touch ${cdoout}

python ${codedir}monthly_hours_since.py 1980-01-02T00:00 2002-01 2005-12 > ${codedir}monthly_hours_since_1980010200.txt

arr=()
while IFS= read -r line; do
   arr+=("$line")
done < ${codedir}monthly_hours_since_1980010200.txt
nstarts=$(expr ${#arr[@]} - 1)

i=0
#while [ $i -lt 1 ]
while [ $i -lt ${nstarts} ]
do
    
    START1="$(date +%s)"

    hs=${arr[i]}
    hn=${arr[$((i+1))]}
    he=$((hn-1)) 
    
    fnlist=$(seq -f pgbf%.0f.gdas2.1980010200 ${hs} ${he} | tr '\n' ' ') 
    fnarray=($fnlist) 
    nhours=${#fnarray[@]} 
    
    if [ ! -f "${tmpdir}${fnarray[0]}_3D.nc" ]; then
        echo "Converting GRIB files to NetCDF..."
        j=0
        while [ $j -lt ${nhours} ]; do
	        echo ${indir}${fnarray[j]}" "${tmpdir}${fnarray[j]}"_air.nc" >> ${cdoout}
	        cdo -f nc -copy -selname,var33,var34,var39,var7,var11,var51 -selltype,100 -sellevel,70000,75000,77500,80000,82500,85000,87500,90000,92500,95000,97500,100000 ${indir}${fnarray[j]} ${tmpdir}${fnarray[j]}"_3D.nc" >> ${cdoout} 2>&1
	        j=$((j+1))
        done
    fi

    echo "Calculating monthly mean diurnal cycle..."
    python save_atmo3D_monthlyDiurnalCycle_forJess.py points ${hs} ${he} ${tmpdir} ${outdir} "${fnlist}"
    pyran=$?

    i=$((i+1))

    END1="$(date +%s)"
    DUR1=$[ ${END1} - ${START1} ]
    echo "Time taken to process monthly data: ${DUR1} seconds"

done

if [ $pyran == 0 ]; then
    echo "Tidying up temporary files..."
    rm ${tmpdir}*"nc"
fi
