#!/bin/bash

module load cdo

indir="/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/CFSm501hr1980010200/DATA/"
outdir="/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/"
tmpdir="/lustre/f2/dev/ncep/Jack.Reeveseyre/diurnal_cycle/tmp/"
codedir="/ncrc/home2/Jack.Reeveseyre/cfs-analysis/diurnal_cycle/data_processing/"

cdoout=${codedir}cdoout
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
    
    fnlist=$(seq -f flxf%.0f.gdas2.1980010200 ${hs} ${he} | tr '\n' ' ') 
    fnarray=($fnlist) 
    nhours=${#fnarray[@]} 
    
    if [ ! -f "${tmpdir}${fnarray[0]}.nc" ]; then
        echo "Converting GRIB files to NetCDF..."
        j=0
        while [ $j -lt ${nhours} ]; do
	        echo ${indir}${fnarray[j]}" "${tmpdir}${fnarray[j]}".nc" >> ${cdoout}
	        cdo -f nc -copy -selname,var33,var34 -selltype,105 ${indir}${fnarray[j]} ${tmpdir}${fnarray[j]}".nc" >> ${cdoout} 2>&1
	        j=$((j+1))
        done
    fi

    echo "Calculating monthly mean diurnal cycle..."
    python save_atmo_monthlyDiurnalCycle.py grid ${hs} ${he} ${tmpdir} ${outdir} "${fnlist}"
    pyran=$?

    if [ $pyran == 0 ]; then
	echo "Tidying up temporary files..."
        rm ${tmpdir}*"nc"
    fi

    i=$((i+1))

    END1="$(date +%s)"
    DUR1=$[ ${END1} - ${START1} ]
    echo "Time taken to process monthly data: ${DUR1} seconds"

done
