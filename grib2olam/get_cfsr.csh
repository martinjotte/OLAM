#!/bin/csh -f

@ start_year  = 2006
@ start_month = 01
@ start_day   = 01

@ end_year  = 2006
@ end_month = 01
@ end_day   = 02

set times = ( 00 06 12 18 )
# set times = ( 00 12 )
# set times = ( 00 )

# This shell script downloads CFSRv1 (Jan 1979 - April, 2011)
# or CFSRv2 (from April, 20011) grib reanalysis files from NOAA.
# Set the start and end dates and the desired times.
# All files between the start and end times will be downloaded.
#
# This script uses "wget" to download the files.
#
#################################################################

set head1 = "pgbh00.gdas."
set head2 = "cdas1.t"

set tail1 = ".grb2"
set tail2 = "z.pgrbh00.grib2"

set url1 = "https://www.ncei.noaa.gov/data/climate-forecast-system/access/reanalysis/6-hourly-by-pressure-level"
set url2 = "https://www.ncei.noaa.gov/data/climate-forecast-system/access/operational-analysis/6-hourly-by-pressure"

set ndays = ( 31 28 31 30 31 30 31 31 30 31 30 31 )

@ year  = ${start_year}
@ month = ${start_month}
@ day   = ${start_day}

while ( $year <= $end_year )

    set zyr = `printf "%04d" $year`

    if ( $year == $start_year ) then
        @ month = $start_month
    else
        @ month = 01
    endif

    if ( $year == $end_year ) then
        @ emonth = $end_month
    else
        @ emonth = 12
    endif

    while ( $month <= $emonth )

        set zmn = `printf "%02d" $month`

        if ( ($year == $start_year) && ($month == $start_month) ) then
            @ day = $start_day
        else
            @ day = 01
        endif

        if ( ($year == $end_year) && ($month == $end_month) ) then
            @ eday = $end_day
        else
            @ eday = $ndays[$month]
        endif

        while ( $day <= $eday )

            set zdy = `printf "%02d" $day`

            foreach tim ( $times )

                set ztm = `printf "%02d" $tim`

                if ( (($year == 2011) && ($month >= 4)) || ($year > 2011) ) then
                    set file   = ${head2}${ztm}${tail2}
                    set toturl = ${url2}/${zyr}/${zyr}${zmn}/${zyr}${zmn}${zdy}/$file
                    echo wget -4 -nv $toturl -O ${file:r}.${zyr}${zmn}${zdy}.grib2
                    wget -4 -nv $toturl -O ${file:r}.${zyr}${zmn}${zdy}.grib2
                else
                    set file   = ${head1}${zyr}${zmn}${zdy}${ztm}${tail1}
                    set toturl = ${url1}/${zyr}/${zyr}${zmn}/${zyr}${zmn}${zdy}/$file
                    echo wget -4 -nv $toturl
                    wget -4 -nv $toturl
                endif

            end

            @ day++
        end

        @ month++
    end

    @ year++
end
