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

# This shell script downloads CFSRv1 (before April, 2011) or
# CFSRv2 (from April, 20011) grib analysis files from NOAA. 
# Just set the start and end dates and the desired times.
# All files between the start and end times will be downloaded.
# Analyses are available from Jan 01, 1979 till about a week
# before the current date.
#
# This script uses "wget" to download the files. 
#
#################################################################

set head1 = "pgbh00.gdas."
set head2 = "cdas1.t"

set tail1 = ".grb2"
set tail2 = "z.pgrbh00.grib2"

set url1 = "https://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh"
set url2 = "https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh"

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

                if ( ($year >= 2011) && ($month >= 4) ) then
                    set file   = ${head2}${ztm}${tail2}
                    set toturl = ${url2}/${zyr}/${zyr}${zmn}/${zyr}${zmn}${zdy}/$file
                else
                    set file   = ${head1}${zyr}${zmn}${zdy}${ztm}${tail1}
                    set toturl = ${url1}/${zyr}/${zyr}${zmn}/${zyr}${zmn}${zdy}/$file
                endif

                wget -4 -nv $toturl

            end

            @ day++
        end

        @ month++
    end

    @ year++
end
