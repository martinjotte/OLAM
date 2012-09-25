!==========================================================================================!
!==========================================================================================!
!     This is the main driver for file output in ED.                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
                    ,annual_time,writing_dail,writing_mont,writing_dcyc,history_time       &
                    ,dcycle_time,the_end)


   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in)  :: the_end
   logical, intent(in)  :: analysis_time
   logical, intent(in)  :: dail_analy_time
   logical, intent(in)  :: mont_analy_time
   logical, intent(in)  :: dcyc_analy_time
   logical, intent(in)  :: writing_dail
   logical, intent(in)  :: writing_mont
   logical, intent(in)  :: writing_dcyc
   logical, intent(in)  :: history_time
   logical, intent(in)  :: dcycle_time
   logical, intent(in)  :: new_day
   logical, intent(in)  :: annual_time


   return
end subroutine ed_output
