program clouds_lookup

  use hdf5_utils
  implicit none

  integer, parameter :: n_rh =  56
  integer, parameter :: n_qc = 100
  
  real, parameter :: rh_start = 0.30
  real, parameter :: rh_delta = 0.02
  real, parameter :: rh_end   = rh_start + (n_rh-1) * rh_delta

  real, parameter :: qc_start = 0.025
  real, parameter :: qc_delta = 0.025
  real, parameter :: qc_end   = qc_start + (n_qc-1) * qc_delta

  real :: cldf(n_qc,n_rh)

  real :: rh(n_qc)   ! relative humidity
  real :: qsub(n_qc) ! q_cloud / q_air

  integer :: cnvg(n_qc)
  integer :: irh, iqc
  integer :: ndims, idims(3)

  ! The table is made linear in sqrt(qsub) to better show smaller values
  
  do iqc = 1, n_qc
     qsub(iqc) = ( qc_start + real(iqc - 1) * qc_delta ) ** 2
  enddo

  do irh = 1, n_rh
     rh(:) = rh_start + real(irh - 1) * rh_delta
     call clouds_gno_40(n_qc, rh, qsub, cldf(:,irh), cnvg)
  enddo

  call shdf5_open("bony_eman_table.h5", 'W', 1)
  
  ndims=1 ; idims(1)=1
  call shdf5_orec(ndims, idims, "rh_start", rvars=rh_start)
  call shdf5_orec(ndims, idims, "rh_end"  , rvars=rh_end)
  call shdf5_orec(ndims, idims, "rh_size" , ivars=n_rh)

  call shdf5_orec(ndims, idims, "qc_start", rvars=qc_start)
  call shdf5_orec(ndims, idims, "qc_end"  , rvars=qc_end)
  call shdf5_orec(ndims, idims, "qc_size" , ivars=n_qc)

  ndims=2 ; idims(1:2) = (/ n_qc, n_rh /)
  call shdf5_orec(ndims, idims, "cldfrac" , rvara=cldf)

  call shdf5_close()

end program clouds_lookup
  
  

subroutine clouds_gno_40(nd, relhum, qsub_ov_q, cldf, cnvg)

  implicit none

!==============================================================================
!
!				CLOUDS_GNO 	version 4.0
!
! Purpose:
! --------
!
!   Parameterization of the cloudiness (cloud amount, cloud water content)
!   associated with cumulus convection.
!
! Principle:
! ----------
!
!   This cloud parameterization predicts the cloudiness that is associated with 
!   the presence of condensation within a large-scale domain: this condensation 
!   may be produced at the subgrid-scale by cumulus convection and at the 
!   large-scale by super-saturation.
!
!   IMPORTANT: in the present version of the scheme, the only source of 
!   subgrid-scale condensation that is considered is cumulus convection 
!   (condensation associated with boundary layer turbulence, for instance, 
!   is not considered). 
!
!   The cloud fraction and the in-cloud water content are predicted by a
!   statistical approach. The subgrid-scale variability of total water 
!   (vapor + condensed) within the gridbox is described by a generalized 
!   log-normal Probability Distribution Function (PDF) whose mean, variance 
!   and skewness coefficient are predicted. The predictors are: 
!     1) the local concentration of condensed water that is produced at
!        the subgrid-scale by convection (output of the convection scheme) 
!     2) the saturation deficit or excess of the environment 
!     3) the domain-averaged mixing ratio of total water
!   Note that we impose the distribution of total water to be bounded by zero. 
!   On the other hand, no upper bound of the distribution is considered in this 
!   version of the scheme.
!
!   If no subgrid-scale condensation occurs within the domain, the scheme
!   becomes equivalent to an "all-or-nothing" large-scale saturation scheme.
!
! Inputs:
! -------
!
!  ND----------: Number of vertical levels
!  R--------ND-: Domain-averaged mixing ratio of total water 
!  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
!  QSUB-----ND-: Mixing ratio of condensed water within clouds associated
!                with SUBGRID-SCALE condensation processes (here, it is
!                predicted by the convection scheme)
! Outputs:
! --------
!
!  CLDF-----ND-: cloud fractional area (0-1)
!  CLDQ-----ND-: in-cloud mixing ratio of condensed water (kg/kg)
!
! CALL command:
! -------------
!
!     CALL CLOUDS_GNO(ND,R,RS,QSUBGRID,CLDF,CLDQ)
!
! Reference:
! ----------
!
!   Bony, S and K A Emanuel, 2001: A parameterization of the cloudiness
!       associated with cumulus convection; Evaluation using TOGA COARE data.
!	    J. Atmos. Sci., accepted.
!
! Written by:
! -----------
!
!  Sandrine Bony (MIT & LMD/CNRS; bony@wind.mit.edu) -  July 2000
!
!  Difference with version 1.0:
!	numerical method of resolution of equation 9
!	version 1.0: use a Gaussian PDF when erf(v)->1
!       version 2.0: use an asymptotic expression of erf(v) instead of a 
!       Gaussian PDF
!==============================================================================
!
! -- input/output arguments of the subroutine:

  integer, intent(in)  :: nd
  integer, intent(out) :: cnvg(nd)
  real,    intent(in)  :: relhum(nd), qsub_ov_q(nd)
  real,    intent(out) :: cldf(nd)

! -- parameters controlling the iteration:
! --    nmax    : maximum nb of iterations (hopefully never reached!)
! --    epsilon : accuracy of the numerical resolution (here 2.0%)
! --    vmax    : v-value above which we use an asymptotic expression for ERF(v)

  integer, PARAMETER :: nmax = 200
  REAL,    parameter :: epsilon = 0.001
  REAL,    parameter :: vmax0   = 2.0
  REAL,    parameter :: vmax0sq = vmax0 * vmax0
  real               :: vmax

! -- misc:

  INTEGER :: K, n, m, niter
  REAL    :: mu, qsat, delta, beta 
  REAL    :: xx, aux, coeff, blokk, dist, fprime, det
  REAL    :: erfu, erfv
  real    :: u, v

  real, parameter :: pi = acos(-1.)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! -- loop over vertical levels :

  cnvg(:) = 0

  DO K = 1, ND

!     mu = R(K)
!     mu = MAX(mu,min_mu)
!     qsat = RS(K) 
!     qsat = MAX(qsat,min_mu)
!     delta = log(mu/qsat)
     
     delta = log(relhum(k))

!     IF ( QSUB(K) .lt. min_Q ) THEN

!===========================================================================
!  If no condensation is produced at the subgrid-scale:
!
!  -> the scheme becomes equivalent to a "large-scale condensation scheme"
!     ie: cldf = H(mu-qsat) and cldq = (mu-qsat)*H(mu-qsat)
!     where H is the Heaviside function.
!     (in the absence of subgrid-scale condensation, the generalized
!      log-normal PDF becomes equivalent to a gaussian PDF of variance
!      zero, i.e. it becomes equivalent to a Dirac function and the 
!      cumulative distribution function becomes an Heaviside function).
!===========================================================================

!       CLDQ(K) = MAX( 0.d0, mu-qsat )
!       CLDF(K) = CLDQ(K) / MAX( CLDQ(K), real(min_mu) )
!        cldf(k) = 1.0
!        lambda = mu
!        alpha  = 0.0
!        kew    = 0.0   
!        skew   = 0.0   
!        sigs   = 0.0   

!     ELSE 

!===========================================================================
!  Some condensation is produced at the subgrid-scale: 
!  (presence of subgrid-scale variability):
!
! Use the (iterative) numerical method of Newton to determine the parameters 
! that characterize the PDF of total water.
!
! Remark 1: the accuracy of the resolution is controlled by "epsilon"
! Remark 2: in GCMs, this numerical method may be too much CPU-time consuming.
! In that case, it may be more appropriate to substitute it by a tabulation 
! of equations 9 and 11 (see the Bony-Emanuel article cited in introduction).
!===========================================================================

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        PDF = generalized log-normal distribution (GNO):
!   (k<0 if a lower bound is considered for the PDF of total water)
! 
!   -> determine x (the parameter k of the GNO PDF) 
!      such that the contribution of subgrid-scale processes to the 
!      in-cloud water content is equal to QSUB(K)
!
! NB: the "error function" is called ERF or DERF (in double precision)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!       beta = QSUB(K)/mu + EXP( -MIN(0.d0,delta) )
!!        beta = qsub_ov_q(k) + EXP( -MIN(0.,delta) )
!!
     if (delta >= 0.0) then
        beta = qsub_ov_q(k) + 1.0
     else
        beta = qsub_ov_q(k) + EXP( -delta )
     endif

     vmax = vmax0

     if (delta < 0.0) then
        xx = -0.5*SQRT(log(2.))
     else
        xx = -SQRT(2.)*vmax0*(1.-SQRT(1.+delta/vmax0sq))
        if ( xx .GE. 0. ) xx = -0.5*SQRT(2.)*vmax0*(1.+SQRT(1.+delta/vmax0sq))
     endif

     iter: do n = 1, nmax ! iteration loop

        u = delta/(xx*sqrt(2.)) + xx/(2.*sqrt(2.))
        v = delta/(xx*sqrt(2.)) - xx/(2.*sqrt(2.))

        IF ( v .GT. vmax ) THEN 
           IF ( ABS(u) .GT. vmax .AND. delta .LT. 0. ) THEN
! -- use asymptotic expression of erf for u and v large:
! ( -> analytic solution for xx )
              aux = 2.*delta*(1.-beta*EXP(delta)) &
                    /(1.+beta*EXP(delta))
              xx = -SQRT(aux)
              blokk = EXP(-v*v) / v / sqrt(pi)
              dist = 0.
              fprime = 1.
           ELSE
! -- erfv -> 1.0, use an asymptotic expression of erfv for v large:
              erfu = 1. - erf(u)
              aux = sqrt(pi) * erfu * EXP(v*v)
              coeff = 1. - 0.5/(v**2) + 3./4./(v**4)
              blokk = coeff * EXP(-v*v) / v / sqrt(pi)
              dist = v * aux / coeff - beta
              fprime = 2. / xx * (v**2) &
                     * ( coeff*EXP(-delta) - u * aux ) &
                     / coeff / coeff         
           ENDIF ! ABS(u)
        ELSE
! -- general case:
           erfu = 1. - ERF(u)
           erfv = 1. - ERF(v)
           blokk = erfv
           dist = erfu / erfv - beta
           fprime = 2. / sqrt(pi) / xx / erfv**2 &
                  * (   erfv*v*EXP(-u*u) &
                      - erfu*u*EXP(-v*v) )
        ENDIF

! -- numerical convergence reached?
        if ( ABS(dist/beta) .LT. epsilon ) then
           cnvg(k) = 1
           exit iter
        endif


        xx = xx - dist/fprime

     ENDDO iter

     cldf(k) = 0.5 * blokk

  ENDDO
END subroutine clouds_gno_40
