subroutine clouds_gno_40(ks, ke, nd, r, rs, qsub, cldf)
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
!
! CALL command:
! -------------
!
!     CALL CLOUDS_GNO(ND,R,RS,QSUBGRID,CLDF)
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

  integer, intent(in)    :: ks, ke, nd
  real,    intent(in)    :: r(nd), rs(nd), qsub(nd)
  real,    intent(inout) :: cldf(nd)

! -- lower bound of the PDF of total water:

  real, parameter :: pb = 0.0

! -- parameters controlling the iteration:
! --    nmax    : maximum nb of iterations (hopefully never reached!)
! --    epsilon : accuracy of the numerical resolution (here 2.0%)

  integer, parameter :: nmax    = 5
  real,    parameter :: epsilon = 0.05

! -- gardes-fou:

  real, PARAMETER :: min_mu = 1.e-12
  real, PARAMETER :: min_Q  = 1.e-12

! -- misc:

  INTEGER :: K, n
  REAL    :: mu, qsat, delta, beta 
  REAL    :: xx, dist, fprime
  REAL    :: u, v, erfu, erfv

! real    :: alpha, lambda, kew, skew, sigs

  real, parameter :: pi = ACOS(-1.)
  real, parameter :: ff = 2.0 / sqrt(pi)
  real, parameter :: xx0 = -0.5*sqrt(log(2.))
  real, parameter :: uva = 1.0 / sqrt(2.0) / xx0
  real, parameter :: uvb = 0.5 / sqrt(2.0) * xx0

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  ! -- loop over vertical levels :

  DO K = ks, ke

     mu = R(K)
     mu = MAX(mu,min_mu)
     qsat = RS(K) 
     qsat = MAX(qsat,min_mu)

!===========================================================================
!  If no condensation is produced at the subgrid-scale:
!
!  -> the scheme becomes equivalent to a "large-scale condensation scheme"
!     ie: cldf = H(mu-qsat)
!     where H is the Heaviside function.
!     (in the absence of subgrid-scale condensation, the generalized
!      log-normal PDF becomes equivalent to a gaussian PDF of variance
!      zero, i.e. it becomes equivalent to a Dirac function and the 
!      cumulative distribution function becomes an Heaviside function).
!===========================================================================

     IF ( QSUB(K) .lt. min_Q ) cycle

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

     delta = log(mu/qsat)
        
     if (delta > 0.0) then
        beta = QSUB(K)/mu + 1.0
     else
        beta = QSUB(K)/mu + EXP( -delta)
     endif

     xx = xx0
     u = delta * uva + uvb
     v = delta * uva - uvb

     iter: do n = 1, nmax ! iteration loop

        erfu = 1. - ERF(u)
        erfv = 1. - ERF(v)

        erfv = max(erfv, 1.e-16)
        dist  = erfu / erfv - beta

! -- numerical convergence reached?
           
        if ( ABS(dist/beta) .LT. epsilon .or. n == nmax .and. n > 1) then

           exit iter

        else

           fprime = ff / xx / erfv**2 &
                  * ( erfv*v*EXP(-u*u) - erfu*u*EXP(-v*v) )

           xx = xx - dist/fprime

           u = delta/(xx*sqrt(2.)) + xx/(2.*sqrt(2.))
           v = delta/(xx*sqrt(2.)) - xx/(2.*sqrt(2.))

        endif

     ENDDO iter

! deduce the cloud fraction:

     CLDF(K) = max( 0.5 * erfv, 0.1 )

  ENDDO

end subroutine clouds_gno_40
