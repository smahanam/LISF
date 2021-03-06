!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ecmwfreanal_forcingMod
!BOP
! !MODULE: ecmwfreanal_forcingMod
! 
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the bias-corrected ECMWF atmospheric reanalysis
!  data (Berg et al. 2003). The data is global 0.5 degree dataset in latlon
!  projection, and at 6 hourly intervals. The derived
!  data type {\tt ecmwfreanal\_struc}
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[fmodeltime1]
!    The nearest, previous 6 hour instance of the incoming 
!    data (as a real time). 
!  \item[fmodeltime2]
!    The nearest, next 6 hour instance of the incoming 
!    data (as a real time).
!  \item[remask1d]
!    The data used to remask the input data to the LDT mask. 
!  \item[ecmwfreanaldir]
!    Directory containing the input data
!  \item[emaskfile]
!    File containing the 0.5 deg land-sea mask used in the input data. 
!  \item[elevfile]
!    File with the elevation definition for the input data. 
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n12,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
! 
!  Berg, A.A., J.S. Famiglietti, J.P. Walker, and P.R. Houser, 2003: 
!  Impact of bias correction to reanalysis products on simulations of
!  North American soil moisture and hydrological fluxes. Journal of 
!  Geophysical Research, 108, 4490, DOI: 10.1029/2002JD003334. \newline
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_ECMWFREANAL      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ecmwfreanal_struc

!EOP

  type, public :: ecmwfreanal_type_dec
     real                   :: ts
     integer                :: nc, nr   
     real*8                 :: fmodeltime1,fmodeltime2
     integer                :: remask1d(720*360)
     character*100          :: ecmwfreanaldir
     character*100          :: emaskfile  ! 1/2deg ECMWFREANAL Land-Sea Mask File
     character*100          :: elevfile
     character*50           :: elevtransform
     integer                :: mi
     integer                :: findtime1, findtime2

     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)     
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
  end type ecmwfreanal_type_dec

  type(ecmwfreanal_type_dec), allocatable :: ecmwfreanal_struc(:) 

contains
  
!BOP
!
! !ROUTINE: init_ECMWFREANAL
! \label{init_ECMWFREANAL}
!
! !REVISION HISTORY: 
! 26Jan2004: Sujay Kumar; Init_ial Specification
! 
! !INTERFACE:
  subroutine init_ECMWFREANAL(findex)
! !USES: 
   use LDT_coreMod,    only : LDT_rc, LDT_domain
   use LDT_timeMgrMod, only : LDT_update_timestep
   use LDT_logMod,     only : LDT_logunit, LDT_endrun

   implicit none
! !AGRUMENTS: 
   integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for ECMWFREANAL 
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_ecmwfreanal](\ref{readcrd_ecmwfreanal}) \newline
!     reads the runtime options specified for ECMWFREANAL data
!   \item[readmask\_ecmwfreanal](\ref{readmask_ecmwfreanal}) \newline
!     reads the 0.5 degree land sea mask
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_ecmwfreanal\_elev](\ref{read_ecmwfreanal_elev}) \newline
!    reads the native elevation of the ecmwfreanal data to be used
!    for topographic adjustments to the forcing 
!  \end{description}
!EOP
    integer  :: n 
    real, allocatable :: elev(:,:)
    real, allocatable :: dummy(:,:)
    real     :: gridDesci(20)

    allocate(ecmwfreanal_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"MSG: Initializing ECMWF-Reanalysis forcing grid ... "

!- Readin ldt.config file entries:
    call readcrd_ecmwfreanal(findex)

    do n=1, LDT_rc%nnest
       ecmwfreanal_struc(n)%ts = 6*3600 
       call LDT_update_timestep(LDT_rc, n, ecmwfreanal_struc(n)%ts)
    enddo

  ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex)  = "latlon"
    LDT_rc%met_nc(findex) = 720
    LDT_rc%met_nr(findex) = 360

    ecmwfreanal_struc%nc = 720
    ecmwfreanal_struc%nr = 360

 !- AGRMET Grid description:
    gridDesci(1)  = 0
    gridDesci(2)  = ecmwfreanal_struc(1)%nc
    gridDesci(3)  = ecmwfreanal_struc(1)%nr
    gridDesci(4)  = -89.750
    gridDesci(5)  = -179.750
    gridDesci(6)  = 128
    gridDesci(7)  = 89.750
    gridDesci(8)  = 179.750
    gridDesci(9)  = 0.50
    gridDesci(10) = 0.50
    gridDesci(20) = 64

    LDT_rc%met_gridDesc(findex,1:20) = gridDesci(1:20) 

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9
    LDT_rc%met_ts(findex) = 6*3600
    LDT_rc%met_zterp(findex) = .true.

    do n=1,LDT_rc%nnest
       
       call readmask_ecmwfreanal(n)

       ecmwfreanal_struc(n)%findtime1 = 0 
       ecmwfreanal_struc(n)%findtime2 = 0

       ecmwfreanal_struc(n)%mi = ecmwfreanal_struc(n)%nc * &
                                 ecmwfreanal_struc(n)%nr

     ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

        case( "bilinear" )
          allocate(ecmwfreanal_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesci(:),&
               ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
               ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,&
               ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
               ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221)

        case( "budget-bilinear" )

          allocate(ecmwfreanal_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(ecmwfreanal_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesci(:),&
               ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
               ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,&
               ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
               ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221)

          allocate(ecmwfreanal_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(ecmwfreanal_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input(n, gridDesci(:),&
               ecmwfreanal_struc(n)%n112,ecmwfreanal_struc(n)%n122,&
               ecmwfreanal_struc(n)%n212,ecmwfreanal_struc(n)%n222,&
               ecmwfreanal_struc(n)%w112,ecmwfreanal_struc(n)%w122,&
               ecmwfreanal_struc(n)%w212,ecmwfreanal_struc(n)%w222)

       case default
         write(LDT_logunit,*) " ERR: ONLY 'bilinear' OR 'budget-bilinear' ARE"
         write(LDT_logunit,*) "  CURRENTLY AVAILABLE FOR THE ECMWF-REANALYSIS DATASET"
         write(LDT_logunit,*) " Stopping ... "
         call LDT_endrun
       end select

#if 0
     ! Read in elevation file:
       if( LDT_rc%met_ecor(findex).ne."none" ) then 
          allocate( elev(LDT_rc%lnc(n),LDT_rc%lnr(n)) )
          allocate( dummy(nldas2_struc(n)%nc,nldas2_struc(n)%nr) )

          call read_ecmwfreanal_elev(n,findex, elev, dummy)

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(LDT_domain(n)%gindex(c,r).ne.-1) then
                   LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = elev(c,r)
                endif
             enddo
          enddo
          deallocate( elev, dummy )

       endif
#endif

    enddo   ! End nest loop

  end subroutine init_Ecmwfreanal

end module ecmwfreanal_forcingMod
