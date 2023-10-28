implicit none
include 'mpif.h'

!==========================================================!
!~~~ fundamental constants [all numbers are in SI units]~~~!
!==========================================================!
double precision, parameter :: pi=3.141592653589793d0,c=299792458.0d0
double precision, parameter :: sqrt2=1.414213562373095d0,sqrt3=1.732050807568877d0,sqrt6=2.449489742783178d0
double precision, parameter :: sqrt2pi=2.5066282746310002d0
double precision, parameter :: mu0=4.0d-7*pi,eps0=1.0d0/(c*c*mu0)
double precision, parameter :: h=1.054571628d-34
double precision, parameter :: qe=1.6d-19
double complex, parameter :: Im=(0.0, 1.0)
double precision, parameter :: Hz_to_ev=4.1356691d-15
double precision, parameter :: ev_to_radsec=2.0*pi*2.418d14
double precision, parameter :: ev_to_Joule =1.602176565d-19
double precision, parameter :: ev_to_THz=241.7988407d0
double precision, parameter :: eV_to_wavenumber=8065.541154d0
double precision, parameter :: eV_to_au=1.0d0/27.2117d0
double precision, parameter :: Debye_to_Cm=3.33564d-30
!==========================================================!
!~~~=================== physical grid ==================~~~!
!==========================================================!
logical,parameter, dimension(2) :: periods=(/.true.,.false./) !<--- periodic boundaries in X
                                                              !<--- non-periodic boundaries in Y
integer, parameter :: ndims=2
integer, parameter, dimension(2) :: dims=(/6,8/)  !<--- dims(1) is the number of processors along X
                                                  !<--- dims(2) is the number of processors along Y
                                                  !<--- thus 6x8=48 processors
                                                  !<--- always check if dims(1)*dims(2)=nprocs
double precision, parameter :: dx=1.0d-9,dy=dx

integer, parameter :: Nx=dims(1)*140
integer, parameter :: Ny=dims(2)*150

double precision, parameter :: x0=-dx*Nx/2.0,xM=-x0-dx
double precision, parameter :: y0=-dy*Ny/2.0,yM=-y0-dy

integer, parameter :: Nx_loc=Nx/dims(1),Ny_loc=Ny/dims(2)

double precision x,xM2,y,yM2

!==========================================================!
!~~~================== time iterations =================~~~!
!~~~==================  and FFT stuff  =================~~~!
!==========================================================!
double precision, parameter :: dt=dx/(2.0*c)

double precision, parameter :: t_final= 300.0d-15  !time propagation
integer, parameter :: Nt=int(t_final/dt)          !total # of time steps

double precision, parameter :: dt_eps0=dt/eps0,dt_mu0=dt/mu0

integer N_w
double precision, parameter :: omega_min=0.0,omega_max=4.0
integer resolution_fr
double precision, allocatable :: omega_P(:)
integer, parameter :: Nt_skip=30
integer n_new
integer, parameter :: Nt_new=int(Nt/Nt_skip)
double precision wsave(4*Nt_new+15)

!==========================================================!
!~~~============ incident pulse & detection ============~~~!
!~~~================== TF/SF boundary ==================~~~!
!==========================================================!
integer, parameter ::  ms=dims(2)-1,js=68                               !<=== locaion of the source at y = +450 nm
integer, parameter :: mj1=dims(2)-2,j1=111                  !<=== TF/SF boundary, should be below ms,js y = +350 nm
double precision, parameter :: tau=1.0d-15
double precision, parameter :: E0=1.0 !9.22d28 !E field does not couple why?
double precision, parameter :: omega=ev_to_radsec*2.7d0
double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/)
double precision pulse

integer, parameter :: mwR=dims(2)-1,jwR=41                                !<=== locaion of R at y = +430 nm, outside TF/SF
integer, parameter :: mwT=0,jwT=55                                !<=== locaion of T at y = -430 nm
integer, parameter :: mw_reduce=1                                 !<=== for mpi_reduce
double complex Ex_temp(Nx_loc,Nt_new),Hz_temp(Nx_loc,Nt_new)
double complex Ex_temp_inc(Nx_loc,Nt_new),Hz_temp_inc(Nx_loc,Nt_new)
double precision Ex_w,Hz_w,av_x1,av_x2,av_y
double complex sum1,cTMP1,cTMP2
double precision P_inc
double precision, allocatable :: R_loc(:),T_loc(:),R_global(:),T_global(:),P_inc_loc(:),P_inc_global(:)   !(resolution_fr)

double precision, parameter :: shift_Y=0.0,delta_x=5d-9!-100.2525256456d-9,delta_x=5.252524567d-9
double precision circle1,circle2,circle3,circle4,circle5

!=====================================================================================!
!~~~~================================ glass ======================================~~~~!
!=====================================================================================!
double precision, parameter :: eps_glass=1.4*1.4
double precision, parameter :: dt_eps_glass=dt_eps0/eps_glass
logical F_glass_x(Nx_loc,Ny_loc),F_glass_y(Nx_loc,Ny_loc),F_glass(Nx_loc,Ny_loc)

!=====================================================================================!
!~~~~================================ SU8 ======================================~~~~!
!=====================================================================================!
double precision, parameter :: eps_SU8=eps_glass!1.5854*1.5854
double precision, parameter :: level_SU8=-90.252524567d-9
double precision, parameter :: dt_eps_SU8=dt_eps0/eps_SU8
logical F_SU8_x(Nx_loc,Ny_loc),F_SU8_y(Nx_loc,Ny_loc),F_SU8(Nx_loc,Ny_loc)

!====================================================================================!
!~~~~======================== Drude-Lorentz model for Au ========================~~~~!
!====================================================================================!
!
!~~~ Drude pole ~~~!
!
double precision, parameter :: eps_r=1.0,omegaD=ev_to_radsec*(9.03*0.8718),GammaD=ev_to_radsec*0.053
double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
!
!~~~ Lorentz poles ~~~!
!
double precision, parameter :: omegaP1=ev_to_radsec*0.415
double precision, parameter :: omegaP2=ev_to_radsec*0.830
double precision, parameter :: omegaP3=ev_to_radsec*2.969
double precision, parameter :: omegaP4=ev_to_radsec*4.304
double precision, parameter :: omegaP5=ev_to_radsec*13.32

double precision, parameter :: GammaP1=ev_to_radsec*0.241
double precision, parameter :: GammaP2=ev_to_radsec*0.345
double precision, parameter :: GammaP3=ev_to_radsec*0.870
double precision, parameter :: GammaP4=ev_to_radsec*2.494
double precision, parameter :: GammaP5=ev_to_radsec*2.214

double precision, parameter :: omegaP=ev_to_radsec*9.03

double precision, parameter :: AP1=omegaP*omegaP*0.024
double precision, parameter :: AP2=omegaP*omegaP*0.010
double precision, parameter :: AP3=omegaP*omegaP*0.071
double precision, parameter :: AP4=omegaP*omegaP*0.601
double precision, parameter :: AP5=omegaP*omegaP*4.384

double precision alpha_k(5),beta_k(5),gamma_k(5)

double precision C1,C2,C3,C4
double precision B1_k(5),B2_k(5)

double precision tmpE
double precision tmpPL1,tmpPL2,tmpPL3,tmpPL4,tmpPL5
!~~~~~~~~~~~~~~ X-component ~~~~~~~~~~~~~~~~!
!~~~ Drude current ~~~!
double precision PDx(Nx_loc,Ny_loc)
!~~~ Lorentz currents ~~~!
double precision PLx1(Nx_loc,Ny_loc),PLx1_P(Nx_loc,Ny_loc)
double precision PLx2(Nx_loc,Ny_loc),PLx2_P(Nx_loc,Ny_loc)
double precision PLx3(Nx_loc,Ny_loc),PLx3_P(Nx_loc,Ny_loc)
double precision PLx4(Nx_loc,Ny_loc),PLx4_P(Nx_loc,Ny_loc)
double precision PLx5(Nx_loc,Ny_loc),PLx5_P(Nx_loc,Ny_loc)

!~~~~~~~~~~~~~~ Y-component ~~~~~~~~~~~~~~~~!
!~~~ Drude current ~~~!
double precision PDy(Nx_loc,Ny_loc)
!~~~ Lorentz currents ~~~!
double precision PLy1(Nx_loc,Ny_loc),PLy1_P(Nx_loc,Ny_loc)
double precision PLy2(Nx_loc,Ny_loc),PLy2_P(Nx_loc,Ny_loc)
double precision PLy3(Nx_loc,Ny_loc),PLy3_P(Nx_loc,Ny_loc)
double precision PLy4(Nx_loc,Ny_loc),PLy4_P(Nx_loc,Ny_loc)
double precision PLy5(Nx_loc,Ny_loc),PLy5_P(Nx_loc,Ny_loc)

double precision Ex_P(Nx_loc,Ny_loc),Ey_P(Nx_loc,Ny_loc)

!=============================================================================!
!~~~~============================ metal layer ============================~~~~!
!=============================================================================!
double precision, parameter :: amplitude=100.2525256456d-9
double precision, parameter :: metal_thickness=10.2525256456d-9
logical FBx(Nx_loc,Ny_loc),FBy(Nx_loc,Ny_loc),FB(Nx_loc,Ny_loc)

!=============================================================================!
!~~~~============================ Al2O3 layer ============================~~~~!
!=============================================================================!
double precision, parameter :: eps_Al2O3=1.7650*1.7650
double precision, parameter :: dt_eps_Al2O3=dt_eps0/eps_Al2O3
double precision, parameter :: Al2O3_thickness=20.2525256456d-9 !varies 5 - 25 nm
logical FAl2O3x(Nx_loc,Ny_loc),FAl2O3y(Nx_loc,Ny_loc),FAl2O3(Nx_loc,Ny_loc)

!=============================================================================!
!~~~~======================== molecular layer ============================~~~~!
!=============================================================================!
double precision, parameter :: molecules_thickness=dy*10.0 !10.2525256456d-9 !should be very small but enough to be continous
logical FMx(Nx_loc,Ny_loc),FMy(Nx_loc,Ny_loc),FM(Nx_loc,Ny_loc)

double precision, parameter :: n0=1.0d26                                        !<------ number density
double precision, parameter :: dp=Debye_to_Cm*10.0 !(default)                              !<---- transition dipole
double precision, parameter :: gamma1=2.0*pi/1000.0D-15                     !<---- relaxation rate
double precision, parameter :: gamma_dp=2.0*pi/500.0D-15                    !<---- pure dephasing rate
double precision, parameter :: gamma2=gamma1/2.0+gamma_dp
double precision, parameter :: Omega0=ev_to_radsec*1.578 !(resonance)                         !<---- transition frequency
double precision, parameter :: cx=1.0/(4.0*eps0),cy=cx
double precision, parameter :: dt_epsM=dt/eps0

double precision, parameter :: const=1.0/(h*Omega0)
double precision, parameter :: coupling=2.0*Omega0*dp*dp/(3.0*h)

double complex a_1,a_2

double complex, dimension(3,3) :: Hamiltonian,rho,muy
double precision, dimension(3,3) :: mux
double complex, dimension(Nx_loc,Ny_loc) :: rho11,rho12,rho13
double complex, dimension(Nx_loc,Ny_loc) :: rho22,rho23
double complex, dimension(Nx_loc,Ny_loc) :: rho33
! HTC: population observation
double precision, dimension(Nt/Nt_skip) :: population
double precision, dimension(:), allocatable :: total_population
double precision excited_population, reduced_excited_population
double precision, parameter :: initial_population=0.00000
! HTC: introduce random initial phase 
logical, parameter :: RandomInitialPhase=.false.
double precision phase(3)
double complex exp_phase(3)

double complex, dimension(3,3) :: tmp_P
double complex kk1(3,3),kk2(3,3),kk3(3,3),kk4(3,3)

double complex Omega_plus_n,Omega_minus_n
double complex Omega_plus_n05,Omega_minus_n05
double complex Omega_plus_n1,Omega_minus_n1
double precision, dimension(Nx_loc,Ny_loc) :: dPx,dPy,dPx_old,dPy_old
double precision, dimension(Nx_loc,Ny_loc) :: Px_old,Py_old
double precision Px_loc,Py_loc
double precision Px_av,Py_av
double precision avX_n,avX_n1,avX_n05,avY_n,avY_n05,avY_n1
double precision dPx_send(Ny_loc),dPx_recv(Ny_loc),ExM_send(Ny_loc,2),ExM_recv(Ny_loc,2)
double precision dPy_send(Nx_loc),dPy_recv(Nx_loc),EyM_send(Nx_loc,2),EyM_recv(Nx_loc,2)
double precision Ex_old(Nx_loc,Ny_loc),Ey_old(Nx_loc,Ny_loc)

!=============================================================================!
!~~~~============================= PMMA layer ============================~~~~!
!=============================================================================!
double precision, parameter :: eps_PMMA=1.4881*1.4881
double precision, parameter :: dt_eps_PMMA=dt_eps0/eps_PMMA
double precision, parameter :: PMMA_thickness=50.2525256456d-9 !varies 5 - 25 nm
logical FPMMAx(Nx_loc,Ny_loc),FPMMAy(Nx_loc,Ny_loc),FPMMA(Nx_loc,Ny_loc)

!============================================================================!
!~~~~==================== =electromagnetic field= =======================~~~~!
!============================================================================!
double precision Ex(Nx_loc,Ny_loc),Ey(Nx_loc,Ny_loc),Hz(Nx_loc,Ny_loc)
double precision Ex_send(Nx_loc),Ex_recv(Nx_loc),Ey_send(Ny_loc),Ey_recv(Ny_loc)
double precision Hzx_send(Nx_loc),Hzx_recv(Nx_loc),Hzy_send(Ny_loc),Hzy_recv(Ny_loc)

double precision Ex_inc(Ny_loc),Hz_inc(Ny_loc)
double precision Ex_inc_send,Ex_inc_recv,Hz_inc_send,Hz_inc_recv

!============================================================================!
!~~~~============================ cpml stuff  ===========================~~~~!
!============================================================================!
double precision, parameter :: eps_delectric=1.0
integer, parameter :: npml=50,m=3,ma=1                        !<--- npml should be less than Nx_loc&Ny_loc
double precision, parameter :: alphaCPML=0.05,kappaCPML=5.0   !<--- optimal parameters per Taflove
double precision sigmaCPML
double precision be_y(npml),ce_y(npml),alphae_y(npml),sige_y(npml),kappae_y(npml)
double precision bh_y(npml-1),ch_y(npml-1),alphah_y(npml-1),sigh_y(npml-1),kappah_y(npml-1)
double precision den_ey(Ny_loc),den_hy(Ny_loc)
double precision, parameter :: den_ex=1.0/dx,den_hx=1.0/dx
double precision psi_Hzy_inc(npml-1),psi_Exy_inc(npml)
double precision psi_Hzy(Nx_loc,npml-1),psi_Exy(Nx_loc,npml)

double precision sigmaCPML_air
double precision be_y_air(npml),ce_y_air(npml),alphae_y_air(npml),sige_y_air(npml),kappae_y_air(npml)
double precision bh_y_air(npml-1),ch_y_air(npml-1),alphah_y_air(npml-1),sigh_y_air(npml-1),kappah_y_air(npml-1)
double precision den_ey_air(Ny_loc),den_hy_air(Ny_loc)

!============================================================================!
!~~~~========================== some variables  =========================~~~~!
!============================================================================!
integer i,j,i_global,j_global,n,nn,jj,rank_tmp1,rank_tmp2,coords_tmp(ndims),k
double precision tmp1,tmp2,tmp,t,cpu1,cpu2,temp_send(2),temp_recv(2)

!============================================================================!
!~~~~====================== IO PART,    OUTPUT FIELD=====================~~~~!
!============================================================================!
character(len=24) :: filename, folderName
character(len=10) :: filehead
integer io_coords(ndims), folderStatus, n_print

!============================================================================!
!~~~~====================== MPI part, don't modify  =====================~~~~!
!============================================================================!
integer ierr,nprocs,myrank,nbr_left,nbr_right,nbr_up,nbr_down
integer sendtag1,recvtag1,sendtag2,recvtag2,sendtag3,recvtag3,sendtag4,recvtag4
integer sendtag5,recvtag5,sendtag6,recvtag6,sendtag7,recvtag7,sendtag8,recvtag8
integer old_comm,cartesian_comm
logical, parameter :: reorder=.false.   !<--- probably shouldn't reoder processors
integer coords(ndims)
integer :: istatus(MPI_STATUS_SIZE)     !<--- needed for send/recv calls

!~~~~=== setting and initializing MPI ===~~~~!
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)      !<--- nprocs=total number of processors
call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)      !<--- myrank now is set to unique integer for each processor

old_comm=MPI_COMM_WORLD

!~~~~=== creating Cartesian topology ===~~~~!
call mpi_cart_create( &
      old_comm,         &        !<--- original communicator
      ndims,            &        !<--- ndims = integer, number of dimensions
      dims,             &        !<--- dims(ndims) = number of processors in each dimension
      periods,          &        !<--- periods(ndims) = logical array defining boundary conditions
      reorder,          &        !<--- reorder = logical, reorder or not processors [set to .false.]
      cartesian_comm,   &        !<--- new topology defined
      ierr              &        !<--- ierror = integer, error
                     )
!<---- now I should have a Cartesian grid, first index is for X and the second one is for Y

!~~~~=== get Cartesian coordinates ===~~~~!
call mpi_cart_coords( &
      cartesian_comm, &
      myrank,         &
      ndims,          &
      coords,         &     !<--- coordinates(ndims) = gives integers identifying local coordinates of a given block for myrank
      ierr            &     !<--- coords(1) = x, coords(2) = y
                     )

!~~~~=== map of neighbors for X ===~~~~!
call mpi_cart_shift(  &
      cartesian_comm, &
      0,              &     !<~~~ direction along X
      1,              &     !<~~~ displacement along X
      nbr_left,       &     !<~~~ neighbor to the left
      nbr_right,      &     !<~~~ neighbor to the right
      ierr            &
                    )
if(nbr_left<0)then
  nbr_left = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
endif
if(nbr_right<0)then
  nbr_right = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
endif

!~~~~=== map of neighbors for Y ===~~~~!
call mpi_cart_shift(  &
      cartesian_comm, &
      1,              &     !<~~~ direction along Y
      1,              &     !<~~~ displacement along Y
      nbr_down,       &     !<~~~ neighbor down
      nbr_up,         &     !<~~~ neighbor up
      ierr            &
                    )
if(nbr_down<0)then
  nbr_down = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
endif
if(nbr_up<0)then
  nbr_up = MPI_PROC_NULL      !<~~~ making sure we don't send/receive outside of the grid
endif

! HTC: allocate total population
allocate(total_population(Nt/Nt_skip * nprocs))

!~~~ Drude-Lorentz parameters for Ampere Law updates~~~!
alpha_k(1)=(2.0-omegaP1*omegaP1*dt*dt)/(1+0.5*GammaP1*dt)
alpha_k(2)=(2.0-omegaP2*omegaP2*dt*dt)/(1+0.5*GammaP2*dt)
alpha_k(3)=(2.0-omegaP3*omegaP3*dt*dt)/(1+0.5*GammaP3*dt)
alpha_k(4)=(2.0-omegaP4*omegaP4*dt*dt)/(1+0.5*GammaP4*dt)
alpha_k(5)=(2.0-omegaP5*omegaP5*dt*dt)/(1+0.5*GammaP5*dt)

beta_k(1)=(GammaP1*dt-2.0)/(GammaP1*dt+2.0)
beta_k(2)=(GammaP2*dt-2.0)/(GammaP2*dt+2.0)
beta_k(3)=(GammaP3*dt-2.0)/(GammaP3*dt+2.0)
beta_k(4)=(GammaP4*dt-2.0)/(GammaP4*dt+2.0)
beta_k(5)=(GammaP5*dt-2.0)/(GammaP5*dt+2.0)

gamma_k(1)=eps0*AP1*dt/(GammaP1*dt+2.0)
gamma_k(2)=eps0*AP2*dt/(GammaP2*dt+2.0)
gamma_k(3)=eps0*AP3*dt/(GammaP3*dt+2.0)
gamma_k(4)=eps0*AP4*dt/(GammaP4*dt+2.0)
gamma_k(5)=eps0*AP5*dt/(GammaP5*dt+2.0)

C1=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k))
C2=0.5*SUM(gamma_k)/(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k))
C3=1/(dx*(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k)))
C4=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k))

do k=1,5
 B1_k(k)=0.5*(1.0+alpha_k(k))/(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k))
 B2_k(k)=0.5*beta_k(k)/(eps_r*eps0/dt+0.5*A2+0.5*SUM(gamma_k))
enddo



!~~~~=== setting geometry ===~~~!
F_glass_x=.false.
F_glass_y=.false.
F_glass=.false.

F_SU8_x=.false.
F_SU8_y=.false.
F_SU8=.false.

FBx=.false.
FBy=.false.
FB=.false.

FAl2O3x=.false.
FAl2O3y=.false.
FAl2O3=.false.

FMx=.false.
FMy=.false.
FM=.false.

FPMMAx=.false.
FPMMAy=.false.
FPMMA=.false.


do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
    if( &
        (y-shift_Y)<level_SU8 &
      )then
      F_glass_x(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
    if( &
        (y-shift_Y)<level_SU8 &
      )then
      F_glass_y(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
    if( &
        (y-shift_Y)<level_SU8 &
      )then
      F_glass(i,j)=.true.
    endif
  enddo
enddo





do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
    if( &
        (circle1<1.0).or.((y-shift_Y)<0.0.and.(y-shift_Y)>level_SU8) &
      )then
      F_SU8_x(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
    if( &
        (circle1<1.0).or.((y-shift_Y)<0.0.and.(y-shift_Y)>level_SU8) &
      )then
      F_SU8_y(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
    if( &
        (circle1<1.0).or.((y-shift_Y)<0.0.and.(y-shift_Y)>level_SU8) &
      )then
      F_SU8(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
      !circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
    if( &
        (circle1<1.0) &
      )then
      FBx(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
      !circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
    if( &
        (circle1<1.0) &
      )then
      FBy(i,j)=.true.
    endif
  enddo
enddo

do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle1=(x/amplitude)**2+((y-shift_Y)/amplitude)**2
     ! circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
    if( &
        (circle1<1.0) &
      )then
      FB(i,j)=.true.
    endif
  enddo
enddo


!Al2O3
FAl2O3x=.false. !xM2,y grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
      circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
    if( &
    (((circle2>1.0).and.((y-shift_Y-metal_thickness)>0.0)).and. &
     ((circle3<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness)<0.0))).and. &
     (.not.(FBx(i,j))) &
      )then
      FAl2O3x(i,j)=.true.
     else
      FAl2O3x(i,j)=.false.
    endif
  enddo
enddo

FAl2O3y=.false. !x,yM2 grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
    circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
    circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
  if( &
  (((circle2>1.0).and.((y-shift_Y-metal_thickness)>0.0)).and. &
   ((circle3<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness)<0.0))).and. &
   (.not.(FBy(i,j))) &
      )then
      FAl2O3y(i,j)=.true.
     else
      FAl2O3y(i,j)=.false.
    endif
  enddo
enddo

FAl2O3=.false. !x,yM2 grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
    circle2=(x/(amplitude+delta_x))**2+((y-shift_Y-metal_thickness)/(amplitude+delta_x))**2
    circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
  if( &
  (((circle2>1.0).and.((y-shift_Y-metal_thickness)>0.0)).and. &
   ((circle3<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness)<0.0))).and. &
   (.not.(FB(i,j))) &
      )then
      FAl2O3(i,j)=.true.
     else
      FAl2O3(i,j)=.false.
    endif
  enddo
enddo


!molecules
FMx=.false. !xM2,y grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
      circle4=(x/(amplitude+3.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)/(amplitude+3.0*delta_x))**2
    if( &
    (((circle3>1.0).and.((y-shift_Y-metal_thickness-Al2O3_thickness)>0.0)).and. &
      ((circle4<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)<0.0))).and. &
     (.not.(FAl2O3x(i,j))) &
      )then
      FMx(i,j)=.true.
     else
      FMx(i,j)=.false.
    endif
  enddo
enddo

FMy=.false. !x,yM2 grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
      circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
      circle4=(x/(amplitude+3.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)/(amplitude+3.0*delta_x))**2
    if( &
    (((circle3>1.0).and.((y-shift_Y-metal_thickness-Al2O3_thickness)>0.0)).and. &
      ((circle4<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)<0.0))).and. &
     (.not.(FAl2O3y(i,j))) &
      )then
      FMy(i,j)=.true.
     else
      FMy(i,j)=.false.
    endif
  enddo
enddo


FM=.false. !x,y grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle3=(x/(amplitude+2.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness)/(amplitude+2.0*delta_x))**2
      circle4=(x/(amplitude+3.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)/(amplitude+3.0*delta_x))**2
    if( &
    (((circle3>1.0).and.((y-shift_Y-metal_thickness-Al2O3_thickness)>0.0)).and. &
      ((circle4<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)<0.0))).and. &
     (.not.(FAl2O3(i,j))) &
      )then
      FM(i,j)=.true.
     else
      FM(i,j)=.false.
    endif
  enddo
enddo


!PMMA
FPMMAx=.false. !xM2,y grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)+dx/2.0
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)
      circle4=(x/(amplitude+3.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)/(amplitude+3.0*delta_x))**2
      circle5=(x/(amplitude+4.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness-PMMA_thickness)/(amplitude+4.0*delta_x))**2
    if( &
    (((circle4>1.0).and.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)>0.0)).and. &
     ((circle5<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness-PMMA_thickness)<0.0))).and. &
     (.not.(FMx(i,j))) &
      )then
      FPMMAx(i,j)=.true.
     else
      FPMMAx(i,j)=.false.
    endif
  enddo
enddo

FPMMAy=.false. !x,yM2 grid
do i=1,Nx_loc
  i_global=coords(1)*Nx_loc+i
  x=x0+dx*(i_global-1)
  do j=1,Ny_loc
    j_global=coords(2)*Ny_loc+j
    y=y0+dy*(j_global-1)+dy/2.0
      circle4=(x/(amplitude+3.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)/(amplitude+3.0*delta_x))**2
      circle5=(x/(amplitude+4.0*delta_x))**2+((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness-PMMA_thickness)/(amplitude+4.0*delta_x))**2
    if( &
    (((circle4>1.0).and.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness)>0.0)).and. &
     ((circle5<1.0).or.((y-shift_Y-metal_thickness-Al2O3_thickness-molecules_thickness-PMMA_thickness)<0.0))).and. &
     (.not.(FMy(i,j))) &
      )then
      FPMMAy(i,j)=.true.
     else
      FPMMAy(i,j)=.false.
    endif
  enddo
enddo
F_glass_x=.false.
F_glass_y=.false.
F_glass=.false.

F_SU8_x=.false.
F_SU8_y=.false.
F_SU8=.false.

!FBx=.false.
!FBy=.false.
!FB=.false.

FAl2O3x=.false.
FAl2O3y=.false.
FAl2O3=.false.

FMx=.false.
FMy=.false.
FM=.false.

FPMMAx=.false.
FPMMAy=.false.
FPMMA=.false.
!~~~~~~~~~~~~~~ end geometry ~~~~~~~~~~~~~~~~~~~~~!

! !HTC: if we want to remove molecule and PMMA
if (molecules_thickness<dx) then
  FMx=.false.
  FMy=.false.
  FM=.false.
endif
if (PMMA_thickness<dx) then
  FPMMAx=.false.
  FPMMAy=.false.
  FPMMA=.false.
endif 

!<~~~~=== setting Fourier stuff ===~~~~!
N_w=0
do n=2,Nt_new/2
    tmp=dble(n-1)/(dt*Nt_skip*Nt_new)
    tmp=Hz_to_ev*tmp
    if((tmp>=omega_min).and.(tmp<=omega_max)) then
        N_w=N_w+1
    endif
enddo
resolution_fr=N_w

! if (myrank==0) then
!     write(*, *) 'frequency resolution', resolution_fr, '# steps', Nt
! endif

allocate(omega_P(resolution_fr))
allocate(R_loc(resolution_fr))
allocate(T_loc(resolution_fr))
allocate(R_global(resolution_fr))
allocate(T_global(resolution_fr))
allocate(P_inc_loc(resolution_fr))
allocate(P_inc_global(resolution_fr))

omega_P=0.0
Ex_temp=0.0
Hz_temp=0.0
Ex_temp_inc=0.0
Hz_temp_inc=0.0
R_loc=0.0
T_loc=0.0
R_global=0.0
T_global=0.0
P_inc_loc=0.0
P_inc_global=0.0

nn=0
do n=2,Nt/2
    tmp=dble(n-1)/(dt*Nt)
    tmp=Hz_to_ev*tmp
    if((tmp>=omega_min).and.(tmp<=omega_max))then
        nn=nn+1
        omega_P(nn)=ev_to_radsec*tmp
    endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!================== CPML vectors ==================!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sigmaCPML=0.8*(m+1)/(dx*(mu0/eps0*eps_glass)**0.5)

do j=1,npml
 sige_y(j)=sigmaCPML*((npml-j)/(npml-1.0))**m
 alphae_y(j)=alphaCPML*((j-1)/(npml-1.0))**ma
 kappae_y(j)=1.0+(kappaCPML-1.0)*((npml-j)/(npml-1.0))**m
 be_y(j)=exp(-(sige_y(j)/kappae_y(j)+alphae_y(j))*dt/eps0)
 if( &
    (sige_y(j)==0.0).and.&
    (alphae_y(j)==0.0).and. &
    (j==npml) &
   )then
   ce_y(j)=0.0
  else
   ce_y(j)=sige_y(j)*(be_y(j)-1.0)/(sige_y(j)+kappae_y(j)*alphae_y(j))/kappae_y(j)
 endif
enddo

do j=1,npml-1
 sigh_y(j)=sigmaCPML*((npml-j-0.5)/(npml-1.0))**m
 alphah_y(j)=alphaCPML*((j-0.5)/(npml-1.0))**ma
 kappah_y(j)=1.0+(kappaCPML-1.0)*((npml-j-0.5)/(npml-1.0))**m
 bh_y(j)=exp(-(sigh_y(j)/kappah_y(j)+alphah_y(j))*dt/eps0)
 ch_y(j)=sigh_y(j)*(bh_y(j)-1.0)/(sigh_y(j)+kappah_y(j)*alphah_y(j))/kappah_y(j)
enddo

den_hy=1.0/dy
if(coords(2)==0)then                    !<~~~ lower y-boundary
  do j=1,Ny_loc
   if(j<=(npml-1))then
    den_hy(j)=1.0/(kappah_y(j)*dy)
   endif
  enddo
elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary
  jj=npml-1
  do j=1,(Ny_loc-1)
   if(j>=(Ny_loc+1-npml))then
     den_hy(j)=1.0/(kappah_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

den_ey=1.0/dy
if(coords(2)==0)then                    !<~~~ lower y-boundary
  do j=1,Ny_loc
   if(j<=npml)then
    den_ey(j)=1.0/(kappae_y(j)*dy)
   endif
  enddo
 elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary
  jj=npml
  do j=1,(Ny_loc-1)
   if(j>=(Ny_loc+1-npml))then
     den_ey(j)=1.0/(kappae_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

!<===========================================-----     -----===========================================>!
!<===========================================-----#air#-----===========================================>!
!<===========================================-----     -----===========================================>!
sigmaCPML_air=0.8*(m+1)/(dx*(mu0/eps0)**0.5)

do j=1,npml
 sige_y_air(j)=sigmaCPML_air*((npml-j)/(npml-1.0))**m
 alphae_y_air(j)=alphaCPML*((j-1)/(npml-1.0))**ma
 kappae_y_air(j)=1.0+(kappaCPML-1.0)*((npml-j)/(npml-1.0))**m
 be_y_air(j)=exp(-(sige_y_air(j)/kappae_y_air(j)+alphae_y_air(j))*dt/eps0)
 if( &
    (sige_y_air(j)==0.0).and.&
    (alphae_y_air(j)==0.0).and. &
    (j==npml) &
   )then
   ce_y_air(j)=0.0
  else
   ce_y_air(j)=sige_y_air(j)*(be_y_air(j)-1.0)/(sige_y_air(j)+kappae_y_air(j)*alphae_y_air(j))/kappae_y_air(j)
 endif
enddo

do j=1,npml-1
 sigh_y_air(j)=sigmaCPML_air*((npml-j-0.5)/(npml-1.0))**m
 alphah_y_air(j)=alphaCPML*((j-0.5)/(npml-1.0))**ma
 kappah_y_air(j)=1.0+(kappaCPML-1.0)*((npml-j-0.5)/(npml-1.0))**m
 bh_y_air(j)=exp(-(sigh_y_air(j)/kappah_y_air(j)+alphah_y_air(j))*dt/eps0)
 ch_y_air(j)=sigh_y_air(j)*(bh_y_air(j)-1.0)/(sigh_y_air(j)+kappah_y_air(j)*alphah_y_air(j))/kappah_y_air(j)
enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!=========~~~ end of CPML ==========~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!~~~~~ dipole operator ~~~~~~~~!
mux=0.0
muy=(0.0,0.0)

mux(1,2)=-dp/sqrt6
mux(1,3)=dp/sqrt6
mux(2,1)=-dp/sqrt6
mux(3,1)=dp/sqrt6

muy(1,2)=Im*dp/sqrt6
muy(1,3)=Im*dp/sqrt6
muy(2,1)=-Im*dp/sqrt6
muy(3,1)=-Im*dp/sqrt6

Hamiltonian=(0.0,0.0)                           !<--- this [and below] is Hamiltonian/hbar
Hamiltonian(2,2)=Omega0
Hamiltonian(3,3)=Omega0

rho11=(0.0,0.0)
rho12=(0.0,0.0)
rho13=(0.0,0.0)
rho22=(0.0,0.0)
rho23=(0.0,0.0)
rho33=(0.0,0.0)

! HTC: initialize density matrix 
 do j=1,Ny_loc
  do i=1,Nx_loc
    if(FM(i,j))then
      !  rho11(i,j)=(1.0,0.0)
      if(RandomInitialPhase)then
        call RANDOM_SEED()
        call RANDOM_NUMBER(phase)
        exp_phase(1) = (1.0,0.0)*cos(2.0*pi*(phase(1)-phase(2))) + (0.0,1.0)*sin(2.0*pi*(phase(1)-phase(2)))
        exp_phase(2) = (1.0,0.0)*cos(2.0*pi*(phase(1)-phase(3))) + (0.0,1.0)*sin(2.0*pi*(phase(1)-phase(3)))
        exp_phase(3) = (1.0,0.0)*cos(2.0*pi*(phase(2)-phase(3))) + (0.0,1.0)*sin(2.0*pi*(phase(2)-phase(3)))
      else
        exp_phase(1) = (1.0,0.0)
        exp_phase(2) = (1.0,0.0)
        exp_phase(3) = (1.0,0.0)
      endif
    
      rho11(i,j)=(1.0,0.0)*(1.0-initial_population)
      rho22(i,j)=(1.0,0.0)*initial_population/2
      rho33(i,j)=(1.0,0.0)*initial_population/2
      rho12(i,j)=(1.0,0.0)*sqrt(rho11(i,j)*rho22(i,j))*exp_phase(1)
      rho13(i,j)=(1.0,0.0)*sqrt(rho11(i,j)*rho33(i,j))*exp_phase(2)
      rho23(i,j)=(1.0,0.0)*sqrt(rho22(i,j)*rho33(i,j))*exp_phase(3)
	  else
      rho11(i,j)=(0.0,0.0)
    endif
  enddo
 enddo


Ex_inc=0.0
Hz_inc=0.0

psi_Hzy_inc=0.0
psi_Exy_inc=0.0

Ex_inc_send=0.0
Ex_inc_recv=0.0
Hz_inc_send=0.0
Hz_inc_recv=0.0

Ex=0.0
Ey=0.0
Hz=0.0

Ex_P=0.0
Ey_P=0.0

tmpE=0.0
tmpPL1=0.0
tmpPL2=0.0
tmpPL3=0.0
tmpPL4=0.0
tmpPL5=0.0

PDx=0.0D0
PLx1=0.0
PLx1_P=0.0
PLx2=0.0
PLx2_P=0.0
PLx3=0.0
PLx3_P=0.0
PLx4=0.0
PLx4_P=0.0
PLx5=0.0
PLx5_P=0.0

PDy=0.0
PLy1=0.0
PLy1_P=0.0
PLy2=0.0
PLy2_P=0.0
PLy3=0.0
PLy3_P=0.0
PLy4=0.0
PLy4_P=0.0
PLy5=0.0
PLy5_P=0.0

psi_Hzy=0.0
psi_Exy=0.0

Ex_send=0.0
Ex_recv=0.0
Hzx_send=0.0
Hzx_recv=0.0
Hzy_send=0.0

if(myrank==0)then
 call cpu_time(cpu1)
endif

n_print=INT(Nt*0.995)
!~~~~=== time iterations ===~~~!
do n=1,Nt
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::!<--- 1D incident field propagation --->!:::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
   t=dt*dble(n)
   if(t<=tau)then
     pulse=aBH(1)+ &
  				 aBH(2)*cos(2.0*pi*t/tau)+ &
  				 aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
  				 aBH(4)*cos(2.0*pi*3.0*t/tau)
    else
     pulse=0.0
   endif

   pulse=E0*sin(omega*t) ! cw

!<=======================================~~~~ Hz incident ~~~~=======================================!
  !<~~~ first we need to send Ex_inc(1) down ~~~
  Ex_inc_send=Ex_inc(1)
  sendtag1=100+n             !<~~~~ some tag integers
  recvtag1=sendtag1           !...
  call mpi_sendrecv(Ex_inc_send,          &     !<=== sending
                    1,                    &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_down,             &     !<=== where sending
                    sendtag1,             &     !<=== sending tag
                    Ex_inc_recv,          &     !<=== receiving
                    1,                    &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_up,               &     !<=== receiving from where
                    recvtag1,             &     !<=== receiving tag
                    cartesian_comm,       &     !<=== handle of Cartesian coordinates
                    istatus,              &     !<=== istatus
                    ierr)                       !<=== error code

    do j=1,Ny_loc
      if(j==Ny_loc)then
          tmp1=Ex_inc_recv
        else
          tmp1=Ex_inc(j+1)
      endif
      Hz_inc(j)=Hz_inc(j)+dt_mu0*(tmp1-Ex_inc(j))*den_hy(j)
    enddo
  !<~~~~ now apply PML at the bottom and at the top all in air
  if(coords(2)==0)then                      !<~~~ lower y-boundary
      do j=1,npml-1
        psi_Hzy_inc(j)=bh_y_air(j)*psi_Hzy_inc(j)+ch_y_air(j)*(Ex_inc(j+1)-Ex_inc(j))/dy
        Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_inc(j)
      enddo
    elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary
     jj=npml-1
      do j=Ny_loc+1-npml,Ny_loc-1
        psi_Hzy_inc(jj)=bh_y_air(jj)*psi_Hzy_inc(jj)+ch_y_air(jj)*(Ex_inc(j+1)-Ex_inc(j))/dy
        Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_inc(jj)
        jj=jj-1
      enddo
  endif

     ! scattered/total field updates
      if(coords(2)==mj1)then
       do i=1,Nx_loc
        Ex(i,j1)=Ex(i,j1)+dt_eps0*Hz_inc(j1)/dy
       enddo
      endif

  !<=======================================~~~~ Ex incident ~~~~=======================================!
  !<~~~ first we need to send Hz_inc(Ny_loc) up ~~~
  Hz_inc_send=Hz_inc(Ny_loc)
  sendtag2=200+n             !<~~~~ some tag integers
  recvtag2=sendtag2           !...
  call mpi_sendrecv(Hz_inc_send,          &     !<=== sending
                    1,                    &     !<=== size
                    mpi_double_precision, &     !<=== type
                    nbr_up,               &     !<=== where sending
                    sendtag2,             &     !<=== sending tag
                    Hz_inc_recv,          &     !<=== receiving
                    1,                    &     !<=== size
                    mpi_double_precision, &     !<=== type
                    nbr_down,             &     !<=== receiving from where
                    recvtag2,             &     !<=== receiving tag
                    cartesian_comm,       &     !<=== handle of Cartesian coordinates
                    istatus,              &     !<=== istatus
                    ierr)                       !<=== error code

  do j=1,Ny_loc
   if((coords(2)==ms).and.(j==js))then !laser pulse
     Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))/dy+ &
  	                    pulse
    elseif(j==1)then
 	    Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc_recv)*den_ey(j)
    else
      Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
   endif
  enddo

  !<~~~~ now apply PML at the bottom and at the top
  if(coords(2)==0)then                      !<~~~ lower y-boundary
      do j=2,npml
        psi_Exy_inc(j)=be_y_air(j)*psi_Exy_inc(j)+ce_y_air(j)*(Hz_inc(j)-Hz_inc(j-1))/dy
        Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_inc(j)
      enddo
    elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary
      jj=npml
      do j=Ny_loc+1-npml,Ny_loc-1
       psi_Exy_inc(jj)=be_y_air(jj)*psi_Exy_inc(jj)+ce_y_air(jj)*(Hz_inc(j)-Hz_inc(j-1))/dy
       Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_inc(jj)
       jj=jj-1
      enddo
  endif


    ! scattered/total field updates
     if(coords(2)==mj1)then
      do i=1,Nx_loc
       Hz(i,j1)=Hz(i,j1)+dt_mu0*Ex_inc(j1)/dy
      enddo
     endif

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::!<--- total 2D field propagation --->!::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  !<~~~== first we need to send Ey(1,1:Ny_loc)-->left & Ex(1:Nx_loc,1)-->down ==~~~
  Ey_send(1:Ny_loc)=Ey(1,1:Ny_loc)
  sendtag3=300+n              !<~~~~ some tag integers
  recvtag3=sendtag3           !...
  call mpi_sendrecv(Ey_send,              &     !<=== sending
                    Ny_loc,               &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_left,             &     !<=== where sending
                    sendtag3,             &     !<=== sending tag
                    Ey_recv,              &     !<=== receiving
                    Ny_loc,               &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_right,            &     !<=== receiving from where
                    recvtag3,             &     !<=== receiving tag
                    cartesian_comm,       &     !<=== handle of Cartesian coordinates
                    istatus,              &     !<=== istatus
                    ierr)                       !<=== error code

  Ex_send(1:Nx_loc)=Ex(1:Nx_loc,1)
  sendtag4=400+n              !<~~~~ some tag integers
  recvtag4=sendtag4           !...
  call mpi_sendrecv(Ex_send,              &     !<=== sending
                    Nx_loc,               &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_down,             &     !<=== where sending
                    sendtag4,             &     !<=== sending tag
                    Ex_recv,              &     !<=== receiving
                    Nx_loc,               &     !<=== the size
                    mpi_double_precision, &     !<=== type
                    nbr_up,               &     !<=== receiving from where
                    recvtag4,             &     !<=== receiving tag
                    cartesian_comm,       &     !<=== handle of Cartesian coordinates
                    istatus,              &     !<=== istatus
                    ierr)                       !<=== error code

  do i=1,Nx_loc
    do j=1,Ny_loc
      if(i==Nx_loc)then
          tmp1=Ey_recv(j)
        else
          tmp1=Ey(i+1,j)
      endif

      if(j==Ny_loc)then
          tmp2=Ex_recv(i)
        else
          tmp2=Ex(i,j+1)
      endif

      Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-tmp1)*den_hx+ &
 			                        (tmp2-Ex(i,j))*den_hy(j))
    enddo
  enddo

  !<~~~~ now apply PML at the bottom and at the top
  if(coords(2)==0)then                      !<~~~ lower y-boundary in air
    do i=1,Nx_loc
      do j=1,npml-1
        psi_Hzy(i,j)=bh_y_air(j)*psi_Hzy(i,j)+ch_y_air(j)*(Ex(i,j+1)-Ex(i,j))/dy
        Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy(i,j)
      enddo
    enddo
  elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary in air
    do i=1,Nx_loc
     jj=npml-1
      do j=Ny_loc+1-npml,Ny_loc-1
        psi_Hzy(i,jj)=bh_y_air(jj)*psi_Hzy(i,jj)+ch_y_air(jj)*(Ex(i,j+1)-Ex(i,j))/dy
        Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy(i,jj)
        jj=jj-1
      enddo
    enddo
  endif

        Ex_old=Ex
        Ey_old=Ey

   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
   !
   dPx_send(1:Ny_loc)=dPx(1,1:Ny_loc)+dPx_old(1,1:Ny_loc)
   sendtag3=300+n              !<~~~~ some tag integers
   recvtag3=sendtag3           !...
   call mpi_sendrecv(dPx_send,             &     !<=== sending
                     Ny_loc,               &     !<=== the size
                     mpi_double_precision, &     !<=== type
                     nbr_left,             &     !<=== where sending
                     sendtag3,             &     !<=== sending tag
                     dPx_recv,             &     !<=== receiving
                     Ny_loc,               &     !<=== the size
                     mpi_double_precision, &     !<=== type
                     nbr_right,            &     !<=== receiving from where
                     recvtag3,             &     !<=== receiving tag
                     cartesian_comm,       &     !<=== handle of Cartesian coordinates
                     istatus,              &     !<=== istatus
                     ierr)                       !<=== error code

   !<~~~== first we need to send Hz(1:Nx_loc,Ny_loc)-->up ==~~~
   Hzx_send(1:Nx_loc)=Hz(1:Nx_loc,Ny_loc)
   sendtag5=500+n             !<~~~~ some tag integers
   recvtag5=sendtag5           !...
   call mpi_sendrecv(Hzx_send,             &     !<=== sending
                     Nx_loc,               &     !<=== the size
                     mpi_double_precision, &     !<=== type
                     nbr_up,               &     !<=== where sending
                     sendtag5,             &     !<=== sending tag
                     Hzx_recv,             &     !<=== receiving
                     Nx_loc,               &     !<=== the size
                     mpi_double_precision, &     !<=== type
                     nbr_down,             &     !<=== receiving from where
                     recvtag5,             &     !<=== receiving tag
                     cartesian_comm,       &     !<=== handle of Cartesian coordinates
                     istatus,              &     !<=== istatus
                     ierr)                       !<=== error code

   do i=1,Nx_loc
     do j=1,Ny_loc
       if(j==1)then
           tmp1=Hzx_recv(i)
         else
           tmp1=Hz(i,j-1)
       endif

       if(i==Nx_loc)then
           tmp2=dPx_recv(j)
         else
           tmp2=dPx(i+1,j)+dPx_old(i+1,j)
       endif

       if(FBx(i,j))then
         tmpE=C1*Ex(i,j)+C2*Ex_P(i,j)+C3*(Hz(i,j)-tmp1)- &
              C4*PDx(i,j)- &
              (B1_k(1)*PLx1(i,j)+B2_k(1)*PLx1_P(i,j)+ &
               B1_k(2)*PLx2(i,j)+B2_k(2)*PLx2_P(i,j)+ &
               B1_k(3)*PLx3(i,j)+B2_k(3)*PLx3_P(i,j)+ &
               B1_k(4)*PLx4(i,j)+B2_k(4)*PLx4_P(i,j)+ &
               B1_k(5)*PLx5(i,j)+B2_k(5)*PLx5_P(i,j))

         PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))

         tmpPL1=alpha_k(1)*PLx1(i,j)+beta_k(1)*PLx1_P(i,j)+gamma_k(1)*(tmpE-Ex_P(i,j))
         tmpPL2=alpha_k(2)*PLx2(i,j)+beta_k(2)*PLx2_P(i,j)+gamma_k(2)*(tmpE-Ex_P(i,j))
         tmpPL3=alpha_k(3)*PLx3(i,j)+beta_k(3)*PLx3_P(i,j)+gamma_k(3)*(tmpE-Ex_P(i,j))
         tmpPL4=alpha_k(4)*PLx4(i,j)+beta_k(4)*PLx4_P(i,j)+gamma_k(4)*(tmpE-Ex_P(i,j))
         tmpPL5=alpha_k(5)*PLx5(i,j)+beta_k(5)*PLx5_P(i,j)+gamma_k(5)*(tmpE-Ex_P(i,j))

         Ex_P(i,j)=Ex(i,j)
         Ex(i,j)=tmpE

         PLx1_P(i,j)=PLx1(i,j)
         PLx1(i,j)=tmpPL1

         PLx2_P(i,j)=PLx2(i,j)
         PLx2(i,j)=tmpPL2

         PLx3_P(i,j)=PLx3(i,j)
         PLx3(i,j)=tmpPL3

         PLx4_P(i,j)=PLx4(i,j)
         PLx4(i,j)=tmpPL4

         PLx5_P(i,j)=PLx5(i,j)
         PLx5(i,j)=tmpPL5
        elseif(F_glass_x(i,j))then
          Ex(i,j)=Ex(i,j)+dt_eps_glass*(Hz(i,j)-tmp1)*den_ey(j)
        elseif(F_SU8_x(i,j))then
          Ex(i,j)=Ex(i,j)+dt_eps_SU8*(Hz(i,j)-tmp1)*den_ey(j)
        elseif(FAl2O3x(i,j))then
          Ex(i,j)=Ex(i,j)+dt_eps_Al2O3*(Hz(i,j)-tmp1)*den_ey(j)
        elseif(FPMMAx(i,j))then
          Ex(i,j)=Ex(i,j)+dt_eps_PMMA*(Hz(i,j)-tmp1)*den_ey(j)
        elseif(FMx(i,j))then
          Px_av=dPx(i,j)+dPx_old(i,j)+tmp2
          Ex(i,j)=Ex(i,j)+dt_epsM*(Hz(i,j)-tmp1)*den_ey(j)- &
     		                    cx*Px_av
        else
          Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-tmp1)*den_ey(j)
      endif
     enddo
   enddo

   !<~~~~ now apply PML at the bottom and at the top
   if(coords(2)==0)then                      !<~~~ lower y-boundary in air
     do i=1,Nx_loc
       do j=2,npml
         psi_Exy(i,j)=be_y_air(j)*psi_Exy(i,j)+ce_y_air(j)*(Hz(i,j)-Hz(i,j-1))/dy
         Ex(i,j)=Ex(i,j)+dt_eps_glass*psi_Exy(i,j)
       enddo
     enddo
   elseif(coords(2)==(dims(2)-1))then     !<~~~ upper y-boundary in air
     do i=1,Nx_loc
      jj=npml
       do j=Ny_loc+1-npml,Ny_loc-1
         psi_Exy(i,jj)=be_y_air(jj)*psi_Exy(i,jj)+ce_y_air(jj)*(Hz(i,j)-Hz(i,(j-1)))/dy
         Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy(i,jj)
         jj=jj-1
       enddo
     enddo
   endif

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !
    dPy_send(1:Nx_loc)=dPy(1:Nx_loc,1)+dPy_old(1:Nx_loc,1)
    sendtag5=500+n             !<~~~~ some tag integers
    recvtag5=sendtag5           !...
    call mpi_sendrecv(dPy_send,             &     !<=== sending
                      Nx_loc,               &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_down,             &     !<=== where sending
                      sendtag5,             &     !<=== sending tag
                      dPy_recv,             &     !<=== receiving
                      Nx_loc,               &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_up,               &     !<=== receiving from where
                      recvtag5,             &     !<=== receiving tag
                      cartesian_comm,       &     !<=== handle of Cartesian coordinates
                      istatus,              &     !<=== istatus
                      ierr)                       !<=== error code

    !<~~~== first we need to send Hz(Nx_loc,1:Ny_loc)-->right ==~~~
    Hzy_send(1:Ny_loc)=Hz(Nx_loc,1:Ny_loc)
    sendtag6=600+n             !<~~~~ some tag integers
    recvtag6=sendtag6           !...
    call mpi_sendrecv(Hzy_send,             &     !<=== sending
                      Ny_loc,               &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_right,            &     !<=== where sending
                      sendtag6,             &     !<=== sending tag
                      Hzy_recv,             &     !<=== receiving
                      Ny_loc,               &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_left,             &     !<=== receiving from where
                      recvtag6,             &     !<=== receiving tag
                      cartesian_comm,       &     !<=== handle of Cartesian coordinates
                      istatus,              &     !<=== istatus
                      ierr)                       !<=== error code

    do i=1,Nx_loc
      do j=1,Ny_loc
        if(i==1)then
            tmp1=Hzy_recv(j)
          else
            tmp1=Hz(i-1,j)
        endif

               if(j==Ny_loc)then
                   tmp2=dPy_recv(i)
                 else
                   tmp2=dPy(i,j+1)+dPy_old(i,j+1)
               endif

        if(FBy(i,j))then
          tmpE=C1*Ey(i,j)+C2*Ey_P(i,j)+C3*(tmp1-Hz(i,j))- &
               C4*PDy(i,j)- &
               (B1_k(1)*PLy1(i,j)+B2_k(1)*PLy1_P(i,j)+ &
                B1_k(2)*PLy2(i,j)+B2_k(2)*PLy2_P(i,j)+ &
                B1_k(3)*PLy3(i,j)+B2_k(3)*PLy3_P(i,j)+ &
                B1_k(4)*PLy4(i,j)+B2_k(4)*PLy4_P(i,j)+ &
                B1_k(5)*PLy5(i,j)+B2_k(5)*PLy5_P(i,j))

          PDy(i,j)=A1*PDy(i,j)+A2*(tmpE+Ey(i,j))

          tmpPL1=alpha_k(1)*PLy1(i,j)+beta_k(1)*PLy1_P(i,j)+gamma_k(1)*(tmpE-Ey_P(i,j))
          tmpPL2=alpha_k(2)*PLy2(i,j)+beta_k(2)*PLy2_P(i,j)+gamma_k(2)*(tmpE-Ey_P(i,j))
          tmpPL3=alpha_k(3)*PLy3(i,j)+beta_k(3)*PLy3_P(i,j)+gamma_k(3)*(tmpE-Ey_P(i,j))
          tmpPL4=alpha_k(4)*PLy4(i,j)+beta_k(4)*PLy4_P(i,j)+gamma_k(4)*(tmpE-Ey_P(i,j))
          tmpPL5=alpha_k(5)*PLy5(i,j)+beta_k(5)*PLy5_P(i,j)+gamma_k(5)*(tmpE-Ey_P(i,j))

          Ey_P(i,j)=Ey(i,j)
          Ey(i,j)=tmpE

          PLy1_P(i,j)=PLy1(i,j)
          PLy1(i,j)=tmpPL1

          PLy2_P(i,j)=PLy2(i,j)
          PLy2(i,j)=tmpPL2

          PLy3_P(i,j)=PLy3(i,j)
          PLy3(i,j)=tmpPL3

          PLy4_P(i,j)=PLy4(i,j)
          PLy4(i,j)=tmpPL4

          PLy5_P(i,j)=PLy5(i,j)
          PLy5(i,j)=tmpPL5
         elseif(F_glass_y(i,j))then
           Ey(i,j)=Ey(i,j)+dt_eps_glass*(tmp1-Hz(i,j))*den_ex
         elseif(F_SU8_y(i,j))then
           Ey(i,j)=Ey(i,j)+dt_eps_SU8*(tmp1-Hz(i,j))*den_ex
         elseif(FAl2O3y(i,j))then
           Ey(i,j)=Ey(i,j)+dt_eps_Al2O3*(tmp1-Hz(i,j))*den_ex
         elseif(FPMMAy(i,j))then
           Ey(i,j)=Ey(i,j)+dt_eps_PMMA*(tmp1-Hz(i,j))*den_ex
         elseif(FMy(i,j))then
           Py_av=dPy(i,j)+dPy_old(i,j)+tmp2
           Ey(i,j)=Ey(i,j)+dt_epsM*(tmp1-Hz(i,j))*den_ex- &
      		                    cy*Py_av
         else
           Ey(i,j)=Ey(i,j)+dt_eps0*(tmp1-Hz(i,j))*den_ex
       endif
      enddo
    enddo

    !~~@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@~~!
    					!~~~ Liouville-von Neumann equation with damping ~~~!
    !~~@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@~~!
    dPx_old=dPx
    dPy_old=dPy

    ExM_send(1:Ny_loc,1)=Ex_old(Nx_loc,1:Ny_loc)
    ExM_send(1:Ny_loc,2)=Ex(Nx_loc,1:Ny_loc)
    sendtag3=300+n              !<~~~~ some tag integers
    recvtag3=sendtag3           !...
    call mpi_sendrecv(ExM_send,             &     !<=== sending
                      2*Ny_loc,             &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_right,            &     !<=== where sending
                      sendtag3,             &     !<=== sending tag
                      ExM_recv,             &     !<=== receiving
                      2*Ny_loc,             &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_left,             &     !<=== receiving from where
                      recvtag3,             &     !<=== receiving tag
                      cartesian_comm,       &     !<=== handle of Cartesian coordinates
                      istatus,              &     !<=== istatus
                      ierr)                       !<=== error code

    EyM_send(1:Nx_loc,1)=Ey_old(1:Nx_loc,Ny_loc)
    EyM_send(1:Nx_loc,2)=Ey(1:Nx_loc,Ny_loc)
    sendtag5=500+n             !<~~~~ some tag integers
    recvtag5=sendtag5           !...
    call mpi_sendrecv(EyM_send,             &     !<=== sending
                      2*Nx_loc,             &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_up,               &     !<=== where sending
                      sendtag5,             &     !<=== sending tag
                      EyM_recv,             &     !<=== receiving
                      2*Nx_loc,             &     !<=== the size
                      mpi_double_precision, &     !<=== type
                      nbr_down,             &     !<=== receiving from where
                      recvtag5,             &     !<=== receiving tag
                      cartesian_comm,       &     !<=== handle of Cartesian coordinates
                      istatus,              &     !<=== istatus
                      ierr)                       !<=== error code

     do i=1,Nx_loc
      do j=1,Ny_loc
       if(FM(i,j))then

         if(i==1)then
                 avX_n=(ExM_recv(j,1)+Ex_old(i,j))/2.0
                 avX_n1=(ExM_recv(j,2)+Ex(i,j))/2.0
             else
                 avX_n=(Ex_old(i-1,j)+Ex_old(i,j))/2.0
                 avX_n1=(Ex(i-1,j)+Ex(i,j))/2.0
         endif

         if(j==1)then
               avY_n=(EyM_recv(i,1)+Ey_old(i,j))/2.0
               avY_n1=(EyM_recv(i,2)+Ey(i,j))/2.0
          else
    	         avY_n=(Ey_old(i,j-1)+Ey_old(i,j))/2.0
    	         avY_n1=(Ey(i,j-1)+Ey(i,j))/2.0
         endif

    	avX_n05=(avX_n1+avX_n)/2.0
    	avY_n05=(avY_n1+avY_n)/2.0

        Omega_plus_n=dp*(avX_n+Im*avY_n)/(h*sqrt6)
        Omega_minus_n=dp*(avX_n-Im*avY_n)/(h*sqrt6)

        Omega_plus_n05=dp*(avX_n05+Im*avY_n05)/(h*sqrt6)
        Omega_minus_n05=dp*(avX_n05-Im*avY_n05)/(h*sqrt6)

        Omega_plus_n1=dp*(avX_n1+Im*avY_n1)/(h*sqrt6)
        Omega_minus_n1=dp*(avX_n1-Im*avY_n1)/(h*sqrt6)

    !first RK step!
        rho(1,1)=rho11(i,j)
        rho(1,2)=rho12(i,j)
        rho(1,3)=rho13(i,j)

        rho(2,1)=conjg(rho12(i,j))
        rho(2,2)=rho22(i,j)
        rho(2,3)=rho23(i,j)

        rho(3,1)=conjg(rho13(i,j))
        rho(3,2)=conjg(rho23(i,j))
        rho(3,3)=rho33(i,j)

        Hamiltonian(1,2)=Omega_minus_n
        Hamiltonian(1,3)=-Omega_plus_n
        Hamiltonian(2,1)=conjg(Hamiltonian(1,2))
        Hamiltonian(3,1)=conjg(Hamiltonian(1,3))

        kk1=-Im*(matmul(Hamiltonian,rho)-matmul(rho,Hamiltonian))

        kk1(1,1)=kk1(1,1)+ &
                    gamma1*(rho(2,2)+rho(3,3))

        kk1(1,2)=kk1(1,2)- &
                    gamma2*rho(1,2)

        kk1(1,3)=kk1(1,3)- &
                    gamma2*rho(1,3)

        kk1(2,2)=kk1(2,2)- &
                    gamma1*rho(2,2)

        kk1(2,3)=kk1(2,3)- &
                    gamma2*rho(2,3)

        kk1(3,3)=kk1(3,3)- &
                    gamma1*rho(3,3)

        kk1=dt*kk1

    !second RK step!
        rho(1,1)=rho11(i,j)+kk1(1,1)/2.0
        rho(1,2)=rho12(i,j)+kk1(1,2)/2.0
        rho(1,3)=rho13(i,j)+kk1(1,3)/2.0

        rho(2,1)=conjg(rho12(i,j)+kk1(1,2)/2.0)
        rho(2,2)=rho22(i,j)+kk1(2,2)/2.0
        rho(2,3)=rho23(i,j)+kk1(2,3)/2.0

        rho(3,1)=conjg(rho13(i,j)+kk1(1,3)/2.0)
        rho(3,2)=conjg(rho23(i,j)+kk1(2,3)/2.0)
        rho(3,3)=rho33(i,j)+kk1(3,3)/2.0

        Hamiltonian(1,2)=Omega_minus_n05
        Hamiltonian(1,3)=-Omega_plus_n05
        Hamiltonian(2,1)=conjg(Hamiltonian(1,2))
        Hamiltonian(3,1)=conjg(Hamiltonian(1,3))

        kk2=-Im*(matmul(Hamiltonian,rho)-matmul(rho,Hamiltonian))

        kk2(1,1)=kk2(1,1)+ &
                    gamma1*(rho(2,2)+rho(3,3))

        kk2(1,2)=kk2(1,2)- &
                    gamma2*rho(1,2)

        kk2(1,3)=kk2(1,3)- &
                    gamma2*rho(1,3)

        kk2(2,2)=kk2(2,2)- &
                    gamma1*rho(2,2)

        kk2(2,3)=kk2(2,3)- &
                    gamma2*rho(2,3)

        kk2(3,3)=kk2(3,3)- &
                    gamma1*rho(3,3)

        kk2=dt*kk2

    !third RK step!
        rho(1,1)=rho11(i,j)+kk2(1,1)/2.0
        rho(1,2)=rho12(i,j)+kk2(1,2)/2.0
        rho(1,3)=rho13(i,j)+kk2(1,3)/2.0

        rho(2,1)=conjg(rho12(i,j)+kk2(1,2)/2.0)
        rho(2,2)=rho22(i,j)+kk2(2,2)/2.0
        rho(2,3)=rho23(i,j)+kk2(2,3)/2.0

        rho(3,1)=conjg(rho13(i,j)+kk2(1,3)/2.0)
        rho(3,2)=conjg(rho23(i,j)+kk2(2,3)/2.0)
        rho(3,3)=rho33(i,j)+kk2(3,3)/2.0

        Hamiltonian(1,2)=Omega_minus_n05
        Hamiltonian(1,3)=-Omega_plus_n05
        Hamiltonian(2,1)=conjg(Hamiltonian(1,2))
        Hamiltonian(3,1)=conjg(Hamiltonian(1,3))

        kk3=-Im*(matmul(Hamiltonian,rho)-matmul(rho,Hamiltonian))

        kk3(1,1)=kk3(1,1)+ &
                    gamma1*(rho(2,2)+rho(3,3))

        kk3(1,2)=kk3(1,2)- &
                    gamma2*rho(1,2)

        kk3(1,3)=kk3(1,3)- &
                    gamma2*rho(1,3)

        kk3(2,2)=kk3(2,2)- &
                    gamma1*rho(2,2)

        kk3(2,3)=kk3(2,3)- &
                    gamma2*rho(2,3)

        kk3(3,3)=kk3(3,3)- &
                    gamma1*rho(3,3)

        kk3=dt*kk3

    !fourth RK step!
        rho(1,1)=rho11(i,j)+kk3(1,1)
        rho(1,2)=rho12(i,j)+kk3(1,2)
        rho(1,3)=rho13(i,j)+kk3(1,3)

        rho(2,1)=conjg(rho12(i,j)+kk3(1,2))
        rho(2,2)=rho22(i,j)+kk3(2,2)
        rho(2,3)=rho23(i,j)+kk3(2,3)

        rho(3,1)=conjg(rho13(i,j)+kk3(1,3))
        rho(3,2)=conjg(rho23(i,j)+kk3(2,3))
        rho(3,3)=rho33(i,j)+kk3(3,3)

        Hamiltonian(1,2)=Omega_minus_n1
        Hamiltonian(1,3)=-Omega_plus_n1
        Hamiltonian(2,1)=conjg(Hamiltonian(1,2))
        Hamiltonian(3,1)=conjg(Hamiltonian(1,3))

        kk4=-Im*(matmul(Hamiltonian,rho)-matmul(rho,Hamiltonian))

        kk4(1,1)=kk4(1,1)+ &
                    gamma1*(rho(2,2)+rho(3,3))

        kk4(1,2)=kk4(1,2)- &
                    gamma2*rho(1,2)

        kk4(1,3)=kk4(1,3)- &
                    gamma2*rho(1,3)

        kk4(2,2)=kk4(2,2)- &
                    gamma1*rho(2,2)

        kk4(2,3)=kk4(2,3)- &
                    gamma2*rho(2,3)

        kk4(3,3)=kk4(3,3)- &
                    gamma1*rho(3,3)

        kk4=dt*kk4

    ! RK next step correction
        rho11(i,j)=rho11(i,j)+(kk1(1,1)+2.0*(kk2(1,1)+kk3(1,1))+kk4(1,1))/6.0
        rho12(i,j)=rho12(i,j)+(kk1(1,2)+2.0*(kk2(1,2)+kk3(1,2))+kk4(1,2))/6.0
        rho13(i,j)=rho13(i,j)+(kk1(1,3)+2.0*(kk2(1,3)+kk3(1,3))+kk4(1,3))/6.0
        rho22(i,j)=rho22(i,j)+(kk1(2,2)+2.0*(kk2(2,2)+kk3(2,2))+kk4(2,2))/6.0
        rho23(i,j)=rho23(i,j)+(kk1(2,3)+2.0*(kk2(2,3)+kk3(2,3))+kk4(2,3))/6.0
        rho33(i,j)=rho33(i,j)+(kk1(3,3)+2.0*(kk2(3,3)+kk3(3,3))+kk4(3,3))/6.0

    !polarization currents
        rho(1,1)=rho11(i,j)
        rho(1,2)=rho12(i,j)
        rho(1,3)=rho13(i,j)

        rho(2,1)=conjg(rho12(i,j))
        rho(2,2)=rho22(i,j)
        rho(2,3)=rho23(i,j)

        rho(3,1)=conjg(rho13(i,j))
        rho(3,2)=conjg(rho23(i,j))
        rho(3,3)=rho33(i,j)

        tmp_P=(0.0,0.0)
        tmp_P=matmul(rho,mux)
        Px_loc=n0*dble(tmp_P(1,1)+tmp_P(2,2)+tmp_P(3,3))

        tmp_P=(0.0,0.0)
        tmp_P=matmul(rho,muy)
        Py_loc=n0*dble(tmp_P(1,1)+tmp_P(2,2)+tmp_P(3,3))

        dPx(i,j)=Px_loc-Px_old(i,j)
        dPy(i,j)=Py_loc-Py_old(i,j)

        Px_old(i,j)=Px_loc
        Py_old(i,j)=Py_loc

        endif
      enddo
     enddo



  if (mod(n,Nt_skip)==0) then
      n_new=int(n/Nt_skip)
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !:::::::::::::::::::::::::::::::::::::::! detection !:::::::::::::::::::::::::::::::::::::::!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    if(coords(2)==mwR)then
      j=jwR       !<=== make it always >1 to avoid send/receive up/down
      !<~~~== first we need to send Ex(Nx_loc,j),Hz(Nx_loc,j)+Hz(Nx_loc,j-1)-->right ==~~~
      temp_send(1)=Ex(Nx_loc,j)
      temp_send(2)=Hz(Nx_loc,j)+Hz(Nx_loc,j-1)
      sendtag7=700+n              !<~~~~ some tag integers
      recvtag7=sendtag7           !...
      call mpi_sendrecv(temp_send,            &     !<=== sending
                        2,                    &     !<=== the size
                        mpi_double_precision, &     !<=== type
                        nbr_right,            &     !<=== where sending
                        sendtag7,             &     !<=== sending tag
                        temp_recv,            &     !<=== receiving
                        2,                    &     !<=== the size
                        mpi_double_precision, &     !<=== type
                        nbr_left,             &     !<=== receiving from where
                        recvtag7,             &     !<=== receiving tag
                        cartesian_comm,       &     !<=== handle of Cartesian coordinates
                        istatus,              &     !<=== istatus
                        ierr)                       !<=== error code

       do i=1,Nx_loc
         if(i==1)then
            Ex_temp(i,n_new)=(temp_recv(1)+Ex(i,j))/2.0
            Hz_temp(i,n_new)=(temp_recv(2)+Hz(i,j)+Hz(i,j-1))/4.0
           else
            Ex_temp(i,n_new)=(Ex(i-1,j)+Ex(i,j))/2.0
            Hz_temp(i,n_new)=(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0
          endif
       enddo
    elseif(coords(2)==mwT)then
      j=jwT       !<=== make it always >1
      !<~~~== first we need to send Ex(Nx_loc,j),Hz(Nx_loc,j)+Hz(Nx_loc,j-1)-->right ==~~~
      temp_send(1)=Ex(Nx_loc,j)
      temp_send(2)=Hz(Nx_loc,j)+Hz(Nx_loc,j-1)
      sendtag8=800+n              !<~~~~ some tag integers
      recvtag8=sendtag8           !...
      call mpi_sendrecv(temp_send,            &     !<=== sending
                        2,                    &     !<=== the size
                        mpi_double_precision, &     !<=== type
                        nbr_right,            &     !<=== where sending
                        sendtag8,             &     !<=== sending tag
                        temp_recv,            &     !<=== receiving
                        2,                    &     !<=== the size
                        mpi_double_precision, &     !<=== type
                        nbr_left,             &     !<=== receiving from where
                        recvtag8,             &     !<=== receiving tag
                        cartesian_comm,       &     !<=== handle of Cartesian coordinates
                        istatus,              &     !<=== istatus
                        ierr)                       !<=== error code

       do i=1,Nx_loc
         if(i==1)then
            Ex_temp(i,n_new)=(temp_recv(1)+Ex(i,j))/2.0
            Hz_temp(i,n_new)=(temp_recv(2)+Hz(i,j)+Hz(i,j-1))/4.0
           else
            Ex_temp(i,n_new)=(Ex(i-1,j)+Ex(i,j))/2.0
            Hz_temp(i,n_new)=(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0
          endif

          Ex_temp_inc(i,n_new)=Ex_inc(j)
          Hz_temp_inc(i,n_new)=(Hz_inc(j-1)+Hz_inc(j))/2.0
       enddo
    endif

    ! -------------------------------------------------------------------------!
    !                           population observation 
    ! -------------------------------------------------------------------------!
    excited_population = 0.0d0
    population(n/Nt_skip) = 0.0d0
    do i=1,Nx_loc
      do j=1,Ny_loc
        if(FM(i,j))then
          excited_population = excited_population + abs(rho22(i,j)) + abs(rho33(i,j))
          population(n/Nt_skip) = population(n/Nt_skip) + abs(rho22(i,j)) + abs(rho33(i,j))
        endif
      enddo
    enddo
    ! write(*,*) n/Nt_skip,myrank,population(n/Nt_skip)

    call MPI_REDUCE(excited_population,         &
                    reduced_excited_population, &
                    1,                          &
                    MPI_DOUBLE_PRECISION,       &
                    MPI_SUM,                    &
                    mw_reduce,                  &
                    MPI_COMM_WORLD,             &
                    ierr)
    if(myrank==mw_reduce)then
      write(*,*) dt*n*1.0d15,reduced_excited_population
    endif
 endif

  ! if(mod(n,5000)==0.and.myrank==0)then
  !   call cpu_time(cpu2)
  !   write(*,*) 'cpu time [minutes]',(cpu2-cpu1)*(Nt/5000)/60.0
  !   cpu1=cpu2
  !   write(*,*) n
  !   write(*,*) Ex(Nx_loc/2,Ny_loc/2),Ey(Nx_loc/2,Ny_loc/2)
  !   write(*,*) Ex_inc(Ny_loc/2),Hz_inc(Ny_loc/2)
  ! endif

  if (mod(n-n_print,Nt_skip)==0 .and. n > n_print) then
    write(filehead, "(A5,I5.5)") "Ex_t=", INT(n-n_print)
    write(*,"(A10)") filehead
    call WriteMatrixToFile(myrank,filehead, Ex, Nx_loc, Ny_loc)
        
    write(filehead, "(A5,I5.5)") "Ey_t=", INT(n-n_print)
    call WriteMatrixToFile(myrank,filehead, Ey, Nx_loc, Ny_loc)
    
    write(filehead, "(A5,I5.5)") "Hz_t=", INT(n-n_print)
    call WriteMatrixToFile(myrank,filehead, Hz, Nx_loc, Ny_loc) 
  endif

 enddo!<-- end of time iterations

!~~~~=== post-processing ===~~~~!
if(coords(2)==mwR)then
  call zffti(Nt_new,wsave)

  do i=1,Nx_loc
    call zfftf(Nt_new,Ex_temp(i,:),wsave)
    call zfftf(Nt_new,Hz_temp(i,:),wsave)
  enddo

  Ex_temp=Ex_temp/sqrt(1.0d0*Nt_new)
  Hz_temp=Hz_temp/sqrt(1.0d0*Nt_new)

  nn=0
  do n=2,Nt_new/2
      tmp=dble(n-1)/(dt*Nt_skip*Nt_new)
      tmp=Hz_to_ev*tmp
      if ((tmp>=omega_min).and.(tmp<=omega_max)) then
          nn=nn+1
          sum1=(0.0D0,0.0D0)
          do i=1,Nx_loc
            cTMP1=Ex_temp(i,n)
            cTMP2=Hz_temp(i,n)
            sum1=sum1+cTMP1*conjg(cTMP2)*Nt_skip
          enddo
          R_loc(nn)=dreal(sum1)/(Nx_loc)
          R_loc(nn)=abs(R_loc(nn))
      endif
  enddo
endif
call mpi_reduce(R_loc,R_global,resolution_fr,mpi_double_precision,mpi_sum,mw_reduce,MPI_COMM_WORLD,ierr)


if(coords(2)==mwT)then
  call zffti(Nt_new,wsave)

  do i=1,Nx_loc
    call zfftf(Nt_new,Ex_temp(i,:),wsave)
    call zfftf(Nt_new,Hz_temp(i,:),wsave)
  enddo

  Ex_temp=Ex_temp/sqrt(1.0d0*Nt_new)
  Hz_temp=Hz_temp/sqrt(1.0d0*Nt_new)

  nn=0
  do n=2,Nt_new/2
      tmp=dble(n-1)/(dt*Nt_skip*Nt_new)
      tmp=Hz_to_ev*tmp
      if ((tmp>=omega_min).and.(tmp<=omega_max)) then
          nn=nn+1
          sum1=(0.0D0,0.0D0)
          do i=1,Nx_loc
            cTMP1=Ex_temp(i,n)
            cTMP2=Hz_temp(i,n)
            sum1=sum1+cTMP1*conjg(cTMP2)*Nt_skip
          enddo
          T_loc(nn)=dreal(sum1)/(Nx_loc)
          T_loc(nn)=abs(T_loc(nn))
      endif
  enddo

  do i=1,Nx_loc
    call zfftf(Nt_new,Ex_temp_inc(i,:),wsave)
    call zfftf(Nt_new,Hz_temp_inc(i,:),wsave)
  enddo

  Ex_temp_inc=Ex_temp_inc/sqrt(1.0d0*Nt_new)
  Hz_temp_inc=Hz_temp_inc/sqrt(1.0d0*Nt_new)

  nn=0
  do n=2,Nt_new/2
      tmp=dble(n-1)/(dt*Nt_skip*Nt_new)
      tmp=Hz_to_ev*tmp
      if ((tmp>=omega_min).and.(tmp<=omega_max)) then
          nn=nn+1
          sum1=(0.0D0,0.0D0)
          do i=1,Nx_loc
            cTMP1=Ex_temp_inc(i,n)
            cTMP2=Hz_temp_inc(i,n)
            sum1=sum1+cTMP1*conjg(cTMP2)*Nt_skip
          enddo
          P_inc_loc(nn)=dreal(sum1)/(Nx_loc)
          P_inc_loc(nn)=abs(P_inc_loc(nn))
      endif
  enddo
endif

call mpi_reduce(T_loc,T_global,resolution_fr,mpi_double_precision,mpi_sum,mw_reduce,MPI_COMM_WORLD,ierr)
call mpi_reduce(P_inc_loc,P_inc_global,resolution_fr,mpi_double_precision,mpi_sum,mw_reduce,MPI_COMM_WORLD,ierr)


if(myrank==mw_reduce)then
  !  open(file='T_R-DL-mol_1p578eV_1e26-cylinders_delta5nm-amplitude_50nm-glass1p4_ONLY_Au50nm-Al2O3_20nm-PMMA.dat',unit=32)
  open(file='spectrum.out',unit=32)
   do nn=1,N_w
    write(32,*) omega_P(nn)/ev_to_radsec,T_global(nn)/P_inc_global(nn),R_global(nn)/P_inc_global(nn)
   enddo
   close(unit=32)
endif

!! print out cart coordinates for each {myrank} processor
if ( myrank == mw_reduce) then
  open(file='./result/coords.dat', unit=323)
  write(323, "(A)", advance="no") ""
  close(unit=323)
  open(file='./result/coords.dat', unit=323, Access="append")
  do  i=0, nprocs-1   
    call mpi_cart_coords( &
          cartesian_comm, &
          i,              &
          ndims,          &
          io_coords,      &     !<--- coordinates(ndims) = gives integers identifying local coordinates of a given block for myrank
          ierr            &     !<--- coords(1) = x, coords(2) = y
                       )

    write(323, *) i, io_coords(1), io_coords(2)
  enddo
  close(unit=323)
endif
!! print out fields
write(filename, "(A15, I2.2, A4)") "./result/metal_", myrank, ".dat"
open(file=filename, unit=323)
!write(323, *) FB(:, :)
! Print the 2D array
do j = 1, Ny_loc
  do i = 1, Nx_loc
    if( FB(i,j)) then
      write(323, "(I2)", advance='no') 1
    else
      write(323, "(I2)", advance='no') 0
    end if
  end do
  write(323, *)  ! Newline after each row
end do

close(unit=323)

call WriteMatrixToFile(myrank,"Ex___final", Ex, Nx_loc, Ny_loc)
call WriteMatrixToFile(myrank,"Ey___final", Ey, Nx_loc, Ny_loc)
call WriteMatrixToFile(myrank,"Hz___final", Hz, Nx_loc, Ny_loc)




! ! Output total population 
! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
! call MPI_GATHER(population, Nt/Nt_skip, MPI_DOUBLE_PRECISION, total_population, Nt/Nt_skip, MPI_DOUBLE_PRECISION, mw_reduce, MPI_COMM_WORLD, ierr)

! if(myrank==mw_reduce)then
!   do i = 1, (nprocs-1)
!     do n = 1, Nt/Nt_skip
!       total_population(n) = total_population(n) + total_population(i * Nt/Nt_skip + n)
!     enddo
!   enddo
  
!   open(file='population.dat',unit=33)
!   do n=1,Nt/Nt_skip
!     write(33,*) dt*n*1.0d15*Nt_skip,total_population(n)
!   enddo
!   close(unit=33)
! endif
!============================================================================!
!----------------------------------------------------------------------------!
 call MPI_FINALIZE(ierr)
!----------------------------------------------------------------------------!
!============================================================================!
end
subroutine WriteMatrixToFile(myrank, filehead, Ey, Nx_loc, Ny_loc)
  implicit none

  ! Input arguments
  integer, intent(in) :: myrank
  integer, intent(in) :: Nx_loc, Ny_loc
  double precision, intent(in) :: Ey(Nx_loc, Ny_loc)
  character(len=10), intent(in) :: filehead

  ! Local variables
  character(len=34) :: filename  ! length much match exactly to input str len
  integer :: i, j

  ! Construct the filename
  write(filename, "(A9,A10,A1,I2.2, A4)") "./result/", &
                filehead,"_", myrank, ".dat"
  write(*,*) filename
  ! Open the file for writing
  open(file=filename, unit=328)

  ! Write the 2D array to the file
  ! Print the 2D array
  do j = 1, Ny_loc
    do i = 1, Nx_loc
      write(328, "(ES20.8)", advance='no') Ey(i, j)
      write(328, "(A)", advance='no') " "
    end do
    write(328, *)  ! Newline after each row
  end do

  ! Close the file
  close(328)

end subroutine WriteMatrixToFile

SUBROUTINE ZFFTI (N,WSAVE)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       WSAVE(1)
IF (N .EQ. 1) RETURN
IW1 = N+N+1
IW2 = IW1+N+N
CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
RETURN
END

SUBROUTINE ZFFTF (N,C,WSAVE)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       C(1)       ,WSAVE(1)
IF (N .EQ. 1) RETURN
IW1 = N+N+1
IW2 = IW1+N+N
CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
RETURN
END

SUBROUTINE ZFFTB (N,C,WSAVE)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       C(1)       ,WSAVE(1)
IF (N .EQ. 1) RETURN
IW1 = N+N+1
IW2 = IW1+N+N
CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
RETURN
END

SUBROUTINE CFFTI1 (N,WA,IFAC)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
NL = N
NF = 0
J = 0
101 J = J+1
IF (J-4) 102,102,103
102 NTRY = NTRYH(J)
GO TO 104
103 NTRY = NTRY+2
104 NQ = NL/NTRY
NR = NL-NTRY*NQ
IF (NR) 101,105,101
105 NF = NF+1
IFAC(NF+2) = NTRY
NL = NQ
IF (NTRY .NE. 2) GO TO 107
IF (NF .EQ. 1) GO TO 107
DO 106 I=2,NF
   IB = NF-I+2
   IFAC(IB+2) = IFAC(IB+1)
106 CONTINUE
IFAC(3) = 2
107 IF (NL .NE. 1) GO TO 104
IFAC(1) = N
IFAC(2) = NF
TPI =  6.28318530717958647692D0
ARGH = TPI/FLOAT(N)
I = 2
L1 = 1
DO 110 K1=1,NF
   IP = IFAC(K1+2)
   LD = 0
   L2 = L1*IP
   IDO = N/L2
   IDOT = IDO+IDO+2
   IPM = IP-1
   DO 109 J=1,IPM
      I1 = I
      WA(I-1) = 1.0D0
      WA(I) = 0.0D0
      LD = LD+L1
      FI = 0.0D0
      ARGLD = FLOAT(LD)*ARGH
      DO 108 II=4,IDOT,2
         I = I+2
         FI = FI+1.D0
         ARG = FI*ARGLD
         WA(I-1) = COS(ARG)
         WA(I) = SIN(ARG)
108       CONTINUE
      IF (IP .LE. 5) GO TO 109
      WA(I1-1) = WA(I-1)
      WA(I1) = WA(I)
109    CONTINUE
   L1 = L2
110 CONTINUE
RETURN
END

SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
NF = IFAC(2)
NA = 0
L1 = 1
IW = 1
DO 116 K1=1,NF
   IP = IFAC(K1+2)
   L2 = IP*L1
   IDO = N/L2
   IDOT = IDO+IDO
   IDL1 = IDOT*L1
   IF (IP .NE. 4) GO TO 103
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IF (NA .NE. 0) GO TO 101
   CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
   GO TO 102
101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
102    NA = 1-NA
   GO TO 115
103    IF (IP .NE. 2) GO TO 106
   IF (NA .NE. 0) GO TO 104
   CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
   GO TO 105
104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
105    NA = 1-NA
   GO TO 115
106    IF (IP .NE. 3) GO TO 109
   IX2 = IW+IDOT
   IF (NA .NE. 0) GO TO 107
   CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
   GO TO 108
107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
108    NA = 1-NA
   GO TO 115
109    IF (IP .NE. 5) GO TO 112
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IX4 = IX3+IDOT
   IF (NA .NE. 0) GO TO 110
   CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
   GO TO 111
110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
111    NA = 1-NA
   GO TO 115
112    IF (NA .NE. 0) GO TO 113
   CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
   GO TO 114
113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
114    IF (NAC .NE. 0) NA = 1-NA
115    L1 = L2
   IW = IW+(IP-1)*IDOT
116 CONTINUE
IF (NA .EQ. 0) RETURN
N2 = N+N
DO 117 I=1,N2
   C(I) = CH(I)
117 CONTINUE
RETURN
END

SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
NF = IFAC(2)
NA = 0
L1 = 1
IW = 1
DO 116 K1=1,NF
   IP = IFAC(K1+2)
   L2 = IP*L1
   IDO = N/L2
   IDOT = IDO+IDO
   IDL1 = IDOT*L1
   IF (IP .NE. 4) GO TO 103
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IF (NA .NE. 0) GO TO 101
   CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
   GO TO 102
101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
102    NA = 1-NA
   GO TO 115
103    IF (IP .NE. 2) GO TO 106
   IF (NA .NE. 0) GO TO 104
   CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
   GO TO 105
104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
105    NA = 1-NA
   GO TO 115
106    IF (IP .NE. 3) GO TO 109
   IX2 = IW+IDOT
   IF (NA .NE. 0) GO TO 107
   CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
   GO TO 108
107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
108    NA = 1-NA
   GO TO 115
109    IF (IP .NE. 5) GO TO 112
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IX4 = IX3+IDOT
   IF (NA .NE. 0) GO TO 110
   CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
   GO TO 111
110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
111    NA = 1-NA
   GO TO 115
112    IF (NA .NE. 0) GO TO 113
   CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
   GO TO 114
113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
114    IF (NAC .NE. 0) NA = 1-NA
115    L1 = L2
   IW = IW+(IP-1)*IDOT
116 CONTINUE
IF (NA .EQ. 0) RETURN
N2 = N+N
DO 117 I=1,N2
   C(I) = CH(I)
117 CONTINUE
RETURN
END

SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                CH2(IDL1,IP)
IDOT = IDO/2
NT = IP*IDL1
IPP2 = IP+2
IPPH = (IP+1)/2
IDP = IP*IDO
IF (IDO .LT. L1) GO TO 106
DO 103 J=2,IPPH
   JC = IPP2-J
   DO 102 K=1,L1
      DO 101 I=1,IDO
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
101       CONTINUE
102    CONTINUE
103 CONTINUE
DO 105 K=1,L1
   DO 104 I=1,IDO
      CH(I,K,1) = CC(I,1,K)
104    CONTINUE
105 CONTINUE
GO TO 112
106 DO 109 J=2,IPPH
   JC = IPP2-J
   DO 108 I=1,IDO
      DO 107 K=1,L1
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
107       CONTINUE
108    CONTINUE
109 CONTINUE
DO 111 I=1,IDO
   DO 110 K=1,L1
      CH(I,K,1) = CC(I,1,K)
110    CONTINUE
111 CONTINUE
112 IDL = 2-IDO
INC = 0
DO 116 L=2,IPPH
   LC = IPP2-L
   IDL = IDL+IDO
   DO 113 IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
      C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
113    CONTINUE
   IDLJ = IDL
   INC = INC+IDO
   DO 115 J=3,IPPH
      JC = IPP2-J
      IDLJ = IDLJ+INC
      IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
      WAR = WA(IDLJ-1)
      WAI = WA(IDLJ)
      DO 114 IK=1,IDL1
         C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
         C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
114       CONTINUE
115    CONTINUE
116 CONTINUE
DO 118 J=2,IPPH
   DO 117 IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
117    CONTINUE
118 CONTINUE
DO 120 J=2,IPPH
   JC = IPP2-J
   DO 119 IK=2,IDL1,2
      CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
      CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
      CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
      CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
119    CONTINUE
120 CONTINUE
NAC = 1
IF (IDO .EQ. 2) RETURN
NAC = 0
DO 121 IK=1,IDL1
   C2(IK,1) = CH2(IK,1)
121 CONTINUE
DO 123 J=2,IP
   DO 122 K=1,L1
      C1(1,K,J) = CH(1,K,J)
      C1(2,K,J) = CH(2,K,J)
122    CONTINUE
123 CONTINUE
IF (IDOT .GT. L1) GO TO 127
IDIJ = 0
DO 126 J=2,IP
   IDIJ = IDIJ+2
   DO 125 I=4,IDO,2
      IDIJ = IDIJ+2
      DO 124 K=1,L1
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
124       CONTINUE
125    CONTINUE
126 CONTINUE
RETURN
127 IDJ = 2-IDO
DO 130 J=2,IP
   IDJ = IDJ+IDO
   DO 129 K=1,L1
      IDIJ = IDJ
      DO 128 I=4,IDO,2
         IDIJ = IDIJ+2
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
128       CONTINUE
129    CONTINUE
130 CONTINUE
RETURN
END

SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
                WA1(1)
IF (IDO .GT. 2) GO TO 102
DO 101 K=1,L1
   CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
   CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
   CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
   CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
      TR2 = CC(I-1,1,K)-CC(I-1,2,K)
      CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
      TI2 = CC(I,1,K)-CC(I,2,K)
      CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
      CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
                WA1(1)     ,WA2(1)
!     *** TAUI IS -SQRT(3)/2 ***
DATA TAUR,TAUI /-0.5D0,-0.86602540378443864676D0/
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TR2 = CC(1,2,K)+CC(1,3,K)
   CR2 = CC(1,1,K)+TAUR*TR2
   CH(1,K,1) = CC(1,1,K)+TR2
   TI2 = CC(2,2,K)+CC(2,3,K)
   CI2 = CC(2,1,K)+TAUR*TI2
   CH(2,K,1) = CC(2,1,K)+TI2
   CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
   CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
   CH(1,K,2) = CR2-CI3
   CH(1,K,3) = CR2+CI3
   CH(2,K,2) = CI2+CR3
   CH(2,K,3) = CI2-CR3
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TR2 = CC(I-1,2,K)+CC(I-1,3,K)
      CR2 = CC(I-1,1,K)+TAUR*TR2
      CH(I-1,K,1) = CC(I-1,1,K)+TR2
      TI2 = CC(I,2,K)+CC(I,3,K)
      CI2 = CC(I,1,K)+TAUR*TI2
      CH(I,K,1) = CC(I,1,K)+TI2
      CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
      CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
      DR2 = CR2-CI3
      DR3 = CR2+CI3
      DI2 = CI2+CR3
      DI3 = CI2-CR3
      CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
      CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
      CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
      CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
                WA1(1)     ,WA2(1)     ,WA3(1)
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TI1 = CC(2,1,K)-CC(2,3,K)
   TI2 = CC(2,1,K)+CC(2,3,K)
   TR4 = CC(2,2,K)-CC(2,4,K)
   TI3 = CC(2,2,K)+CC(2,4,K)
   TR1 = CC(1,1,K)-CC(1,3,K)
   TR2 = CC(1,1,K)+CC(1,3,K)
   TI4 = CC(1,4,K)-CC(1,2,K)
   TR3 = CC(1,2,K)+CC(1,4,K)
   CH(1,K,1) = TR2+TR3
   CH(1,K,3) = TR2-TR3
   CH(2,K,1) = TI2+TI3
   CH(2,K,3) = TI2-TI3
   CH(1,K,2) = TR1+TR4
   CH(1,K,4) = TR1-TR4
   CH(2,K,2) = TI1+TI4
   CH(2,K,4) = TI1-TI4
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI1 = CC(I,1,K)-CC(I,3,K)
      TI2 = CC(I,1,K)+CC(I,3,K)
      TI3 = CC(I,2,K)+CC(I,4,K)
      TR4 = CC(I,2,K)-CC(I,4,K)
      TR1 = CC(I-1,1,K)-CC(I-1,3,K)
      TR2 = CC(I-1,1,K)+CC(I-1,3,K)
      TI4 = CC(I-1,4,K)-CC(I-1,2,K)
      TR3 = CC(I-1,2,K)+CC(I-1,4,K)
      CH(I-1,K,1) = TR2+TR3
      CR3 = TR2-TR3
      CH(I,K,1) = TI2+TI3
      CI3 = TI2-TI3
      CR2 = TR1+TR4
      CR4 = TR1-TR4
      CI2 = TI1+TI4
      CI4 = TI1-TI4
      CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
      CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
      CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
      CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
      CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
      CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
!     *** TR11=COS(2*PI/5), TI11=-SIN(2*PI/5)
!     *** TR12=-COS(4*PI/5), TI12=-SIN(4*PI/5)
DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0, &
     -0.95105651629515357212D0, &
     -0.8090169943749474241D0, -0.58778525229247312917D0/
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TI5 = CC(2,2,K)-CC(2,5,K)
   TI2 = CC(2,2,K)+CC(2,5,K)
   TI4 = CC(2,3,K)-CC(2,4,K)
   TI3 = CC(2,3,K)+CC(2,4,K)
   TR5 = CC(1,2,K)-CC(1,5,K)
   TR2 = CC(1,2,K)+CC(1,5,K)
   TR4 = CC(1,3,K)-CC(1,4,K)
   TR3 = CC(1,3,K)+CC(1,4,K)
   CH(1,K,1) = CC(1,1,K)+TR2+TR3
   CH(2,K,1) = CC(2,1,K)+TI2+TI3
   CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
   CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
   CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
   CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
   CR5 = TI11*TR5+TI12*TR4
   CI5 = TI11*TI5+TI12*TI4
   CR4 = TI12*TR5-TI11*TR4
   CI4 = TI12*TI5-TI11*TI4
   CH(1,K,2) = CR2-CI5
   CH(1,K,5) = CR2+CI5
   CH(2,K,2) = CI2+CR5
   CH(2,K,3) = CI3+CR4
   CH(1,K,3) = CR3-CI4
   CH(1,K,4) = CR3+CI4
   CH(2,K,4) = CI3-CR4
   CH(2,K,5) = CI2-CR5
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI5 = CC(I,2,K)-CC(I,5,K)
      TI2 = CC(I,2,K)+CC(I,5,K)
      TI4 = CC(I,3,K)-CC(I,4,K)
      TI3 = CC(I,3,K)+CC(I,4,K)
      TR5 = CC(I-1,2,K)-CC(I-1,5,K)
      TR2 = CC(I-1,2,K)+CC(I-1,5,K)
      TR4 = CC(I-1,3,K)-CC(I-1,4,K)
      TR3 = CC(I-1,3,K)+CC(I-1,4,K)
      CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
      CH(I,K,1) = CC(I,1,K)+TI2+TI3
      CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
      CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
      CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
      CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
      CR5 = TI11*TR5+TI12*TR4
      CI5 = TI11*TI5+TI12*TI4
      CR4 = TI12*TR5-TI11*TR4
      CI4 = TI12*TI5-TI11*TI4
      DR3 = CR3-CI4
      DR4 = CR3+CI4
      DI3 = CI3+CR4
      DI4 = CI3-CR4
      DR5 = CR2+CI5
      DR2 = CR2-CI5
      DI5 = CI2-CR5
      DI2 = CI2+CR5
      CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
      CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
      CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
      CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
      CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
      CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
      CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
      CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
                C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                CH2(IDL1,IP)
IDOT = IDO/2
NT = IP*IDL1
IPP2 = IP+2
IPPH = (IP+1)/2
IDP = IP*IDO
IF (IDO .LT. L1) GO TO 106
DO 103 J=2,IPPH
   JC = IPP2-J
   DO 102 K=1,L1
      DO 101 I=1,IDO
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
101       CONTINUE
102    CONTINUE
103 CONTINUE
DO 105 K=1,L1
   DO 104 I=1,IDO
      CH(I,K,1) = CC(I,1,K)
104    CONTINUE
105 CONTINUE
GO TO 112
106 DO 109 J=2,IPPH
   JC = IPP2-J
   DO 108 I=1,IDO
      DO 107 K=1,L1
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
107       CONTINUE
108    CONTINUE
109 CONTINUE
DO 111 I=1,IDO
   DO 110 K=1,L1
      CH(I,K,1) = CC(I,1,K)
110    CONTINUE
111 CONTINUE
112 IDL = 2-IDO
INC = 0
DO 116 L=2,IPPH
   LC = IPP2-L
   IDL = IDL+IDO
   DO 113 IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
      C2(IK,LC) = WA(IDL)*CH2(IK,IP)
113    CONTINUE
   IDLJ = IDL
   INC = INC+IDO
   DO 115 J=3,IPPH
      JC = IPP2-J
      IDLJ = IDLJ+INC
      IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
      WAR = WA(IDLJ-1)
      WAI = WA(IDLJ)
      DO 114 IK=1,IDL1
         C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
         C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
114       CONTINUE
115    CONTINUE
116 CONTINUE
DO 118 J=2,IPPH
   DO 117 IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
117    CONTINUE
118 CONTINUE
DO 120 J=2,IPPH
   JC = IPP2-J
   DO 119 IK=2,IDL1,2
      CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
      CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
      CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
      CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
119    CONTINUE
120 CONTINUE
NAC = 1
IF (IDO .EQ. 2) RETURN
NAC = 0
DO 121 IK=1,IDL1
   C2(IK,1) = CH2(IK,1)
121 CONTINUE
DO 123 J=2,IP
   DO 122 K=1,L1
      C1(1,K,J) = CH(1,K,J)
      C1(2,K,J) = CH(2,K,J)
122    CONTINUE
123 CONTINUE
IF (IDOT .GT. L1) GO TO 127
IDIJ = 0
DO 126 J=2,IP
   IDIJ = IDIJ+2
   DO 125 I=4,IDO,2
      IDIJ = IDIJ+2
      DO 124 K=1,L1
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
124       CONTINUE
125    CONTINUE
126 CONTINUE
RETURN
127 IDJ = 2-IDO
DO 130 J=2,IP
   IDJ = IDJ+IDO
   DO 129 K=1,L1
      IDIJ = IDJ
      DO 128 I=4,IDO,2
         IDIJ = IDIJ+2
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
128       CONTINUE
129    CONTINUE
130 CONTINUE
RETURN
END

SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &
                WA1(1)
IF (IDO .GT. 2) GO TO 102
DO 101 K=1,L1
   CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
   CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
   CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
   CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
      TR2 = CC(I-1,1,K)-CC(I-1,2,K)
      CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
      TI2 = CC(I,1,K)-CC(I,2,K)
      CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
      CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &
                WA1(1)     ,WA2(1)
!     *** TAUI IS SQRT(3)/2 ***
DATA TAUR,TAUI /-0.5D0,0.86602540378443864676D0/
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TR2 = CC(1,2,K)+CC(1,3,K)
   CR2 = CC(1,1,K)+TAUR*TR2
   CH(1,K,1) = CC(1,1,K)+TR2
   TI2 = CC(2,2,K)+CC(2,3,K)
   CI2 = CC(2,1,K)+TAUR*TI2
   CH(2,K,1) = CC(2,1,K)+TI2
   CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
   CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
   CH(1,K,2) = CR2-CI3
   CH(1,K,3) = CR2+CI3
   CH(2,K,2) = CI2+CR3
   CH(2,K,3) = CI2-CR3
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TR2 = CC(I-1,2,K)+CC(I-1,3,K)
      CR2 = CC(I-1,1,K)+TAUR*TR2
      CH(I-1,K,1) = CC(I-1,1,K)+TR2
      TI2 = CC(I,2,K)+CC(I,3,K)
      CI2 = CC(I,1,K)+TAUR*TI2
      CH(I,K,1) = CC(I,1,K)+TI2
      CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
      CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
      DR2 = CR2-CI3
      DR3 = CR2+CI3
      DI2 = CI2+CR3
      DI3 = CI2-CR3
      CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
      CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
      CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
      CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &
                WA1(1)     ,WA2(1)     ,WA3(1)
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TI1 = CC(2,1,K)-CC(2,3,K)
   TI2 = CC(2,1,K)+CC(2,3,K)
   TR4 = CC(2,4,K)-CC(2,2,K)
   TI3 = CC(2,2,K)+CC(2,4,K)
   TR1 = CC(1,1,K)-CC(1,3,K)
   TR2 = CC(1,1,K)+CC(1,3,K)
   TI4 = CC(1,2,K)-CC(1,4,K)
   TR3 = CC(1,2,K)+CC(1,4,K)
   CH(1,K,1) = TR2+TR3
   CH(1,K,3) = TR2-TR3
   CH(2,K,1) = TI2+TI3
   CH(2,K,3) = TI2-TI3
   CH(1,K,2) = TR1+TR4
   CH(1,K,4) = TR1-TR4
   CH(2,K,2) = TI1+TI4
   CH(2,K,4) = TI1-TI4
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI1 = CC(I,1,K)-CC(I,3,K)
      TI2 = CC(I,1,K)+CC(I,3,K)
      TI3 = CC(I,2,K)+CC(I,4,K)
      TR4 = CC(I,4,K)-CC(I,2,K)
      TR1 = CC(I-1,1,K)-CC(I-1,3,K)
      TR2 = CC(I-1,1,K)+CC(I-1,3,K)
      TI4 = CC(I-1,2,K)-CC(I-1,4,K)
      TR3 = CC(I-1,2,K)+CC(I-1,4,K)
      CH(I-1,K,1) = TR2+TR3
      CR3 = TR2-TR3
      CH(I,K,1) = TI2+TI3
      CI3 = TI2-TI3
      CR2 = TR1+TR4
      CR4 = TR1-TR4
      CI2 = TI1+TI4
      CI4 = TI1-TI4
      CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
      CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
      CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
      CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
      CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
      CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
103    CONTINUE
104 CONTINUE
RETURN
END

SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &
                WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
!     *** TR11=COS(2*PI/5), TI11=SIN(2*PI/5)
!     *** TR12=COS(4*PI/5), TI12=SIN(4*PI/5)
DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0, &
     0.95105651629515357212D0, &
     -0.8090169943749474241D0,0.58778525229247312917D0/
IF (IDO .NE. 2) GO TO 102
DO 101 K=1,L1
   TI5 = CC(2,2,K)-CC(2,5,K)
   TI2 = CC(2,2,K)+CC(2,5,K)
   TI4 = CC(2,3,K)-CC(2,4,K)
   TI3 = CC(2,3,K)+CC(2,4,K)
   TR5 = CC(1,2,K)-CC(1,5,K)
   TR2 = CC(1,2,K)+CC(1,5,K)
   TR4 = CC(1,3,K)-CC(1,4,K)
   TR3 = CC(1,3,K)+CC(1,4,K)
   CH(1,K,1) = CC(1,1,K)+TR2+TR3
   CH(2,K,1) = CC(2,1,K)+TI2+TI3
   CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
   CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
   CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
   CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
   CR5 = TI11*TR5+TI12*TR4
   CI5 = TI11*TI5+TI12*TI4
   CR4 = TI12*TR5-TI11*TR4
   CI4 = TI12*TI5-TI11*TI4
   CH(1,K,2) = CR2-CI5
   CH(1,K,5) = CR2+CI5
   CH(2,K,2) = CI2+CR5
   CH(2,K,3) = CI3+CR4
   CH(1,K,3) = CR3-CI4
   CH(1,K,4) = CR3+CI4
   CH(2,K,4) = CI3-CR4
   CH(2,K,5) = CI2-CR5
101 CONTINUE
RETURN
102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI5 = CC(I,2,K)-CC(I,5,K)
      TI2 = CC(I,2,K)+CC(I,5,K)
      TI4 = CC(I,3,K)-CC(I,4,K)
      TI3 = CC(I,3,K)+CC(I,4,K)
      TR5 = CC(I-1,2,K)-CC(I-1,5,K)
      TR2 = CC(I-1,2,K)+CC(I-1,5,K)
      TR4 = CC(I-1,3,K)-CC(I-1,4,K)
      TR3 = CC(I-1,3,K)+CC(I-1,4,K)
      CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
      CH(I,K,1) = CC(I,1,K)+TI2+TI3
      CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
      CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
      CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
      CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
      CR5 = TI11*TR5+TI12*TR4
      CI5 = TI11*TI5+TI12*TI4
      CR4 = TI12*TR5-TI11*TR4
      CI4 = TI12*TI5-TI11*TI4
      DR3 = CR3-CI4
      DR4 = CR3+CI4
      DI3 = CI3+CR4
      DI4 = CI3-CR4
      DR5 = CR2+CI5
      DR2 = CR2-CI5
      DI5 = CI2-CR5
      DI2 = CI2+CR5
      CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
      CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
      CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
      CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
      CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
      CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
      CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
      CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
103    CONTINUE
104 CONTINUE
RETURN
END
