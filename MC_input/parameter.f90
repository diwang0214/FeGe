 module para
 implicit none
 integer,parameter::dimx=8,dimy=8,dimz=8,natom=3
 integer,parameter::nlatt=dimx*dimy*dimz
 integer,parameter::nloop=dimx*dimy*dimz                                  !the dimension of the lattice
                  
  !lattice neighbor parameter 
 !integer,dimension(dimx)::ui,di
 !integer,dimension(dimy)::uj,dj
 !integer,dimension(dimz)::uk,dk
 integer,dimension(nlatt)::u,d,l,r,lu,ld,ru,rd,uz,dz,uuz,udz,ddz,duz,luz,ldz,ruz,rdz,luuz,ludz,lddz,lduz,ruuz,rudz,rduz,rddz 
 ! all 26 neighbors of one lattice


 double precision,parameter::pi=3.14159265358979323
 double precision,dimension(1:3,1:nloop,1:natom)::spin,newspin,dspin
 double precision::sdelta=1.0
 double precision,dimension(1:nloop,1:natom)::de
 integer,dimension(1:dimx,1:dimy,1:dimz)::lattice


 !control parameter
 character(len=4)::sflag='flat'                            
! integer,parameter::uflag=0                                 !u_flag=0:no phonon,  1:continious phonon
 integer,parameter::steplen=200
 character(len=3)::fliptype='ran'
 integer,parameter::ultime=1
 integer,parameter::ubtime=1
 integer,parameter::ulen=5000,slen=5000
 integer::uadap=1,sadap=1


 !output parameter ######################################################
 !total magnetization parameter
 double precision::totm,totm2,tm,tmsqu,sus
 double precision::totmx,totmy,totmz
 !real,dimension(1:ntem)::magave,magdif


 !total energy parameter
 double precision::te,tote,linkde,dot1e,dot2e,dot1de,dot2de
 double precision::toten!,!linkde,dot1e,dot2e,dot1de,dot2de
! double precision::old1_energy,new1_energy,m2de,old2_energy, new2_energy,m3de
! double precision::old3_energy,new3_energy,m4de,old4_energy,&
! &new4_energy,m5de,m6de,old_energy,new_energy
 double precision::spheat,esqu                                !the specific heat and square of energy
 !##########################################################################
 !lattice neighbor parameter


  !five type of exchange interaction
 !double precision,dimension(1:1)::jexc
 double precision::jxy,jz
 double precision,dimension(1:3)::k31,k32,D3,D4    ! single ion anisotropy energy 

 !umode's parameter
 double precision::j31,umin=-2.0,umax=2.0,omig
 integer,parameter::nmod=2
 double precision,dimension(1:3)::hfield
 double precision::magmom

 !the parameter for calculating accepting rate
 integer(kind=4)::ulcount,ubcount,ucountl,uacceptl,ucountb,uacceptb,scount,sacpt
 integer(kind=4)::totcount,totacpt
 integer(kind=4),parameter::countstep=1000000000
 double precision::totrate

 !integer::ux,uy,uz,um
 double precision::ude2
 double precision::deltau,tu,totu,tusqu,psus
 double precision,dimension(1:nloop,1:nmod)::udis,unew
 double precision,dimension(1:nloop,1:nmod)::lamda
 double precision::dlamda

 

 !distribution parameter

  double precision::eave,eeave,mave,mmave,uave,uuave,emmave,meave,tmsque,tme,susot
  double precision::deltaul=1.0,muul=0.0,deltaub=1.0,muub=0.0
  double precision::pu=0.3,pl=0.2,beta0=2.0,beta1=2.0,probb,probl
  double precision::deltas=1.0,mus=0.0
  double precision::psmax=0.5,psmin=0.4,probs

!################################### parameters of replica exchange #############################

!parameter of MPI
integer::nproc,my_rank
integer::ierr,itag
integer::seed1,seedc,seedx
integer,allocatable::seed(:)

!parameters of replica exchange
integer::mybeta_num,bcast_rank
integer,parameter::sendsize=2,recvsize=2
double precision,allocatable::beta(:)
integer(kind=4),allocatable::exacpt(:)
integer(kind=4),allocatable::excount(:)
integer,allocatable::process(:)
integer,allocatable::betanum(:)
integer,allocatable::exbeta(:,:)
double precision,dimension(0:sendsize-1)::nsend,send_init,recv_init
double precision,dimension(0:recvsize-1)::nrecv
double precision,parameter::betamin=0.0386333d0,betamax=2.318d0
integer,parameter::nscheme=2,nrec=1
integer,parameter::nswp1=1000,nexch1=200
integer,parameter::nswp2=1000,nexch2=200
integer::nswp,nexch
integer,parameter::nequi1=200000,nequi2=1000
integer::irec,iswp,iexch,ischeme,ibeta,iequi,tarbeta1,tarbeta2
integer::exflag,exflag1,stat,irank
integer::betanum_temp,process_temp
 
!parameters of recursion
double precision::wk,wk_temp
double precision,allocatable::wb(:),wb_temp(:)
integer,parameter::initial_beta=1
integer,parameter::exchange_flag=1
integer,parameter::proflag=1
integer,parameter::initialsflag=0
integer,parameter::initialuflag=0
integer,parameter::initial_bedis=0
integer,parameter::sam_flag=0
integer,parameter::anis_flag=1
integer,parameter::uflag=0                   !u_flag=0:no phonon,  1:continious phonon
integer,parameter::hflag=0

!sample parameters
integer,parameter::nbeta_sam=1
integer,dimension(1:nbeta_sam)::beta_sam=(/279/)
double precision,dimension(1:3,1:nloop,1:natom)::spin_out
double precision,dimension(1:nloop,1:nmod)::udis_out

! 1. add for new lattice in 20190620
integer,parameter::mmldx=1, mmldy=1,mmldz=1
integer,dimension(-mmldx:mmldx,1:dimx)::di
integer,dimension(-mmldy:mmldy,1:dimy)::dj
integer,dimension(-mmldz:mmldz,1:dimz)::dk
integer,dimension(-mmldz:mmldz,-mmldy:mmldy,-mmldx:mmldx, nlatt)::mld

! 2. add for initial energy in 20190620
real(kind=8),dimension(1:natom,1:natom,-mmldz:mmldz,-mmldy:mmldy,-mmldx:mmldx)::jexc
real(kind=8),dimension(1:natom)::dote,doten!,do1e,do2e,do3e,do4e,do5e,do6e,do7e,do8e,do9e,do10e

end module

subroutine readincar
use para
implicit none
integer::i,j,l1,m,n!,temp11
open(2000,file='incar')
do i=1,natom,1
  do j=1,natom,1
    do l1=-mmldx,mmldx,1
      do m=-mmldy,mmldy,1
        do n=-mmldz,mmldz,1
        read(2000, "(1F9.2)") jexc(i,j,n,m,l1)
        write(*,*) i,j,n,m,l1
        end do !n
      end do !m
    end do !l1
  end do !j
end do !i

if(my_rank==0) write(*,*) "Read jexc is ok"

read(2000,*)D4
if(hflag==1)then
read(2000,*)hfield
read(2000,*)magmom
end if

close(2000)
end subroutine readincar
