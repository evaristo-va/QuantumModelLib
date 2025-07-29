!===========================================================================
!===========================================================================
!This file is part of ModelLib.
!
!    ModelLib is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ModelLib is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ModelLib.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2016 David Lauvergnat [1]
!      with contributions of:
!        Félix MOUHAT [2]
!        Liang LIANG [3]
!        Emanuele MARSILI [1,4]
!        Evaristo Villaseco Arribas [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
![3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
![4]: Durham University, Durham, UK
![5]: Department of Physics, Rutgers University, Newark, New Jersey 07102, USA
!* Originally, it has been developed during the Quantum-Dynamics E-CAM project :
!     https://www.e-cam2020.eu/quantum-dynamics
!
!===========================================================================
!===========================================================================
!
!> @brief Module which makes the initialization, calculation of the Benzene potentials (value and gradient).
!! @brief Reference: Thomas J. Penfold and Graham A. Worth. A model Hamiltonian to simulate the complex photochemistry of benzene II. J. Chem. Phys. 131, 064303 (2009), https://doi.org/10.1063/1.3197555
!> @author Evaristo Villaseco Arribas
!! @date 09/02/2023
!!
MODULE QML_Benzene_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Benzene parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Benzene_t
    PRIVATE

    real(kind=Rkind), allocatable :: E(:)
    real(kind=Rkind), allocatable :: K(:,:)
    real(kind=Rkind), allocatable :: EPS(:,:)
    real(kind=Rkind), allocatable :: L(:,:,:)
    real(kind=Rkind), allocatable :: G(:,:,:)
    real(kind=Rkind), allocatable :: MU(:,:,:,:)
    real(kind=Rkind), allocatable :: d0(:,:)
    real(kind=Rkind), allocatable :: a0(:,:)
    real(kind=Rkind), allocatable :: q(:,:)
    real(kind=Rkind), allocatable :: e0(:,:)
    real(kind=Rkind), allocatable :: F(:)
    real(kind=Rkind), allocatable :: M(:)
    real(kind=Rkind), allocatable :: Xref(:)
    real(kind=Rkind), allocatable :: Ktransf(:,:)


   CONTAINS
    PROCEDURE :: EvalPot_QModel    => EvalPot_QML_Benzene
    PROCEDURE :: Write_QModel      => Write_QML_Benzene
    PROCEDURE :: Cart_to_Q_QModel  => Cart_to_Q_QML_Benzene
  END TYPE QML_Benzene_t

  PUBLIC :: QML_Benzene_t,Init_QML_Benzene

  CONTAINS
!> @brief Function which makes the initialization of the Benzene parameters.
!!
!! @param QModel             TYPE(QML_Benzene_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):   type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Benzene(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat, file_open2
    USE QMLLib_UtilLib_m, ONLY : make_QMLInternalFileName
    IMPLICIT NONE

    TYPE (QML_Benzene_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    character(len=:),            allocatable     :: FileName
    character(len=:),            allocatable     :: FileName2
    integer                                      :: nio_fit,i,j,k,l,s,num
    character(len=30)                            :: temp

    !----- converstions --------------------------------------------------
    real(kind=Rkind), parameter :: eV_to_au    = 0.03674931_Rkind
    real(kind=Rkind), parameter :: amu_to_au   = 1822.888486192_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Benzene'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 7
    QModel%ndim     = 36
    QModel%pot_name = 'Benzene'
    QModel%ndimCart = 36
    QModel%ndimQ    = 36


    IF (allocated(QModel%E)) deallocate(QModel%E)
    IF (allocated(QModel%K)) deallocate(QModel%K)
    IF (allocated(QModel%EPS)) deallocate(QModel%EPS)
    IF (allocated(QModel%L)) deallocate(QModel%L)
    IF (allocated(QModel%G)) deallocate(QModel%G)
    IF (allocated(QModel%MU)) deallocate(QModel%MU)
    IF (allocated(QModel%d0)) deallocate(QModel%d0)
    IF (allocated(QModel%a0)) deallocate(QModel%a0)
    IF (allocated(QModel%q)) deallocate(QModel%q)
    IF (allocated(QModel%e0)) deallocate(QModel%e0)

    IF (allocated(QModel%F)) deallocate(QModel%F)
    IF (allocated(QModel%M)) deallocate(QModel%M)
    IF (allocated(QModel%Xref)) deallocate(QModel%Xref)
    IF (allocated(QModel%Ktransf)) deallocate(QModel%Ktransf)

    allocate(QModel%E(QModel%nsurf))
    allocate(QModel%K(QModel%nsurf,QModel%ndim))
    allocate(QModel%EPS(QModel%nsurf,QModel%ndim))
    allocate(QModel%L(QModel%nsurf,QModel%nsurf,QModel%ndim))
    allocate(QModel%G(QModel%nsurf,QModel%ndim,QModel%ndim))
    allocate(QModel%MU(QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim))
    allocate(QModel%d0(QModel%nsurf,QModel%ndim))
    allocate(QModel%a0(QModel%nsurf,QModel%ndim))
    allocate(QModel%q(QModel%nsurf,QModel%ndim))
    allocate(QModel%e0(QModel%nsurf,QModel%ndim))
    allocate(QModel%F(QModel%ndim))
    allocate(QModel%M(QModel%ndim))
    allocate(QModel%Xref(QModel%ndim))
    allocate(QModel%Ktransf(QModel%ndim,QModel%ndim))

    QModel%E=0.0_Rkind
    QModel%K=0.0_Rkind
    QModel%EPS=0.0_Rkind
    QModel%L=0.0_Rkind
    QModel%G=0.0_Rkind
    QModel%MU=0.0_Rkind
    QModel%d0=0.0_Rkind
    QModel%a0=0.0_Rkind
    QModel%q=0.0_Rkind
    QModel%e0=0.0_Rkind
    QModel%F=0.0_Rkind
    QModel%M=0.0_Rkind
    QModel%Xref=0.0_Rkind
    QModel%Ktransf=0.0_Rkind

    !-------------READ BENZENE PARAMETERS FROM EXTERNAL INPUT------------!
    FileName = make_QMLInternalFileName('InternalData/BENZENE/R0.txt')
    CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)

      !------Ref geometry and Masses----!
      DO i=1,QModel%ndim
         read(nio_fit,*) QModel%Xref(i), QModel%M(i)
      ENDDO

    close(nio_fit)
    deallocate(FileName)

    FileName2 = make_QMLInternalFileName('InternalData/BENZENE/V0.txt')
    CALL file_open2(name_file=FileName2,iunit=nio_fit,old=.TRUE.)

      !------NMCtoCART transf Matrix----!
      DO i=1,QModel%ndim
         read(nio_fit,*) (QModel%Ktransf(i,j), j=1,QModel%ndim)
      ENDDO

    close(nio_fit)
    deallocate(FileName2)

    FileName = make_QMLInternalFileName('InternalData/BENZENE/param_VCM.dat')
    CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)

      !------------Ref energies---------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, QModel%E(j+1)
      ENDDO

      !--------Kappa term (KPi*Q)-------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, QModel%K(j+1,k)
      ENDDO

      !----Kquartic term (EPS*Q^4)---------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, QModel%EPS(j+1,k)
      ENDDO

      !-----Lambda term (Lij*Q)-----------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, l, QModel%L(j+1,k+1,l)
         QModel%L(k+1,j+1,l)=QModel%L(j+1,k+1,l)
      ENDDO

      !---Gamma term (GM^i_{ab}*Q_a*Q_b)-----------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, l, QModel%G(j+1,k,l)
      ENDDO

      !----------Morse parameters---------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, QModel%d0(j+1,k), QModel%a0(j+1,k), QModel%q(j+1,k), QModel%e0(j+1,k)
      ENDDO

      !---MU term (MU^{ij}_{ab}*Q_a*Q_b)-----------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, l, s, QModel%MU(j+1,k+1,l,s)
      ENDDO

      !-----------NM frequencies----------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) k, QModel%F(k)
      ENDDO

    close(nio_fit)
    deallocate(FileName)

    !--------------CONVERSION TO ATOMIC UNITS OF PARAMETERS-------------!
    QModel%M   = Qmodel%M   * amu_to_au

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Benzene
!> @brief Subroutine wich prints the QML_Benzene parameters.
!!
!! @param QModel            CLASS(QML_Benzene_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Benzene(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Benzene_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'model name: ', QModel%pot_name
    write(nio,*) 'Orignial model from:'
    write(nio,*) '   Mariana Assmann, Horst Köppel, and Spiridoula Matsika,'
    write(nio,*) '   J. Phys. Chem. A 2015, 119, 866−875.'
    write(nio,*)  '  DOI: 10.1021/jp512221x'
    write(nio,*) 'and modifed in:'
    write(nio,*) '   Patricia Vindel Zandbergen, Spiridoula Matsika, and Neepa T. Maitra'
    write(nio,*) '   J. Phys. Chem. Lett. 2022, 13, 7, 1785–1790'
    write(nio,*) '   https://doi.org/10.1021/acs.jpclett.1c04132'

  END SUBROUTINE Write_QML_Benzene
!> @brief Subroutine wich calculates the Benzene potential with derivatives up to the 1d order.
!!
!! @param QModel             CLASS(QML_Benzene_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         Potential with derivatives,.
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to secify the derivative order:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Benzene(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Benzene_t),  intent(in)    :: QModel
  TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv
    integer                             :: i,j,i_dof,j_dof

    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)        :: exponential

   Mat_OF_PotDia(:,:)=0.0_Rkind

   !Computation of Diagonal terms H_{ii}
   DO i=1,QModel%nsurf
     ! Reference energy + harmonic term (E_i + (1/2) * \Sum_nu F^{nu} * Q^{nu} **2 )
     Mat_OF_PotDia(i,i) = QModel%E(i) + 0.5_Rkind * DOT_PRODUCT(QModel%F(:),dnQ(:)**2) 
     ! Linear kappa term (K_i^{nu} * Q^{nu})
     Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + DOT_PRODUCT(QModel%K(i,:),dnQ(:))
     DO i_dof=1,QModel%ndim
       ! Morse terms
       IF (QModel%d0(i,i_dof) /= 0.0_Rkind) THEN
         exponential = exp( QModel%a0(i,i_dof) * ( dnQ(i_dof) - QModel%q(i,i_dof) ) )
         Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + QModel%d0(i,i_dof) * ( exponential - 1)**2 + QModel%e0(i,i_dof)
       ENDIF
       DO j_dof=1,QModel%ndim
         ! Quadratic on-diagonal and off-diagonal terms ((1/2) * \Sum_nu\Sum_mu G_{i}^{munu} * Q^{nu} * Q^{mu})
         Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + 0.5_Rkind * QModel%G(i,i_dof,j_dof)*dnQ(i_dof)*dnQ(j_dof)
       ENDDO
     ENDDO
     ! Quartic term ((1/24) * EPS_nu^{i} * Q^{nu} ** 4
     Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + (1.0_Rkind/24.0_Rkind) * DOT_PRODUCT(QModel%EPS(i,:),dnQ(:)**4)
   ENDDO

   ! Off diagonal terms H_{ij}
   DO i=1,QModel%nsurf
       DO j=i+1,QModel%nsurf
         ! Linear lamba term (L_{ij}^{nu} * Q^{nu})
         Mat_OF_PotDia(i,j) = DOT_PRODUCT(QModel%L(i,j,:),dnQ(:))
         DO i_dof=1,QModel%ndim
           DO j_dof=1,QModel%ndim
           ! quadratic mu term (MU_{ij}^{numu} * Q^{nu} * Q^{mu})
           Mat_OF_PotDia(i,j) = Mat_OF_PotDia(i,j) + 0.5_Rkind * QModel%MU(i,j,i_dof,j_dof)*dnQ(i_dof)*dnQ(j_dof)
           ENDDO
         ENDDO
         Mat_OF_PotDia(j,i) = Mat_OF_PotDia(i,j)
       ENDDO
   ENDDO

   CALL dealloc_dnS(exponential)

  END SUBROUTINE EvalPot_QML_Benzene

  SUBROUTINE Cart_TO_Q_QML_Benzene(QModel,dnX,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Benzene_t), intent(in)     :: QModel
    TYPE (dnS_t),        intent(in)     :: dnX(:,:)
    TYPE (dnS_t),        intent(inout)  :: dnQ(:)
    integer,             intent(in)     :: nderiv
    integer :: i,j,iat,coord

    DO i=1,QModel%ndim
      dnQ(i)=0.0_Rkind
      DO iat=1,size(dnX,dim=2)
        DO coord=1,3
          j=3*(iat-1)+coord
          dnQ(i)=dnQ(i)+QModel%Ktransf(j,i)*sqrt(Qmodel%M(j))*(dnX(coord,iat)-Qmodel%Xref(j))
        ENDDO
      ENDDO
      dnQ(i)=dnQ(i)*sqrt(QModel%F(i))
    ENDDO


  END SUBROUTINE Cart_to_Q_QML_Benzene


END MODULE QML_Benzene_m
