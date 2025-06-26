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
!> @brief Module which makes the initialization, calculation of the Thiophene potentials (value and gradient).
!! @brief Reference: model fitted by Jiri Suchan
!> @author Evaristo Villaseco Arribas
!! @date 09/02/2023
!!
MODULE QML_Thiophene_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Thiophene parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Thiophene_t
    PRIVATE

    real(kind=Rkind), allocatable :: F(:)
    real(kind=Rkind), allocatable :: M(:)
    real(kind=Rkind), allocatable :: E(:)
    real(kind=Rkind), allocatable :: KP(:,:)
    real(kind=Rkind), allocatable :: L(:,:,:)
    real(kind=Rkind), allocatable :: Xref(:)
    real(kind=Rkind), allocatable :: Ktransf(:,:)


   CONTAINS
    PROCEDURE :: EvalPot_QModel    => EvalPot_QML_Thiophene
    PROCEDURE :: Write_QModel      => Write_QML_Thiophene
    PROCEDURE :: Cart_to_Q_QModel  => Cart_to_Q_QML_Thiophene
  END TYPE QML_Thiophene_t

  PUBLIC :: QML_Thiophene_t,Init_QML_Thiophene

  CONTAINS
!> @brief Function which makes the initialization of the Thiophene parameters.
!!
!! @param QModel             TYPE(QML_Thiophene_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):   type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Thiophene(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat, file_open2
    USE QMLLib_UtilLib_m, ONLY : make_QMLInternalFileName
    IMPLICIT NONE

    TYPE (QML_Thiophene_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    character(len=:),            allocatable     :: FileName
    character(len=:),            allocatable     :: FileName2
    integer                                      :: nio_fit,i,j,k,l,num
    character(len=30)                            :: temp

    !----- converstions --------------------------------------------------
    real(kind=Rkind), parameter :: eV_to_au    = 0.03674931_Rkind
    real(kind=Rkind), parameter :: amu_to_au   = 1822.888486192_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Thiophene'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 25
    QModel%ndim     = 27
    QModel%pot_name = 'Thiophene'
    QModel%ndimCart = 27
    QModel%ndimQ    = 27


    IF (allocated(QModel%M)) deallocate(QModel%M)
    IF (allocated(QModel%F)) deallocate(QModel%F)
    IF (allocated(QModel%E)) deallocate(QModel%E)
    IF (allocated(QModel%KP)) deallocate(QModel%KP)
    IF (allocated(QModel%L)) deallocate(QModel%L)
    IF (allocated(QModel%Xref)) deallocate(QModel%Xref)
    IF (allocated(QModel%Ktransf)) deallocate(QModel%Ktransf)

    allocate(QModel%F(QModel%ndim))
    allocate(QModel%M(QModel%ndim))
    allocate(QModel%E(QModel%nsurf))
    allocate(QModel%KP(QModel%nsurf,QModel%ndim))
    allocate(QModel%L(QModel%nsurf,QModel%nsurf,QModel%ndim))
    allocate(QModel%Xref(QModel%ndim))
    allocate(QModel%Ktransf(QModel%ndim,QModel%ndim))

    QModel%F=0.0_Rkind
    QModel%M=0.0_Rkind
    QModel%E=0.0_Rkind
    QModel%KP=0.0_Rkind
    QModel%L=0.0_Rkind
    QModel%Xref=0.0_Rkind
    QModel%Ktransf=0.0_Rkind

    !-------------READ THIOPHENE PARAMETERS FROM EXTERNAL INPUT------------!
    FileName = make_QMLInternalFileName('InternalData/THIOPHENE/R0.txt')
    CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)

      !------Ref geometry and Masses----!
      DO i=1,QModel%ndim
         read(nio_fit,*) QModel%Xref(i), QModel%M(i)
      ENDDO

    close(nio_fit)
    deallocate(FileName)

    FileName2 = make_QMLInternalFileName('InternalData/THIOPHENE/V0.txt')
    CALL file_open2(name_file=FileName2,iunit=nio_fit,old=.TRUE.)

      !------NMCtoCART transf Matrix----!
      DO i=1,QModel%ndim
         read(nio_fit,*) (QModel%Ktransf(i,j), j=1,QModel%ndim)
      ENDDO

    close(nio_fit)
    deallocate(FileName2)

    FileName = make_QMLInternalFileName('InternalData/THIOPHENE/param_VCM.dat')
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
         read(nio_fit,*) j, k, QModel%KP(j+1,k)
      ENDDO

      !-----Lambda term (Lij*Q)-----------!
      read(nio_fit,*) temp
      read(nio_fit,*) num
      DO i=1,num
         read(nio_fit,*) j, k, l, QModel%L(j+1,k+1,l)
         QModel%L(k+1,j+1,l)=QModel%L(j+1,k+1,l)
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
    !QModel%E  = QModel%E  * eV_to_au
    !QModel%KP = QModel%KP * eV_to_au
    !QModel%L  = QModel%L  * eV_to_au
    QModel%M = QModel%M * amu_to_au


    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Thiophene
!> @brief Subroutine wich prints the QML_Thiophene parameters.
!!
!! @param QModel            CLASS(QML_Thiophene_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Thiophene(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Thiophene_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'model name: ', QModel%pot_name
    write(nio,*) 'model fitted by Jiri Suchan'

  END SUBROUTINE Write_QML_Thiophene
!> @brief Subroutine wich calculates the Thiophene potential with derivatives up to the 1d order.
!!
!! @param QModel             CLASS(QML_Thiophene_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         Potential with derivatives,.
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to secify the derivative order:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Thiophene(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Thiophene_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv
    integer                             :: i,j,i_dof

    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)        :: exponential

   Mat_OF_PotDia(:,:)=0.0_Rkind

   !Computation of H0
   DO i=1,QModel%nsurf
     Mat_OF_PotDia(i,i) = QModel%E(i)
     !Harmonic modes
     DO i_dof=1,QModel%ndim
       Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + 0.5_Rkind * QModel%F(i_dof) * dnQ(i_dof)**2
     ENDDO
   ENDDO

   !Adding W matrix
   DO i=1,QModel%nsurf
     Mat_OF_PotDia(i,i) = Mat_OF_PotDia(i,i) + DOT_PRODUCT(QModel%KP(i,:),dnQ(:)) 
       DO j=i+1,QModel%nsurf
         Mat_OF_PotDia(i,j) = DOT_PRODUCT(QModel%L(i,j,:),dnQ(:))
         Mat_OF_PotDia(j,i) = Mat_OF_PotDia(i,j)
       ENDDO
   ENDDO

   CALL dealloc_dnS(exponential)

  END SUBROUTINE EvalPot_QML_Thiophene

  SUBROUTINE Cart_TO_Q_QML_Thiophene(QModel,dnX,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Thiophene_t), intent(in)     :: QModel
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


  END SUBROUTINE Cart_to_Q_QML_Thiophene


END MODULE QML_Thiophene_m
