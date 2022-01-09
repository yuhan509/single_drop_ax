subroutine getcar_cl(i,s,j,car)
use mod_elasbody
! use acar
integer,intent(in) :: i,j
DOUBLE PRECISION,intent(in) :: s
DOUBLE PRECISION,intent(out) :: car!, dcards, d2cards2
DOUBLE PRECISION :: temp

temp = s-ps1(i)

car = acar_cl(j,1,i)*temp**3 + acar_cl(j,2,i)*temp*temp &
 + acar_cl(j,3,i)*temp
If(j==i) car = car + 1d0

!dcards = 3d0*acar_cl(j,1,i)*temp**2 + 2d0*acar_cl(j,2,i)*temp &
! + acar_cl(j,3,i)
!
!d2cards2 = 6d0*acar_cl(j,1,i)*temp + 2d0*acar_cl(j,2,i)

END SUBROUTINE



subroutine getcar_na(i,s,j,car,dcards,d2cards2)
use mod_elasbody
! use acar
integer,intent(in) :: i,j
DOUBLE PRECISION,intent(in) :: s
DOUBLE PRECISION,intent(out) :: car, dcards, d2cards2
DOUBLE PRECISION :: temp

temp = s-ps1(i)

car = acar_na(j,1,i)*temp**3 + acar_na(j,2,i)*temp*temp &
 + acar_na(j,3,i)*temp
If(j==i) car = car + 1d0

dcards = 3d0*acar_na(j,1,i)*temp**2 + 2d0*acar_na(j,2,i)*temp &
 + acar_na(j,3,i)

d2cards2 = 6d0*acar_na(j,1,i)*temp + 2d0*acar_na(j,2,i)

END SUBROUTINE

!
!subroutine sb_getcar(i,s,j,car)
!use mod_singbody
!! use sb_acar
!integer,intent(in) :: i,j
!DOUBLE PRECISION,intent(in) :: s
!DOUBLE PRECISION,intent(out) :: car
!DOUBLE PRECISION :: temp
!
!temp = s-psb(i)
!car = sb_acar(j,1,i)*temp**3d0+sb_acar(j,2,i)*temp*temp&
!+sb_acar(j,3,i)*temp
!If(j==i) car=car+1d0
!
!END SUBROUTINE
