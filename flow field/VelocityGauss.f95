SUBROUTINE fm_gq_cubs_o(func,i,res)
!! a and b are the lower and upper borders, a=i,b=i+1
USE mod_GAUSSQxw   
USE mod_elasbody 
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
INTEGER,INTENT(IN) :: i!,j
! i for ith element
DOUBLE PRECISION,INTENT(OUT) :: res
INTEGER :: k
DOUBLE PRECISION :: xr

res = 0d0
xr = (ps1(i+1)-ps1(i))*.5d0

Do k=1,qm
    res = res + qw(k)*func(i,k)
ENDDO

res = res*xr
!write(*,*) "regular_gauss xr",xr
END SUBROUTINE




SUBROUTINE fm_gq_cubsj_o(func,i,j,res)
USE mod_GAUSSQxw
USE mod_elasbody 
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
INTEGER,INTENT(IN) :: i,j
! i for ith element
DOUBLE PRECISION,INTENT(OUT) :: res
INTEGER :: k
DOUBLE PRECISION :: xr

res = 0d0
!!// integral range: from ps1(i) to ps1(i+1)
xr = (ps1(i+1)-ps1(i))*.5d0
Do k=1,qm
    res = res + qw(k)*func(j,i,k)
ENDDO

res = res*xr
!write(*,*) "regular_gauss xr",xr
END SUBROUTINE


