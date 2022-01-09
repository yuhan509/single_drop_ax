SUBROUTINE fm_sgq_cubs_o(func,i,i0,ss)
!!! Integrand form: f(x)*log(x)=sum of f(x)
USE mod_GAUSSQxw   
USE mod_elasbody 
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
INTEGER,INTENT(IN) :: i,i0!,j
DOUBLE PRECISION,INTENT(OUT) :: ss
INTEGER :: ks

ss=0d0

Do ks = 1,sqm
    IF(i==i0) THEN
        ss = ss + sqw(ks) * func(i,ks)
    ELSE IF(i+1==i0) THEN
        ss = ss + sqw(ks) * func(i,ks+sqm)
    END IF     
ENDDO

ss = ss*(ps1(i+1)-ps1(i))
END SUBROUTINE



SUBROUTINE fm_sgq_cubsj_o(func,i,i0,j,ss)
!!! Integrand form: f(x)*log(x)=sum of f(x)
USE mod_GAUSSQxw   
USE mod_elasbody 
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL :: func
INTEGER,INTENT(IN) :: i0,i,j
DOUBLE PRECISION,INTENT(OUT) :: ss
INTEGER :: ks

ss=0d0

Do ks = 1,sqm
    IF(i==i0) THEN
        ss = ss + sqw(ks) * func(j,i,ks)
    ELSE IF(i+1==i0) THEN
        ss = ss + sqw(ks) * func(j,i,ks+sqm)
    END IF     
ENDDO

ss = ss*(ps1(i+1)-ps1(i))
END SUBROUTINE
