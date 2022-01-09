DOUBLE PRECISION FUNCTION for_r_1o(i,k)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,k
DOUBLE PRECISION :: dlds,Mrr,Mrz,sforce_r,sforce_z

dlds = dlds1o(k,i) 
sforce_r = sforce_r1o(k,i)
sforce_z = sforce_z1o(k,i)
Mrr = Mrr1o(k,i)
Mrz = Mrz1o(k,i)  

for_r_1o = (sforce_r*Mrr+sforce_z*Mrz)*dlds 

END FUNCTION


DOUBLE PRECISION FUNCTION for_z_1o(i,k)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,k
DOUBLE PRECISION :: dlds,Mzr,Mzz,sforce_r,sforce_z

dlds = dlds1o(k,i) 
sforce_r = sforce_r1o(k,i)
sforce_z = sforce_z1o(k,i)
Mzr = Mzr1o(k,i)
Mzz = Mzz1o(k,i)  

for_z_1o = (sforce_r*Mzr+sforce_z*Mzz)*dlds 

END FUNCTION


DOUBLE PRECISION FUNCTION for_r_1o_s(i,k)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,k
DOUBLE PRECISION :: dlds,sMrr,sMrz,sforce_r_s,sforce_z_s

dlds = dlds1o_s(k,i) 
sforce_r_s = sforce_r1o_s(k,i)
sforce_z_s = sforce_z1o_s(k,i)
sMrr = Mrr1o_s(k,i)
sMrz = Mrz1o_s(k,i)  

for_r_1o_s = (sforce_r_s * sMrr + sforce_z_s * sMrz)*dlds 

END FUNCTION


DOUBLE PRECISION FUNCTION for_z_1o_s(i,k)
! use r1,z1,PI
! doulbe-layer potential integrand
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,k
DOUBLE PRECISION :: dlds,sMzr,sMzz,sforce_r_s,sforce_z_s

dlds = dlds1o_s(k,i) 
sforce_r_s = sforce_r1o_s(k,i)
sforce_z_s = sforce_z1o_s(k,i)
sMzr = Mzr1o_s(k,i)
sMzz = Mzz1o_s(k,i)  

for_z_1o_s = (sforce_r_s * sMzr + sforce_z_s * sMzz) * dlds 
!write(82,*) 'i=',i,'k=',k,'sMzr=',sMzr,'sMzz=',sMzz,'sforce_r_s=',sforce_r_s,&
!'sforce_z_s=',sforce_z_s,'dlds=',dlds
END FUNCTION

!!///////////////// U Node J Related ///////////////////////////

DOUBLE PRECISION FUNCTION Ur_coefUrj1o(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,Mrr,Mrz,tmp1,tmp2

dlds = dlds1o(k,i) 
Mrr = Mrr1o(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
Mrz = Mrz1o(k,i)
tmp1 = delForce_r_coefUrj1o(j,k,i)
tmp2 = delForce_z_coefUrj1o(j,k,i)

Ur_coefUrj1o = (tmp1 * Mrr + tmp2 * Mrz)*dlds

END FUNCTION


DOUBLE PRECISION FUNCTION Ur_coefUzj1o(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,Mrr,Mrz,tmp3,tmp4

dlds = dlds1o(k,i) 
Mrr = Mrr1o(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
Mrz = Mrz1o(k,i)
tmp3 = delForce_r_coefUzj1o(j,k,i)
tmp4 = delForce_z_coefUzj1o(j,k,i)

Ur_coefUzj1o = (tmp3 * Mrr + tmp4 * Mrz)*dlds

END FUNCTION


DOUBLE PRECISION FUNCTION Uz_coefUrj1o(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,Mzr,Mzz,tmp1,tmp2

dlds = dlds1o(k,i) 
Mzr = Mzr1o(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
Mzz = Mzz1o(k,i)
tmp1 = delForce_r_coefUrj1o(j,k,i)
tmp2 = delForce_z_coefUrj1o(j,k,i)

Uz_coefUrj1o = (tmp1 * Mzr + tmp2 * Mzz)*dlds

END FUNCTION


DOUBLE PRECISION FUNCTION Uz_coefUzj1o(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,Mzr,Mzz,tmp3,tmp4

dlds = dlds1o(k,i) 
Mzr = Mzr1o(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
Mzz = Mzz1o(k,i)
tmp3 = delForce_r_coefUzj1o(j,k,i)
tmp4 = delForce_z_coefUzj1o(j,k,i)

Uz_coefUzj1o = (tmp3 * Mzr + tmp4 * Mzz)*dlds

END FUNCTION

!!//speical
DOUBLE PRECISION FUNCTION Ur_coefUrj1o_s(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,sMrr,sMrz,tmp1,tmp2

dlds = dlds1o_s(k,i) 
sMrr = Mrr1o_s(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
sMrz = Mrz1o_s(k,i)
tmp1 = delForce_r_coefUrj1o_s(j,k,i)
tmp2 = delForce_z_coefUrj1o_s(j,k,i)

Ur_coefUrj1o_s = (tmp1 * sMrr + tmp2 * sMrz)*dlds

END FUNCTION

!!//speical
DOUBLE PRECISION FUNCTION Ur_coefUzj1o_s(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,sMrr,sMrz,tmp3,tmp4

dlds = dlds1o_s(k,i) 
sMrr = Mrr1o_s(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
sMrz = Mrz1o_s(k,i)
tmp3 = delForce_r_coefUzj1o_s(j,k,i)
tmp4 = delForce_z_coefUzj1o_s(j,k,i)

Ur_coefUzj1o_s = (tmp3 * sMrr + tmp4 * sMrz)*dlds

END FUNCTION


!!//speical
DOUBLE PRECISION FUNCTION Uz_coefUrj1o_s(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,sMzr,sMzz,tmp1,tmp2

dlds = dlds1o_s(k,i) 
sMzr = Mzr1o_s(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
sMzz = Mzz1o_s(k,i)
tmp1 = delForce_r_coefUrj1o_s(j,k,i)
tmp2 = delForce_z_coefUrj1o_s(j,k,i)

Uz_coefUrj1o_s = (tmp1 * sMzr + tmp2 * sMzz)*dlds

END FUNCTION


!!//speical
DOUBLE PRECISION FUNCTION Uz_coefUzj1o_s(j,i,k)
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: j,i,k
DOUBLE PRECISION :: dlds,sMzr,sMzz,tmp3,tmp4

dlds = dlds1o_s(k,i) 
sMzr = Mzr1o_s(k,i) !!// Mrr,Mrz,Mzz,Mzr are depending on i0
sMzz = Mzz1o_s(k,i)
tmp3 = delForce_r_coefUzj1o_s(j,k,i)
tmp4 = delForce_z_coefUzj1o_s(j,k,i)

Uz_coefUzj1o_s = (tmp3 * sMzr + tmp4 * sMzz)*dlds

END FUNCTION

