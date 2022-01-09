SUBROUTINE Ucoef_gax(i0,j,curr1,curz1,cuzr1,cuzz1)
! use N
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0,j
! j for cardinal !i0 for the kernal node
DOUBLE PRECISION,INTENT(OUT) :: curr1,curz1,cuzr1,cuzz1
!z1,r1 corresponds to the single-layer potential 
!the integral over all elements
DOUBLE PRECISION :: res 
INTEGER :: i ! ith element

DOUBLE PRECISION,EXTERNAL :: Ur_coefUrj1o,Ur_coefUzj1o,&
Uz_coefUrj1o,Uz_coefUzj1o,Ur_coefUrj1o_s,Ur_coefUzj1o_s,&
Uz_coefUrj1o_s,Uz_coefUzj1o_s

curr1 = 0d0
curz1 = 0d0
cuzr1 = 0d0
cuzz1 = 0d0

Do i=1,N
    !print*, '   i0',i0,'i',i
    CALL fm_gq_cubsj_o(Ur_coefUrj1o,i,j,res)
    curr1 = curr1 + res
    CALL fm_gq_cubsj_o(Ur_coefUzj1o,i,j,res)
    curz1 = curz1 + res
    CALL fm_gq_cubsj_o(Uz_coefUrj1o,i,j,res)
    cuzr1 = cuzr1 + res
    CALL fm_gq_cubsj_o(Uz_coefUzj1o,i,j,res)
    cuzz1 = cuzz1 + res
    
    IF(i0==i .OR. i0==(i+1)) THEN
        CALL fm_sgq_cubsj_o(Ur_coefUrj1o_s,i,i0,j,res)
        curr1 = curr1 + res
        CALL fm_sgq_cubsj_o(Ur_coefUzj1o_s,i,i0,j,res)
        curz1 = curz1 + res
        CALL fm_sgq_cubsj_o(Uz_coefUrj1o_s,i,i0,j,res)
        cuzr1 = cuzr1 + res
        CALL fm_sgq_cubsj_o(Uz_coefUzj1o_s,i,i0,j,res)
        cuzz1 = cuzz1 + res        
    ENDIF
    
ENDDO

END SUBROUTINE


SUBROUTINE Fcoef_gax(i0,fr1,fz1) 
! use N
use mod_elasbody
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0
! j for cardinal !i0 for the kernal node
DOUBLE PRECISION,INTENT(OUT) :: fr1,fz1
!z1,r1 corresponds to the single-layer potential 
!the integral over all elements
INTEGER :: i ! ith element
DOUBLE PRECISION :: res !,temp
DOUBLE PRECISION,EXTERNAL :: for_r_1o,for_z_1o,&
for_r_1o_s,for_z_1o_s


fz1 = 0d0
fr1 = 0d0

Do i=1,N
    !print*, '   i0',i0,'i'
!    IF (i==39) THEN
!        write(*,*) 'i = 39'
!    ENDIF 
    CALL fm_gq_cubs_o(for_r_1o,i,res)
    fr1 = fr1 + res
    !write(*,*) 'i0',i0,'i',i,'ai',ai,'bi',bi,'ez1',ez1,'fz1',fz1
    CALL fm_gq_cubs_o(for_z_1o,i,res)
    fz1 = fz1 + res
!    if(i0 == N+1) THEN
!       write(80,*) 'i0 =',i0,'i =',i,'fz1_res =',res
!    ENDIF 
    
    IF(i0==i .OR. i0==(i+1)) THEN
        CALL fm_sgq_cubs_o(for_r_1o_s,i,i0,res)
        fr1 = fr1 + res
        CALL fm_sgq_cubs_o(for_z_1o_s,i,i0,res)  
        fz1 = fz1 + res
!        write(81,*) 'i0 =',i0,'i =',i,'fz1_s_res =',res
    ENDIF
    
ENDDO

END SUBROUTINE
