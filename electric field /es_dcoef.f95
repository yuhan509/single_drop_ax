SUBROUTINE es_dcoef_gax(i0,j,c1)
use mod_elasbody
! use N
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0,j
DOUBLE PRECISION,INTENT(OUT) :: c1
DOUBLE PRECISION :: ai,bi,s1,s1s !,temp
DOUBLE PRECISION,EXTERNAL :: es_intg_gax1o,&
es_intg_sgax1
INTEGER :: i !ith element
c1 = 0d0
Do i=1,N
    ai = ps1(i)
    bi = ps1(i+1)
    CALL es_gq_gax_o(es_intg_gax1o,i,ai,bi,j,i0,s1)
    !write(*,*) 'i0',i0,'j',j,'s1',s1,'c1',c1
    c1 = c1 + s1
    !write(*,*) 'after es_gq_gax_o'
    IF(i0==i .OR. i0==(i+1)) THEN        
         CALL es_sgq_gax(es_intg_sgax1,i,ai,bi,j,i0,s1s)
         !write(*,*) 'i0',i0,'j',j,'s1s',s1s,'c1',c1
         c1 = c1 + s1s
    ENDIF
ENDDO
END SUBROUTINE



SUBROUTINE es_dcoef_dgdnax(i0,c2)
use mod_elasbody
! use N
IMPLICIT NONE
INTEGER,INTENT(IN) :: i0
DOUBLE PRECISION,INTENT(OUT) :: c2
DOUBLE PRECISION :: ai,bi,s2,s2s !,temp
DOUBLE PRECISION,EXTERNAL ::es_intg_dgdnax1o&
,es_intg_sdgdnax1
INTEGER :: i !ith element

c2 = 0d0
Do i=1,N
    ai = ps1(i)
    bi = ps1(i+1)
    CALL es_gq_dgdnax_o(es_intg_dgdnax1o,i,ai,bi,i0,s2)!
    !write(*,*) 'i0',i0,'i',i,'ai',ai,'bi',bi,'s2',s2,'c2',c2
    c2 = c2 + s2
    IF(i0==i .OR. i0==(i+1)) THEN        
         CALL es_sgq_dgdnax(es_intg_sdgdnax1,i,ai,bi,i0,s2s)  
         !write(*,*) 'i0',i0,'i',i,'ai',ai,'bi',bi,'s2s',s2s,'c2',c2     
         c2 = c2 + s2s
    ENDIF
ENDDO
END SUBROUTINE

