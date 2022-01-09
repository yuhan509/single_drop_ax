SUBROUTINE fm_EEmarch_normal(dt,istep)
USE mod_elasbody

     IMPLICIT NONE
     INTEGER,INTENT(IN) :: istep
     DOUBLE PRECISION,intent(IN) :: dt
!     DOUBLE PRECISION :: del
     INTEGER i
     
     CALL geom_nodes1d()
     CALL solveVelocity(istep)
        
     do i=1,N+1
!!/  march only in normal direction 
   !!elongational flow   
!        del = (nr(i)*(ur1(i) - omega * 0.5d0 * r1(i)) &
!        + nz(i)*(uz1(i) + omega * z1(i)))*dt
        !write(21,*) 'i=N+1,del',del,
!   !!no external flow
!        del = (nr(i)*ur1(i) + nz(i)*uz1(i))*dt
!
!        z1(i) = del*nz(i) + z1(i)
!        r1(i) = del*nr(i) + r1(i)
!!/

!        !/under elongational flow  
!               
!          r1(i) = (ur1(i) - omega * 0.5d0 * r1(i))*dt + r1(i)
!          z1(i) = (uz1(i) + omega * z1(i))*dt + z1(i)  
     
        !/no external flow
         r1(i) = ur1(i)*dt + r1(i)
         z1(i) = uz1(i)*dt + z1(i)
         
!        r2(i)= r1(i)
!        z2(i)= -z1(i)  
     enddo
     CALL cubic_para()
!        r1(1) = 0d0
!        r1(N+1) = 0d0
!        z1(N/2+1) = 0d0

!        r2(1)= r1(1)
!        r2(N+1) = r1(N+1)
        
END SUBROUTINE 


