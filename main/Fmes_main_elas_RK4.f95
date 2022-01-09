PROGRAM Main
!USE GQRUL_INT
USE mod_gaussqxw 
USE mod_elasbody
IMPLICIT NONE
Double PRECISION :: s,dt,sumdt,start,finish!,xt
Double PRECISION :: mine1=0d0
integer :: j,nstep,istep!,turnoff_step
!s is the checked point 

Call CPU_TIME(start)

s = dble(N)!dble(N/2+1)
nstep = 20000
sumdt = 0d0
dt = 5d-3

beta = 0.20d0
!! if gs = 1 and elongational flow == (z, -r/2) ?
!! then (omega here) == 6*(omega in D.B's papers) ?
gs = 1d0
omega = 48d-2 
!/* surface tension (as a function of surface coverage)
!! surftenMax = k * T * surfcovMax / gamma0
!! for example, k*T = 8.31 J/(K*mol) * 293 K, surfcovMax = 5.3 * 10^-6mol/m^2, 
!! gamma0 = 40 mN/m, then, surftenMax = 12.9/40 = 0.3225
surftenMax = 0.32d0 
!!eta defined as surfcovRef/surfcovMax
eta = 0.0d0
dilaViscosity = 0.0d0
!*/
!write(10,*) 'Ca = electrical stress / elastic stress, kappa = surface tension/ elastic stress'
!write(10,*) 'Ca = ',beta,'kappa = ',kappa!,' Initial gap = ',gap
!write(*,*) 'Ca = ',beta,'kappa = ',kappa!,' Initial gap = ',gap

write(10,*) 'Ca = electrical stress / gamma0, eta = surfcovRef/surfcovMax'
Write(10,*) 'surftenMax = k * T * surfcovMax / gamma0'
write(10,*) 'Ca = ',beta,'eta = ',eta, 'surftenMax =',surftenMax,'dilaViscosity =',dilaViscosity
write(10,*) 'N = ',N,' dt = ',dt
!,' Initial gap = ',gap
write(*,*) 'Ca = ',beta,'eta = ',eta,'surftenMax =',surftenMax,'dilaViscosity =',dilaViscosity
!,' Initial gap = ',gap
!! set up gauss-legendre quadrature & speical quadrature
call gq_4
call sgq_4
!call egq_10
write(10,*) 'QM=',QM,'qx, qw'
Do j = 1,QM
   write(10,*) qx(j), qw(j)
enddo
write(10,*) 'sqm=',sqm,'sqx, sqw'
Do j = 1,sqm
   write(10,*) sqx(j), sqw(j)
enddo


CALL initialContour
CALL recorder(0)
Do istep = 1, nstep
     !CALL fm_RK4_evolve(dt)
     CALL fm_EEmarch_normal(dt,istep)
     CALL recorder(istep)
ENDDO


100 FORMAT(ES25.15,A1,ES25.15)
101 FORMAT(I8,A1,ES25.15,A1,ES25.15)
200 FORMAT(ES25.15,A1,ES25.15,A1,ES25.15)

CALL CPU_TIME(finish)
print '("Time = ",f15.3," seconds.")',finish-start
END PROGRAM
