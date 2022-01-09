Module mod_elasbody
use mod_gaussqxw
IMPLICIT NONE
SAVE
INTEGER,PARAMETER :: N = 64
!! element number N, nodes number = N+1
DOUBLE PRECISION,DIMENSION(N+1) :: r1,z1,dphidn1,ps1
DOUBLE PRECISION,DIMENSION(N+1) :: uz1,ur1,nr,nz
DOUBLE PRECISION,DIMENSION(3,N) :: ar1,az1
DOUBLE PRECISION,DIMENSION(N+1,3,N) :: acar_na, acar_cl
DOUBLE PRECISION :: PI, beta, gap, gs, omega, kappa, dilaViscosity
DOUBLE PRECISION :: surftenMax, eta

!DOUBLE PRECISION,DIMENSION(qm,N) :: r2o,nr2o,z2o,nz2o,dlds2o,div_n2o,&
!sforce_r2o,sforce_z2o,maxwstress_n2o
!DOUBLE PRECISION,DIMENSION(qm,N) :: es_r1o,es_nr1o,es_z1o,es_nz1o,es_dlds1o
!DOUBLE PRECISION,DIMENSION(qm,N) :: es_r2o,es_nr2o,es_z2o,es_nz2o,es_dlds2o

DOUBLE PRECISION,DIMENSION(qm,N) :: ogd11,ogd22,ogu11,ogu22 &
,odgd11ds,odgd22ds,odgu11ds,odgu22ds

DOUBLE PRECISION,DIMENSION(2*sqm,N) :: osgd11,osgd22,osgu11,osgu22 &
,osdgd11ds,osdgd22ds,osdgu11ds,osdgu22ds

DOUBLE PRECISION,DIMENSION(qm,N) :: s1o,r1o,nr1o,z1o,nz1o,dlds1o,&
sforce_r1o,sforce_z1o,tr1o,tz1o,gu11o1,gu22o1,dgu22dso1,d2zds2o1,d2rds2o1,&
meancurv1o,partofcovarintderiv1o,Mzz1o,Mzr1o,Mrz1o,Mrr1o,&
normal_stress1o,tangent_stress1o

DOUBLE PRECISION,DIMENSION(2*sqm,N) :: sforce_r1o_s,sforce_z1o_s,dlds1o_s,&
s1o_s,r1o_s,z1o_s,tr1o_s,tz1o_s,d2rds2o1_s,d2zds2o1_s,gu22o1_s,dgu22dso1_s,&
nr1o_s,nz1o_s,partofcovarintderiv1o_s,meancurv1o_s,Mzz1o_s,Mzr1o_s,Mrz1o_s,&
Mrr1o_s,normal_stress1o_s,tangent_stress1o_s

DOUBLE PRECISION,DIMENSION(N+1,qm,N) :: car_na1o,dcards_na1o,d2cards2_na1o,&
car_cl1o,dcards_cl1o,d2cards2_cl1o,delForce_z_coefUzj1o,delForce_r_coefUzj1o,&
delForce_r_coefUrj1o,delForce_z_coefUrj1o

DOUBLE PRECISION,DIMENSION(N+1,2*sqm,N) :: car_na1o_s,dcards_na1o_s,d2cards2_na1o_s,&
car_cl1o_s,dcards_cl1o_s,d2cards2_cl1o_s,delForce_z_coefUzj1o_s,delForce_r_coefUzj1o_s,&
delForce_r_coefUrj1o_s,delForce_z_coefUrj1o_s


CONTAINS 

SUBROUTINE sphere_nodes_rz()
!! use N,r1,z1,PI
IMPLICIT NONE
INTEGER :: j
!DOUBLE PRECISION :: displ
PI = 4d0*datan(1d0)
!displ = 1d0+gap/2d0
!write(13,*) 'separation=', gap
Do j=1,N/2
    !!Ellipse: z*z/4+r*r=1
    r1(j)=dsin(PI*(j-1)/N)
    r1(N+2-j)=r1(j) !!mirror
    z1(j)=-dcos(PI*(j-1)/N)
    z1(N+2-j)=-z1(j) !!mirror
END DO
r1(N+1) = 0d0
z1(N/2+1) = 0d0
r1(N/2+1) = 1d0    
!r2(N+1) = r1(N+1)
END SUBROUTINE 



SUBROUTINE chordal_para()
     IMPLICIT NONE
     INTEGER j
     DOUBLE PRECISION :: h
     DOUBLE PRECISION,DIMENSION(N+1) :: card
     
     !! chordal parametrization for nodes
     ps1(1)=0d0
     !ps2(1)=0d0
     do j = 1,N
        h = dsqrt((r1(j+1)-r1(j))**2+(z1(j+1)-z1(j))**2)
        ps1(j+1)= ps1(j)+h
        !ps2(j+1)= ps1(j+1)
     enddo
    
     !!generate cubic polynomial parameters in each elements
     call chordal_cubics_cl(N,ps1,z1,az1)
     !call chordal_cubics_cl(N,ps2,z2,az2)
     call chordal_cubics_na(N,ps1,r1,ar1)
     !call chordal_cubics_na(N,ps2,r2,ar2)
     !!generate cubic-spline cardinal functions
     Do j=1,N+1
         card = 0d0
         card(j)=1d0
         Call chordal_cubics_cl(N,ps1,card,acar_cl(j,:,:))
         Call chordal_cubics_na(N,ps1,card,acar_na(j,:,:))
     ENDDO
     !write(*,*) 'acar(1,1,1)=',acar(1,1,1) 
     !write(*,*) 'acar(N+1,1,N)=',acar(N+1,1,N)
     !write(*,*) 'close to 0 ',acar(1,2,N/2)+acar(2,2,N/2)+acar(3,2,N/2)
     !write(*,*) 'card(N+1)',card
     !dont remove otherwise wierd result given
END SUBROUTINE


SUBROUTINE arc_para()
     IMPLICIT NONE
     INTEGER j
     DOUBLE PRECISION :: h
     DOUBLE PRECISION,DIMENSION(N+1) :: card
     

!!!   para from -pi/2 to pi/2  
!     ps1(1)= - PI * 0.5d0
!     h = PI/N
!     do j = 1,N/2
!        ps1(j+1) = ps1(j)+h
!        ps1(N+2-j)= - ps1(j)
!        !ps2(j+1)= ps1(j+1)
!     enddo
!     ps1(N/2+1) = 0d0
!!!

!!   para from 0 to pi
     ps1(1)= 0d0
     h = PI/N
     do j = 1,N
        ps1(j+1) = ps1(j) + h
     enddo
!!


     !!generate cubic polynomial parameters in each elements
     call chordal_cubics_cl(N,ps1,z1,az1)
     !call chordal_cubics_cl(N,ps2,z2,az2)
     call chordal_cubics_na(N,ps1,r1,ar1)
     !call chordal_cubics_na(N,ps2,r2,ar2)
     !!generate cubic-spline cardinal functions
     Do j=1,N+1
         card = 0d0
         card(j)=1d0
         Call chordal_cubics_cl(N,ps1,card,acar_cl(j,:,:))
         Call chordal_cubics_na(N,ps1,card,acar_na(j,:,:))
     ENDDO
     !write(*,*) 'acar(1,1,1)=',acar(1,1,1) 
     !write(*,*) 'acar(N+1,1,N)=',acar(N+1,1,N)
     !write(*,*) 'close to 0 ',acar(1,2,N/2)+acar(2,2,N/2)+acar(3,2,N/2)
     !write(*,*) 'card(N+1)',card
     !dont remove otherwise wierd result given
     
END SUBROUTINE


SUBROUTINE cubic_para()
     IMPLICIT NONE
!     DOUBLE PRECISION :: h
!     DOUBLE PRECISION,DIMENSION(N+1) :: card
     
!          !! chordal parametrization for nodes
!     ps1(1)=0d0
!     !ps2(1)=0d0
!     do j = 1,N
!        h = dsqrt((r1(j+1)-r1(j))**2+(z1(j+1)-z1(j))**2)
!        ps1(j+1)= ps1(j)+h
!        !ps2(j+1)= ps1(j+1)
!     enddo
!     
     !!generate cubic polynomial parameters in each elements
     call chordal_cubics_cl(N,ps1,z1,az1)
     !call chordal_cubics_cl(N,ps2,z2,az2)
     call chordal_cubics_na(N,ps1,r1,ar1)
     !call chordal_cubics_na(N,ps2,r2,ar2)
     !!generate cubic-spline cardinal functions
!!     Do j=1,N+1
!!         card = 0d0
!!         card(j)=1d0
!!         Call chordal_cubics_cl(N,ps1,card,acar(j,:,:))
!!     ENDDO
!     write(*,*) 'acar(1,1,1)=',acar(1,1,1) 
!     write(*,*) 'acar(N+1,1,N)=',acar(N+1,1,N)
!     write(*,*) 'close to 0 ',acar(1,2,N/2)+acar(2,2,N/2)+acar(3,2,N/2)
!     write(*,*) 'card(N+1)',card
     !dont remove otherwise wierd result given
END SUBROUTINE


SUBROUTINE geo1_rest_metrics_arrayo_cubs_gq()

    IMPLICIT NONE
    DOUBLE PRECISION :: s,rs,nrs,zs,nzs,dlds,div_n
    DOUBLE PRECISION :: temp,temp1,&
    dzds,drds,d2zds,d2rds,drdz,d2rdz
    DOUBLE PRECISION :: ps1i,ar11i,ar12i,ar13i,az11i,az12i,az13i,r1i,z1i
    !! derivatives
    DOUBLE PRECISION :: drdl,dzdl,d2lds,d2rdl,d2zdl 
    !! metric tensors for deformed state,g12 and g21 are zero
    !! u means upper index, contravariant; d means lower index, covariant
    DOUBLE PRECISION :: gu11,gu22,gd11,gd22
    DOUBLE PRECISION :: dgu11ds,dgu22ds,dgd11ds,dgd22ds
    INTEGER :: i,k    !!!  kth abscissa points,ith element;

    Do i=1,N 
         ar11i = ar1(1,i)
         ar12i = ar1(2,i) 
         ar13i = ar1(3,i)
         az11i = az1(1,i)
         az12i = az1(2,i)
         az13i = az1(3,i)
         r1i = r1(i)
         z1i = z1(i)
         ps1i = ps1(i)
         
         Do k=1,qm
             s = s1o(k,i)
             temp = s-ps1i
             rs=ar11i*temp**3+ar12i*temp**2+ar13i*temp+r1i
             zs=az11i*temp**3+az12i*temp**2+az13i*temp+z1i
             !first order derivatives
             drds=3d0*ar11i*temp**2+2d0*ar12i*temp+ar13i
             dzds=3d0*az11i*temp**2+2d0*az12i*temp+az13i
             dlds = dsqrt(drds**2+dzds**2)
             
             drdl = drds / dlds
             dzdl = dzds / dlds
             !second order derivatives
             d2rds = 6d0*ar11i*temp+2d0*ar12i
             d2zds = 6d0*az11i*temp+2d0*az12i
             
             d2lds = (drds * d2rds + dzds * d2zds) / dlds         
             d2rdl = d2rds / (dlds*dlds) - drds * d2lds / (dlds**3)
             d2zdl = d2zds / (dlds*dlds) - dzds * d2lds / (dlds**3)        
      
             !normal component,pointing outward
             nrs= dzds/dlds
             nzs= -drds/dlds       
             !For axisymetric curvature
             drdz = drds/dzds
              !second derivative of r in z
             d2rdz = d2rds/(dzds*dzds)-d2zds*drds/(dzds**3)        
             temp1 = dsqrt(1d0+drdz*drdz)
             ! div of unit normal vector
             ! equal to twice the mean curvature
             div_n = 1d0/(rs*temp1)-1d0/temp1**3*d2rdz
             div_n = sign(div_n,drds*d2zds-dzds*d2rds)
             
             ! metric tensors
             gu11 = 1d0 / (rs*rs)
             gd11 = rs * rs
             gu22 = 1d0 / (drds**2 + dzds**2)        
             gd22 = drds**2 + dzds**2
             
             !dervivative of metric tensors
             dgd11ds = 2d0 * rs * drds
             dgd22ds = 2d0 * (drds * d2rds + dzds * d2zds)
             dgu11ds = -2d0 * drds / rs**3
             dgu22ds = -2d0 * (drds * d2rds + dzds * d2zds) / (drds**2 + dzds**2)**2 
             
             !!  store the metric tensors for the rest state        
             ogd11(k,i) = gd11
             ogd22(k,i) = gd22
             ogu11(k,i) = gu11
             ogu22(k,i) = gu22
             odgd11ds(k,i) = dgd11ds
             odgd22ds(k,i) = dgd22ds
             odgu11ds(k,i) = dgu11ds
             odgu22ds(k,i) = dgu22ds
             
             ! write(70,*) 'i = ', i,' k = ', k,'ogd11 = ', ogd11(k,i)  
     
         enddo !! for k 
     enddo !! for i

END SUBROUTINE


SUBROUTINE geom_special_rest_metrics_arrayo_cubs_gq()

    IMPLICIT NONE
    DOUBLE PRECISION :: s,rs,nrs,zs,nzs,dlds,div_n
    DOUBLE PRECISION :: temp,temp1,&
    dzds,drds,d2zds,d2rds,drdz,d2rdz
    DOUBLE PRECISION :: ps1i,ar11i,ar12i,ar13i,az11i,az12i,az13i,r1i,z1i
    !! derivatives
    DOUBLE PRECISION :: drdl,dzdl,d2lds,d2rdl,d2zdl 
    !! metric tensors for deformed state,g12 and g21 are zero
    !! u means upper index, contravariant; d means lower index, covariant
    DOUBLE PRECISION :: gu11,gu22,gd11,gd22
    DOUBLE PRECISION :: dgu11ds,dgu22ds,dgd11ds,dgd22ds
    INTEGER :: i,k    !!!  kth special abscissa points,ith element;

    Do i=1,N 
         ar11i = ar1(1,i)
         ar12i = ar1(2,i)
         ar13i = ar1(3,i)
         az11i = az1(1,i)
         az12i = az1(2,i)
         az13i = az1(3,i)
         r1i = r1(i)
         z1i = z1(i)
         ps1i = ps1(i)
         
         Do k=1,2*sqm
             s = s1o_s(k,i)                       
             temp = s-ps1i    
      
             rs=ar11i*temp**3+ar12i*temp**2+ar13i*temp+r1i
             zs=az11i*temp**3+az12i*temp**2+az13i*temp+z1i
             !first order derivatives
             drds=3d0*ar11i*temp**2+2d0*ar12i*temp+ar13i
             dzds=3d0*az11i*temp**2+2d0*az12i*temp+az13i
             dlds = dsqrt(drds**2+dzds**2)
             
             drdl = drds / dlds
             dzdl = dzds / dlds
             !second order derivatives
             d2rds = 6d0*ar11i*temp+2d0*ar12i
             d2zds = 6d0*az11i*temp+2d0*az12i
             
             d2lds = (drds * d2rds + dzds * d2zds) / dlds         
             d2rdl = d2rds / (dlds*dlds) - drds * d2lds / (dlds**3)
             d2zdl = d2zds / (dlds*dlds) - dzds * d2lds / (dlds**3)        
      
             !normal component,pointing outward
             nrs= dzds/dlds
             nzs= -drds/dlds       
             !For axisymetric curvature
             drdz = drds/dzds
              !second derivative of r in z
             d2rdz = d2rds/(dzds*dzds)-d2zds*drds/(dzds**3)        
             temp1 = dsqrt(1d0+drdz*drdz)
             ! div of unit normal vector
             ! equal to twice the mean curvature
             div_n = 1d0/(rs*temp1)-1d0/temp1**3*d2rdz
             div_n = sign(div_n,drds*d2zds-dzds*d2rds)
             
             ! metric tensors
             gu11 = 1d0 / (rs*rs)
             gd11 = rs * rs
             gu22 = 1d0 / (drds**2 + dzds**2)        
             gd22 = drds**2 + dzds**2
             
             !dervivative of metric tensors
             dgd11ds = 2d0 * rs * drds
             dgd22ds = 2d0 * (drds * d2rds + dzds * d2zds)
             dgu11ds = -2d0 * drds / rs**3
             dgu22ds = -2d0 * (drds * d2rds + dzds * d2zds) / (drds**2 + dzds**2)**2 
             
             !!  store the metric tensors for the rest state        
             osgd11(k,i) = gd11
             osgd22(k,i) = gd22
             osgu11(k,i) = gu11
             osgu22(k,i) = gu22
             osdgd11ds(k,i) = dgd11ds
             osdgd22ds(k,i) = dgd22ds
             osdgu11ds(k,i) = dgu11ds
             osdgu22ds(k,i) = dgu22ds
             !write(71,*) 'i = ', i,' k = ', k,'osdgu22dl(k,i) = ', osdgu22dl(k,i) 
!             write(72,*) 'i=',i,' k=',k,' gu11=',gu11,' gu22=',gu22,' gd11=',gd11,' gd22=',&
!gd22,' dgu11ds=',dgu11ds,' dgu22ds=',dgu22ds,' dgd11ds=',dgd11ds,' dgd22ds=',dgd22ds
         enddo !! for k 
     enddo !! for i

END SUBROUTINE


SUBROUTINE geom_nodes1d()
     !! deliver nr(i),nz(i) for nodes, used in marching
     IMPLICIT NONE
     INTEGER i
     DOUBLE PRECISION drds,dzds,dlds
     Do i=2,N
     drds = ar1(3,i)
     dzds = az1(3,i)
     dlds = dsqrt(drds**2d0+dzds**2d0)
     nr(i)= dzds/dlds
     nz(i)= -drds/dlds !!pointing outward from drop
     Enddo
     !!! 1st and (N+1)th node:
     nr(1)= 0d0
     nz(1)= -1d0
     nr(N+1)= 0d0
     nz(N+1)= 1d0

END SUBROUTINE 



SUBROUTINE fm_dmarch(dt)
     IMPLICIT NONE
     DOUBLE PRECISION,intent(IN) :: dt
     !DOUBLE PRECISION :: del
     INTEGER i
     
     call geom_nodes1d
     do i=1,N+1
!/  march only in normal direction     
!        del = (nr(i)*ur1(i)+nz(i)*uz1(i))*dt
!        !write(21,*) 'i=N+1,del',del,
!        z1(i) = del*nz(i) + z1(i)
!        r1(i) = del*nr(i) + r1(i)
!/
        !/under elongational flow  
!         write(75,*) 'i=',i,'ur=',ur1(i),'uz=',uz1(i) 
               
          r1(i) = (ur1(i) - omega * r1(i))*dt + r1(i)
          z1(i) = (uz1(i) + omega * 2d0 * z1(i))*dt + z1(i)  
     
!        !/no external flow
!         r1(i) = ur1(i)*dt + r1(i)
!         z1(i) = uz1(i)*dt + z1(i)
         
!        r2(i)= r1(i)
!        z2(i)= -z1(i)  
     enddo
        r1(1) = 0d0
        r1(N+1) = 0d0
        z1(N/2+1) = 0d0
!        r2(1)= r1(1)
!        r2(N+1) = r1(N+1)
        
END SUBROUTINE 

end module mod_elasbody
