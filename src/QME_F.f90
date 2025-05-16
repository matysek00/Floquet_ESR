module QME_F
Use declarations
CONTAINS
!
! The Matrix
!
      subroutine coeff_matrix (Ndim, omega, NF, NCF, G, Nmatrix, A, B, A_fast_left,A_fast_right, faster)
      implicit none
      integer,intent (in) :: Ndim, NF, NCF
      integer,intent (in) :: Nmatrix
      integer :: i1, i2, l, j, m, n, u, v, p
      complex (qc), intent (in) :: G (Ndim, Ndim, Ndim, Ndim, NF, 2)
      complex (qc), intent (inout) :: A (Nmatrix, Nmatrix)
      complex (qc), intent (out) :: B(Nmatrix)
      complex (qc) :: A_left (Nmatrix, Nmatrix), A_right (Nmatrix, Nmatrix)
!       complex (qc), intent (inout) :: A_fast(Nmatrix, Nmatrix) ! for a faster calculation (no frequency on the rates)
      complex (qc), intent (inout) :: A_fast_left(Nmatrix, Nmatrix),A_fast_right(Nmatrix, Nmatrix) ! for a faster calculation (no frequency on the rates)
      real (q) :: omega
      logical :: faster
        
    if ((faster).and.((A_fast_right(1,1)+A_fast_left(1,1)).ne.zero)) then
        A = zero
        A = A_fast_left+A_fast_right
        
        ! loop on i1
        do n=1,NF ! n=1 will correspond to the lowest floquet number: 
        ! if NF=5 then n=1 means a Floquet number of -2
        do i1=1,Nmatrix/NF 
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim

    A (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1)) = A (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1))&
    & + Delta(l,j)+ omega*(n-NCF-1)
        
        enddo
        enddo
        
    else
        ! Initialize NF=2*NCF+1
        A = zero
        A_left=zero
        A_right=zero
! loop on i1
        do n=1,NF ! n=1 will correspond to the lowest floquet number: 
        ! if NF=5 then n=1 means a Floquet number of -2
        do i1=1,Nmatrix/NF 
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim 
            
! loop on i2
            do m = max(1,n-2),min(NF,n+2)  ! 2 because of having a cosine square
            do i2=1,Nmatrix/NF 
            v=1+(i2-1)/ndim
            u=i2-(v-1)*ndim

!            
! corresponding term in the rate equation:
!            Gamma_{j,v,v,u:n-m} rho_{l,u,m}
              p = (n-m) + NCF + 1
              i2p= (l-1)*ndim + u
        A_left (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1)) = A_left (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1))&
        & +  ui*((G (j,v,v,u,p,2)))
        A_right (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1)) = A_right (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1))&
        & +  ui*((G (j,v,v,u,p,1)))

!             write(125,*)l,j,v,u,n,m,i1,i2,i2p,hartree*A(i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1))&
!                 &,(G (j,v,v,u,p,1)+G (j,v,v,u,p,2))*Hartree

! second term
!           Gamma_{l,v,v,u:m-n} rho_{u,j,m}
              p = (m-n) + NCF + 1
              i2p=(u-1)*ndim + j
         A_left (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1)) = A_left (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1))&
        & +  ui*CONJG((G (l,v,v,u,p,2)))
        A_right (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1)) = A_right (i1+ndim*ndim*(n-1),i2p+ndim*ndim*(m-1))&
        & +  ui*CONJG((G (l,v,v,u,p,1)))

! third term
!           Gamma_{u,j,l,v:m-n} rho_{v,u,m}
              p = (m-n) + NCF + 1
         A_left (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1)) = A_left (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1))&
        & -  ui*CONJG(G (u,j,l,v,p,2))
        A_right (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1)) = A_right (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1))&
        & -  ui*CONJG(G (u,j,l,v,p,1))

        
! fourth term
!           Gamma_{v,l,j,u:n-m} rho_{v,u,m}
              p = (n-m) + NCF + 1
         A_left (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1)) = A_left (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1))&
        & -  ui*(G (v,l,j,u,p,2))
        A_right (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1)) = A_right (i1+ndim*ndim*(n-1),i2+ndim*ndim*(m-1))&
        & -  ui*(G (v,l,j,u,p,1))

        enddo
        enddo
        

        A (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1)) = A (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1))&
        &+ Delta (l,j) + omega*(n-NCF-1)
        
        
        enddo
        enddo

        A=A+A_left+A_right
        A_fast_left = A_left
        A_fast_right= A_right
        
!         if(faster)then
!         ! loop on i1
!         do n=1,NF ! n=1 will correspond to the lowest floquet number:
!         ! if NF=5 then n=1 means a Floquet number of -2
!         do i1=1,Nmatrix/NF
!
!     A_fast (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1)) = A_fast (i1+ndim*ndim*(n-1),i1+ndim*ndim*(n-1))&
!     & - omega*(n-NCF-1) ! we extracted the diagonal in omega for possible faster calculation
!
!         enddo
!         enddo
!         endif
        
endif
        
        do n=1,NF
        
        A(1+ndim*ndim*(n-1),:)=zero
        
        do l=1,Ndim
        
        i2 = 1 + ndim*ndim*(n-1) + (l-1)*(Ndim+1)
        A(1+ndim*ndim*(n-1),i2)=one
!         print*,n,l,i2,1+ndim*ndim*(n-1)

        enddo
        enddo
        
        ! B definition entails detailed balance
        B = zero
        B (1+ndim*ndim*(NCF)) = one
        
               return
      end subroutine coeff_matrix 
!
! Rate computation
!
   subroutine rates (Ndim, NF, NCF, omega, gamma_R_0,A_R, gamma_L_0,A_L,phi,GammaC,&
   &Temperature,Spin_polarization_R, Spin_polarization_L, G, GC,lambda,orb,seHa)
     implicit none
! Input:
     integer :: Ndim, NF, NCF
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature, omega
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
! Output: the Rates called G (:,:,:,:) here
     complex (qc), intent (out) :: G (:,:,:,:,:,:) ! for QME_T
     complex (qc), intent (out) :: GC (:,:,:,:,:,:) ! for Transport_F
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR1, fermiL1, ufermiR1, ufermiL1
     complex (qc) :: fermiR2, fermiL2, ufermiR2, ufermiL2
     complex (qc) :: fermiR3, fermiL3, ufermiR3, ufermiL3
! Only used in this subroutine
     integer :: v, l, j, u, pn, lorb
!      real (q), intent (in) :: gamma_L_1
     integer, intent(in) :: orb
     real (q), intent (in) :: phi, seHa
!      real (q), intent (in) :: gamma_R_1
     complex (qc), intent (in) :: A_R, A_L,GammaC
     complex (qc) :: g0p_up, g0m_up, g1p_up, g1m_up
     complex (qc) :: g0p_dn, g0m_dn, g1p_dn, g1m_dn
     complex (qc) :: Lvluj, Ljulv

! Driving amplitude by electrode:
!      A_L=gamma_L_1/gamma_L_0
!      A_R=gamma_R_1/gamma_R_0
G=zero
GC=zero
!loop on states:
     do j=1,Ndim
     do u=1,Ndim

     if((lambda(u,j,1).eq.zero).and.(lambda(u,j,2).eq.zero).and.&
     &((lambda(j,u,1).eq.zero).and.(lambda(j,u,2).eq.zero)))then
     
     ! This if will avoid making the integral for rates that are zero always
     
     ! Right electrode
        ! WHAT already set to zero just skip
            G (:,:,j,u,:,1) = zero
            GC (:,:,j,u,:,1) = zero

     ! Left electrode
            G (:,:,j,u,:,2) = zero
            GC (:,:,j,u,:,2) = zero
!             write(*,*) lambda(u,j,1), lambda (u,j,2), lambda (j,u,1),  lambda (j,u,2)
        ! TODO: use continue to skip the rest of the loop no else here 
     else

! Green's functions integrals involving occupation factors
call ExtendedFermiIntegral(Delta(j,u),bias_R,Temperature,Cutoff,dimag(GammaC),N_int,fermiR1,WW,gau)
call ExtendedFermiIntegral(Delta(j,u),bias_L,Temperature,Cutoff,dimag(GammaC),N_int,fermiL1,WW,gau)
call ExtendedFermiIntegral(Delta(j,u)+omega,bias_R,Temperature,Cutoff,dimag(GammaC),N_int,fermiR2,WW,gau)
call ExtendedFermiIntegral(Delta(j,u)+omega,bias_L,Temperature,Cutoff,dimag(GammaC),N_int,fermiL2,WW,gau)
call ExtendedFermiIntegral(Delta(j,u)-omega,bias_R,Temperature,Cutoff,dimag(GammaC),N_int,fermiR3,WW,gau)
call ExtendedFermiIntegral(Delta(j,u)-omega,bias_L,Temperature,Cutoff,dimag(GammaC),N_int,fermiL3,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u),bias_R,Temperature,Cutoff,dble(GammaC),N_int,ufermiR1,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u),bias_L,Temperature,Cutoff,dble(GammaC),N_int,ufermiL1,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u)+omega,bias_R,Temperature,Cutoff,dble(GammaC),N_int,ufermiR2,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u)+omega,bias_L,Temperature,Cutoff,dble(GammaC),N_int,ufermiL2,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u)-omega,bias_R,Temperature,Cutoff,dble(GammaC),N_int,ufermiR3,WW,gau)
call ExtendeduFermiIntegral(Delta(j,u)-omega,bias_L,Temperature,Cutoff,dble(GammaC),N_int,ufermiL3,WW,gau)
!       AHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!!!
       
! FLOQUET n=0
            pn = 0 + NCF + 1
!       So the 5 harmonics are always calculated but everything above NCF and below -NCF
!       is stored outside of the matrix range
!       AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right electrode
        
        ! adding phi azimuthal angle and changing 0.5 in A*A terms to 0.25 (correct value)
        ! no phase shift in floquet number 0
        
      g0p_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(fermiR1+seHa*0.25*fermiR3*(abs(A_R))**2&
                +seHa*0.25*fermiR2*(abs(A_R))**2)
      g0m_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(ufermiR1+seHa*0.25*ufermiR3*(abs(A_R))**2&
                +seHa*0.25*ufermiR2*(abs(A_R))**2)
      g0p_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(fermiR1+seHa*0.25*fermiR3*(abs(A_R))**2&
                +seHa*0.25*fermiR2*(abs(A_R))**2)
      g0m_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(ufermiR1+seHa*0.25*ufermiR3*(abs(A_R))**2&
                +seHa*0.25*ufermiR2*(abs(A_R))**2)

! Left electrode

      g1p_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(fermiL1+seHa*0.25*fermiL3*(abs(A_L))**2&
                +seHa*0.25*fermiL2*(abs(A_L))**2)
      g1m_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(ufermiL1+seHa*0.25*ufermiL3*(abs(A_L))**2&
                +seHa*0.25*ufermiL2*(abs(A_L))**2)
      g1p_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(fermiL1+seHa*0.25*fermiL3*(abs(A_L))**2&
                +seHa*0.25*fermiL2*(abs(A_L))**2)
      g1m_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(ufermiL1+seHa*0.25*ufermiL3*(abs(A_L))**2&
                +seHa*0.25*ufermiL2*(abs(A_L))**2)
        
        do v=1,Ndim
        do l=1,Ndim
            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g0p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g0p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g0m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g0m_dn
        k=k+1
            enddo
            
! Right electrode
            G (v,l,j,u,pn,1) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,1) = 0.5*(Lvluj - Ljulv)

            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g1p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g1p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g1m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g1m_dn
        k=k+1
            enddo
            
! Left electrode
            G (v,l,j,u,pn,2) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,2) = 0.5*(Lvluj - Ljulv)
            
! write(124,*) hartree*G (v,l,j,u,pn,1),hartree*G (v,l,j,u,pn,2),v,l,j,u,pn

        enddo
        enddo

! FLOQUET n=-1
            pn = -1 + NCF + 1

! Right electrode

      g0p_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(0.5*fermiR1*A_R*zexp(ui*phi)+0.5*fermiR3*conjg(A_R)*zexp(ui*phi))
      g0m_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(0.5*ufermiR1*conjg(A_R)*zexp(ui*phi)+0.5*ufermiR2*A_R*zexp(ui*phi))
      g0p_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(0.5*fermiR1*A_R*zexp(ui*phi)+0.5*fermiR3*conjg(A_R)*zexp(ui*phi))
      g0m_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(0.5*ufermiR1*conjg(A_R)*zexp(ui*phi)+0.5*ufermiR2*A_R*zexp(ui*phi))

! Left electrode

      g1p_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(0.5*fermiL1*A_L*zexp(ui*phi)+0.5*fermiL3*conjg(A_L)*zexp(ui*phi))
      g1m_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(0.5*ufermiL1*conjg(A_L)*zexp(ui*phi)+0.5*ufermiL2*A_L*zexp(ui*phi))
      g1p_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(0.5*fermiL1*A_L*zexp(ui*phi)+0.5*fermiL3*conjg(A_L)*zexp(ui*phi))
      g1m_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(0.5*ufermiL1*conjg(A_L)*zexp(ui*phi)+0.5*ufermiL2*A_L*zexp(ui*phi))

        do v=1,Ndim
        do l=1,Ndim
           
           k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g0p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g0p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g0m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g0m_dn
        k=k+1
            enddo
            
! Right electrode
            G (v,l,j,u,pn,1) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,1) = 0.5*(Lvluj - Ljulv)

            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g1p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g1p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g1m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g1m_dn
        k=k+1
            enddo
            
! Left electrode
            G (v,l,j,u,pn,2) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,2) = 0.5*(Lvluj - Ljulv)
        enddo
        enddo

! FLOQUET n=+1
            pn = 1 + NCF + 1

! Right electrode

      g0p_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(0.5*fermiR1*A_R*zexp(-ui*phi)+0.5*fermiR2*conjg(A_R)*zexp(-ui*phi))
      g0m_up = gamma_R_0*(1+Spin_polarization_R)*0.5*(0.5*ufermiR1*conjg(A_R)*zexp(-ui*phi)+0.5*ufermiR3*A_R*zexp(-ui*phi))
      g0p_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(0.5*fermiR1*A_R*zexp(-ui*phi)+0.5*fermiR2*conjg(A_R)*zexp(-ui*phi))
      g0m_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*(0.5*ufermiR1*conjg(A_R)*zexp(-ui*phi)+0.5*ufermiR3*A_R*zexp(-ui*phi))

! Left electrode

      g1p_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(0.5*fermiL1*A_L*zexp(-ui*phi)+0.5*fermiL2*conjg(A_L)*zexp(-ui*phi))
      g1m_up = gamma_L_0*(1+Spin_polarization_L)*0.5*(0.5*ufermiL1*conjg(A_L)*zexp(-ui*phi)+0.5*ufermiL3*A_L*zexp(-ui*phi))
      g1p_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(0.5*fermiL1*A_L*zexp(-ui*phi)+0.5*fermiL2*conjg(A_L)*zexp(-ui*phi))
      g1m_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*(0.5*ufermiL1*conjg(A_L)*zexp(-ui*phi)+0.5*ufermiL3*A_L*zexp(-ui*phi))

        do v=1,Ndim
        do l=1,Ndim
        
            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g0p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g0p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g0m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g0m_dn
        k=k+1
            enddo
            
! Right electrode
            G (v,l,j,u,pn,1) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,1) = 0.5*(Lvluj - Ljulv)

            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g1p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g1p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g1m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g1m_dn
        k=k+1
            enddo
                    
! Left electrode
            G (v,l,j,u,pn,2) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,2) = 0.5*(Lvluj - Ljulv)
            
!         write(124,*) hartree*G (v,l,j,u,pn,1),hartree*G (v,l,j,u,pn,2),v,l,j,u,pn,A_L

        enddo
        enddo

! FLOQUET n=-2
            pn = -2 + NCF + 1

! Right electrode


! int((gau-1)/(1+gau))

      g0p_up = gamma_R_0*(1+Spin_polarization_R)*0.5*seHa*0.25*fermiR3*(abs(A_R)*zexp(ui*phi))**2
      g0m_up = gamma_R_0*(1+Spin_polarization_R)*0.5*seHa*0.25*ufermiR2*(abs(A_R)*zexp(ui*phi))**2
      g0p_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*seHa*0.25*fermiR3*(abs(A_R)*zexp(ui*phi))**2
      g0m_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*seHa*0.25*ufermiR2*(abs(A_R)*zexp(ui*phi))**2

! Left electrode

      g1p_up = gamma_L_0*(1+Spin_polarization_L)*0.5*seHa*0.25*fermiL3*(abs(A_L)*zexp(ui*phi))**2
      g1m_up = gamma_L_0*(1+Spin_polarization_L)*0.5*seHa*0.25*ufermiL2*(abs(A_L)*zexp(ui*phi))**2
      g1p_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*seHa*0.25*fermiL3*(abs(A_L)*zexp(ui*phi))**2
      g1m_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*seHa*0.25*ufermiL2*(abs(A_L)*zexp(ui*phi))**2

        do v=1,Ndim
        do l=1,Ndim
            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g0p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g0p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g0m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g0m_dn
        k=k+1
            enddo
            
! Right electrode
            G (v,l,j,u,pn,1) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,1) = 0.5*(Lvluj - Ljulv)

            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g1p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g1p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g1m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g1m_dn
        k=k+1
            enddo
! Left electrode
            G (v,l,j,u,pn,2) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,2) = 0.5*(Lvluj - Ljulv)
            
        enddo
        enddo

! FLOQUET n=+2
            pn = 2 + NCF + 1

! Right electrode

      g0p_up = gamma_R_0*(1+Spin_polarization_R)*0.5*seHa*0.25*fermiR2*(abs(A_R)*zexp(-ui*phi))**2
      g0m_up = gamma_R_0*(1+Spin_polarization_R)*0.5*seHa*0.25*ufermiR3*(abs(A_R)*zexp(-ui*phi))**2
      g0p_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*seHa*0.25*fermiR2*(abs(A_R)*zexp(-ui*phi))**2
      g0m_dn = gamma_R_0*(1-Spin_polarization_R)*0.5*seHa*0.25*ufermiR3*(abs(A_R)*zexp(-ui*phi))**2

! Left electrode

      g1p_up = gamma_L_0*(1+Spin_polarization_L)*0.5*seHa*0.25*fermiL2*(abs(A_L)*zexp(-ui*phi))**2
      g1m_up = gamma_L_0*(1+Spin_polarization_L)*0.5*seHa*0.25*ufermiL3*(abs(A_L)*zexp(-ui*phi))**2
      g1p_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*seHa*0.25*fermiL2*(abs(A_L)*zexp(-ui*phi))**2
      g1m_dn = gamma_L_0*(1-Spin_polarization_L)*0.5*seHa*0.25*ufermiL3*(abs(A_L)*zexp(-ui*phi))**2

        do v=1,Ndim
        do l=1,Ndim
            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g0p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g0p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g0m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g0m_dn
        k=k+1
            enddo
            
! Right electrode
            G (v,l,j,u,pn,1) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,1) = 0.5*(Lvluj - Ljulv)

            k=0
            Lvluj=zero
            Ljulv=zero
            
            do lorb=1,orb
            
        Lvluj = Lvluj + lambda (v,l,lorb+k)*conjg(lambda(u,j,lorb+k))*g1p_up+  &
                lambda (v,l,lorb+1+k)*conjg(lambda(u,j,lorb+1+k))*g1p_dn
                
        Ljulv = Ljulv + lambda (j,u,lorb+k)*conjg(lambda(l,v,lorb+k))*g1m_up+  &
                lambda (j,u,lorb+1+k)*conjg(lambda(l,v,lorb+1+k))*g1m_dn
        k=k+1
            enddo
! Left electrode
            G (v,l,j,u,pn,2) = 0.5*(Lvluj + Ljulv)
            GC (v,l,j,u,pn,2) = 0.5*(Lvluj - Ljulv)
        enddo
        enddo
    endif
        
      enddo
     enddo

      return
      
   end subroutine rates 
!
! Fermi occupation function
!
    function Fermi (e, T)
      implicit none
      real (q) :: Fermi, e, T
         if ( e > 0._q ) then
            Fermi = exp(-e/T)/(exp(-e/T)+1._q)
         else
           Fermi = 1._q/(exp(e/T)+1._q)
         endif
      return
    end function Fermi
!
! Calculation of energy integration of rates involving the Fermi function
! First type of integral
!
subroutine ExtendedFermiIntegral ( D, V, T, Cutoff, GammaC, N,fermiA,WW,gau)
      implicit none
      real (q) :: D, V, T, Cutoff,imG, GammaC
      real (q) :: e, step_e, WW,gau
      integer :: i, N
      complex (qc):: fermiA
!fermiA is Integral I11 of the Manual
! Trapeze-integration (the best among the better)
        imG=1._q
      step_e = 2*Cutoff/(N-1)
      e= -Cutoff+D
      fermiA=(gau*dexp(-0.5*((e-D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*0.5*Fermi (e-V, T) &
      &*(-ui*GammaC/((e-D)**2+GammaC**2)+(e-D)/((e-D)**2+imG*GammaC**2))
        
        ! factor (1./(WW*sqrt(2*pi_d))) made equal to 1
        
      do i = 2, N-1
      e= -Cutoff +D + (i-1)*step_e
      fermiA=fermiA+(gau*dexp(-0.5*((e-D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*Fermi (e-V, T) &
      &*(-ui*GammaC/((e-D)**2+GammaC**2)+(e-D)/((e-D)**2+imG*GammaC**2))
      enddo
      e = Cutoff +D
      fermiA=fermiA+(gau*dexp(-0.5*((e-D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*0.5*Fermi (e-V, T) &
      &*(-ui*GammaC/((e-D)**2+GammaC**2)+(e-D)/((e-D)**2+imG*GammaC**2))

      fermiA = step_e*ui*fermiA/pi_d
!       print*,fermiA,'f'
      return
      end subroutine ExtendedFermiIntegral
!
! Calculation of energy integration of rates involving 1-Fermi function
!
subroutine ExtendeduFermiIntegral (D,V,T,Cutoff,GammaC,N,ufermiA,WW,gau)
      implicit none
      real (q) :: D, V, T, Cutoff,imG, GammaC
      real (q) :: e, step_e, WW,gau,Ef
      integer :: i, N
      complex (qc):: ufermiA
!ufermiA is Integral I21 of the Manual
! Trapeze-integration (the best among the better)
! we add guassian
        imG=1._q
      step_e = 2*Cutoff/(N-1)
      e= -Cutoff-D
      ufermiA=(gau*dexp(-0.5*((e+D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*0.5*(1-Fermi (e-V, T)) &
      &*(ui*GammaC/((e+D)**2+GammaC**2)+(e+D)/((e+D)**2+imG*GammaC**2))

      do i = 2, N-1
      e= -Cutoff -D + (i-1)*step_e
      ufermiA=ufermiA+(gau*dexp(-0.5*((e+D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*(1-Fermi (e-V, T))&
      &*(ui*GammaC/ ((e+D)**2+GammaC**2)+(e+D)/((e+D)**2+imG*GammaC**2))
      enddo
      e = Cutoff -D
      ufermiA=ufermiA+(gau*dexp(-0.5*((e+D)/WW)**2) &
      &-int((gau-1)/(1+gau)))*0.5*(1-Fermi (e-V, T)) &
      &*(ui*GammaC/((e+D)**2+GammaC**2)+(e+D)/((e+D)**2+imG*GammaC**2))
! print*,int((gau-1)/(1+gau))
      ufermiA = -step_e*ui*ufermiA/pi_d
! print*,ufermiA,'1-f'
      return
      end subroutine ExtendeduFermiIntegral

end module QME_F
