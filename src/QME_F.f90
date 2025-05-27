module QME_F
Use declarations
CONTAINS
!
! Rate computation
!
   subroutine rates (Ndim, orb, frequency, gamma_0, lambda, Spin_polarization,&
        NCF, p_max, Adrive, phi, Bdrive, GammaC, bias, Delta, Cutoff, Temperature,& 
        seHa, WW, gau, N_int, GA, GCA)
     implicit none

     integer :: Ndim, NCF, p_max, orb, N_int
     complex (qc), intent (in), dimension(Ndim,Ndim,orb) :: lambda 
     real (q), intent (in), dimension(Ndim, Ndim) :: Delta
     real (q), intent (in) :: gamma_0, Temperature, frequency, Spin_polarization
     real (q), intent (in) :: Cutoff, WW, gau, phi, seHa, bias, Bdrive
     complex (qc), intent (in) :: GammaC, Adrive

     !complex (qc), intent (out), dimension(Ndim, Ndim, Ndim, Ndim, 2*NCF+1) :: GA, GCA 
     complex (qc), intent (out) :: GA (:,:,:,:,:), GCA (:,:,:,:,:)
!       G_Alpha and GC_Alpha 
     complex (qc), dimension(2*p_max-1) :: fermi, ufermi     
     complex (q), dimension(2*p_max-1) :: Kbess
     complex (qc) :: g_up, g_dn
     complex (qc), dimension(Ndim, Ndim) :: Lvluj, Ljulv
     complex (qc) :: bessel_contribution, ubessel_contribution

     integer :: j, u, nfour, n_index
        
        g_up = 0.5*gamma_0*(1+Spin_polarization)
        g_dn = 0.5*gamma_0*(1-Spin_polarization)

!       Calculate Contribution of Bessel functions
!       K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
        Kbess = Bessel_K(Bdrive/frequency, Adrive, p_max)

!       loop on states:
        level_j: do j=1,Ndim
        level_u: do u=1,Ndim

!       Skip lambda is zero
        if((lambda(u,j,1).eq.zero).and.(lambda(u,j,2).eq.zero).and.&
                &((lambda(j,u,1).eq.zero).and.(lambda(j,u,2).eq.zero)))then
                GA (:,:, j,u,:) = zero
                GCA (:,:, j,u,:) = zero
                cycle 
                ! TODO: should this include other spin channels?
!               Let the record show that this cycle avoided a 400 line if statement 
!               when it was originally added.
        end if

!       Green's functions integrals involving occupation factors
        call ExtendedFermiIntegral(Delta(j,u), frequency, bias, p_max-1, Temperature, Cutoff, &
                        GammaC, N_int, WW, gau, fermi, ufermi)

        call Orbital_overlaps(lambda, j, u, orb, g_up, g_dn, Ndim, Lvluj, Ljulv)
        fourier_component: do nfour = -NCF, NCF
                n_index = nfour + NCF + 1

!               contribution of bessel functions
!               bessel_cont  = sum_p K*_{p-n} K_p  I(p)                              
!               ubessel_cont = sum_p K*_p K_{p+n}  uI(p) 
                call compute_bessel_contribution(Kbess, fermi, ufermi, p_max, nfour, &
                              bessel_contribution, ubessel_contribution)
                
                ! TODO: ADD PHASE
                GA  (:,:,j,u,n_index) = 0.5*(Lvluj*bessel_contribution + Ljulv*ubessel_contribution)
                GCA (:,:,j,u,n_index) = 0.5*(Lvluj*bessel_contribution - Ljulv*ubessel_contribution)
                
        enddo fourier_component
        enddo level_u
        enddo level_j
      return
      
   end subroutine rates 

!
!       Fermi occupation function with limitors to avoid underflow/overflows
!
        function FermiDist (e) result (eFermi)
        implicit none
        real (q) :: eFermi
        real (q) :: eFermi
        real (q), intent(in) :: e
        
        if ( e > 0._q ) then
                eFermi = exp(-e)/(exp(-e)+1._q)
             else
               eFermi = 1._q/(exp(e)+1._q)
        endif
      end function FermiDist

      
        subroutine ExtendedFermiIntegral (D, frequency, bias, p_max, Temperature, Cutoff, GammaC,&
                 N, WW, gau, fermi, ufermi)
        implicit none
        
        integer, intent(in) :: N, p_max
        real (q), intent(in) :: D, bias, Temperature, Cutoff, frequency, WW, gau
        complex (qc), intent(in) :: GammaC
        
        complex (qc), intent(out), dimension(2*p_max+1) :: fermi, ufermi

        real (q) :: e, step_e, gaushift, WWsq, imG, gausian
        real (q) :: rGammaC, iGammaC, rGammaCsq, iGammaCsq
        real (q), dimension(N) :: f, uf        
        integer :: i, p, p_ind

!       The FermiIntegrand function is used to calculate the fermi integral
!       The FermiIntegrand function is used to calculate the fermi integral
!       G(e) = gau*dexp(-0.5*esq/WWsq) - gaushift
!       A(e)  = (e/(esq + imG*Re(GammaC)**2) - ui*Re(GammaC)/(esq + Re(GammaC)**2))
!       uA(e) = (e/(esq + imG*Im(GammaC)**2) + ui*Im(GammaC)/(esq + Im(GammaC)**2))
!       fermi  += int f(e) * A(e) * G(e) de
!       ufermi += int (1-f(e)) * uA(e) * G(e) de
!       f(E,V) = \frac{1}{\exp(\beta (E-V)) + 1} Fermi distribution with V as fermi level

!       TODO: think about how to organize this, probably should splie fermi and ufermi
!       calculate all fermi contributions they are the same for all integrals
        step_e = 2*Cutoff/(N-1)/Temperature ! rescaling to units T = 1 
        e= -(Cutoff - D + bias)/Temperature
        fstep: do i = 1, N
             e = e + step_e
             f(i) = FermiDist (e)
        enddo fstep

        e = -(Cutoff + D + bias)/Temperature
        Ufstep: do i = 1, N
             e = e + step_e
             uf(i) = 1 - FermiDist (e)
        enddo ufstep
        
        step_e = 2*Cutoff/(N-1) ! now in atomic units
        
        imG = 1._q
        gaushift = int((gau-1)/(1+gau))
        
        rGammaC = dble (GammaC)
        iGammaC = dimag(GammaC)
        rGammaCsq = rGammaC**2
        iGammaCsq = iGammaC**2
        WWsq = WW**2

        ploop: do p = -p_max, p_max
                p_ind = p+p_max+1
                p_ind = p+p_max+1
   
                e = -Cutoff + p*frequency
                fermi(p_ind)  = .5*FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        -rGammaC, rGammaCsq, f(1))
                ufermi(p_ind) = .5*FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        iGammaC, iGammaCsq, uf(1))

                istep: do i = 2, N - 1
                e = e + step_e
                fermi(p_ind)  = fermi(p_ind)  + FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        -rGammaC, rGammaCsq, f(i))
                ufermi(p_ind) = ufermi(p_ind) + FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        iGammaC, iGammaCsq, uf(i))
                
                enddo istep
             
             e = Cutoff + p*frequency
             fermi(p_ind)  = fermi(p_ind) +  .5*FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        -rGammaC, rGammaCsq, f(N))
             ufermi(p_ind) = ufermi(p_ind) + .5*FermiIntStep(e, WWsq, gaushift, gau, imG,&
                        iGammaC, iGammaCsq, uf(N))
   
        enddo ploop
        
        fermi  =  step_e*ui*fermi/pi_d
        ufermi = -step_e*ui*ufermi/pi_d

        return
        end subroutine ExtendedFermiIntegral


        function FermiIntStep (e, WWsq, gaushift, gau, imG, riGammaC, riGammaCsq, fd) result (fermiInt)
        implicit none
        real (q), intent(in) :: e, WWsq, gaushift, gau, imG, riGammaC, riGammaCsq, fd
        real (q) :: esq, gausian
        complex (qc) :: fermiInt, integrand
                
                esq = e**2
                gausian = gau*dexp(-0.5*esq/WWsq) - gaushift
                integrand  = e/(esq + imG*riGammaCsq) + ui*riGammaC/(esq + riGammaCsq)
                fermiInt = fd * integrand * gausian
        return
        end function FermiIntegrand
        end function FermiIntegrand


        function Bessel_K(z, Adrive, p_max) result (Kbess)
        implicit none
        real (q), intent (in) :: z
        complex(qc), intent(in) :: Adrive
        integer, intent (in) :: p_max
        complex (qc), dimension(2*p_max-1) :: Kbess
        complex(qc), dimension(2*p_max+1) :: Jbess
        integer :: p
   
                Jbess(p_max+1:) = Bessel_JN(0, p_max, z)
             
                ! J(-p) = (-1)**p J(p) for p = 0,1,2,...
                negative_bessel : do p = 0, p_max -1
                        Jbess(p+1) = ((-1)**(p_max-p))*Jbess(2*p_max+1-p)
                enddo negative_bessel
             
                ! K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
                Kbess = Jbess(2:2*p_max) + 0.5 * Adrive * (Jbess(1:2*p_max-1) + Jbess(3:2+p_max+1))
        
        return
        end function Bessel_K

        subroutine Orbital_overlaps (lambda, j, u, orb, g_up, g_dn, Ndim, overlapvluj, overlapjulv)
        implicit none
        integer, intent(in):: Ndim, j, u, orb
        complex (qc), dimension(Ndim, Ndim, 2), intent (in) :: lambda 
        complex (qc), intent (in) :: g_up, g_dn
        complex (qc), dimension(Ndim, Ndim), intent (out) :: overlapvluj, overlapjulv
        integer :: v, l, lorb 
                
                overlapvluj = zero
                overlapjulv = zero
                
                level_v: do v=1,Ndim
                level_l: do l=1, Ndim
                        
                        orbital: do lorb = 1, orb, 2 
                        overlapvluj(v,l) = lambda(v,l,lorb)*conjg(lambda(u,j,lorb))*g_up+&
                                lambda (v,l,lorb+1)*conjg(lambda(u,j,lorb+1))*g_dn
                        overlapjulv(v,l) = lambda(j,u,lorb)*conjg(lambda(l,v,lorb))*g_up+&
                                lambda (j,u,lorb+1)*conjg(lambda(l,v,lorb+1))*g_dn
                        enddo orbital

                enddo level_l
                enddo level_v
   
        return
        end subroutine Orbital_overlaps


        subroutine compute_bessel_contribution(K, fermi, ufermi, p_max, nfour, result_bessel, result_ubessel)
        integer, intent(in) :: p_max, nfour
        complex (qc), intent(in), dimension(2*p_max-1) :: K
        complex (qc), intent(in), dimension(2*p_max-1) :: fermi, ufermi
        complex (qc), intent(out) :: result_bessel, result_ubessel
   
             integer :: p
             result_bessel = 0.0_qc 
             result_ubessel = 0.0_qc
             
             ! sum_p K*_{p-n} K_p  I(p)
             ! sum_p K*_p K_{p+n}  uI(p)
             ! assumes that for |p|>p_max K_p = 0 
             ! for n<0 goes from p=1-p_max to p=p_max-1-n
             ! for n>0 goes from p=1-p_max to p=p_max-1
             ! thus p, p+n, p-n are all in the range 1-p_max to p_max-1
             bessel: do p = max(1, 1-nfour), min(2*p_max-1, 2*p_max-nfour-1)
                  result_bessel  = result_bessel  + conjg(K(p)) * K(p+nfour) * fermi(p+nfour)
                  result_ubessel = result_ubessel + conjg(K(p)) * K(p+nfour) * ufermi(p)
             enddo bessel
             
             return
        end subroutine compute_bessel_contribution

end module QME_F