module H_QD
Use declarations
Use OpenFiles
Use algebra
Use KrondProd
implicit none

CONTAINS

   subroutine Hamiltonian (N,lambda,Delta,Eigenvalues,H,Ss,spinX,spinY,spinZ,&
   &spin2_T,runs,Nm,hx,hy,hz,ndim,&
   &redimension,orb,bias_R, bias_L, Cutoff)!&
!    &,spin_first)
   implicit none
   integer, intent(out) :: N, Nm,ndim
   integer, intent(in) :: orb
   real(q), intent(in) :: bias_R, bias_L, Cutoff
   real(q) :: bias
   real (q), allocatable, intent (out) :: hx (:), hy(:), hz(:)
   real (q), allocatable, intent (out) ::  Eigenvalues (:), Delta (:,:)
   complex (qc), intent (out), allocatable :: lambda (:,:,:), H(:,:), Ss (:,:,:,:)
   complex (qc),intent(out), allocatable :: spinX(:,:,:),spinY(:,:,:),spinZ(:,:,:)&
   &,spin2_T(:,:)
   logical :: runs, redimension
   real (q), allocatable :: lambdaI (:),lambdaR (:)
   integer, allocatable :: ind1(:),ind2(:),ind3(:)
   integer :: NN, Reason,p
!    real (q) :: spin_first

!
! Create many body configurations for spin excitation
! or ESR dynamics
!
! Brute force diagonalization of a model Hamiltonian
! that contains
!         spins, anisotropic exchange, local magnetic fields, Stephen Operators,
!         an exchange interaction between electron site and first spin
!         finite Hubbard U
!
!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!
! in this code, sites or spins are considered to be molecules, this explains the notation

! input file must exist otherwise stop and politely ask to fill in the info

inquire (file ='H_QD.input', exist = presence)
inquire (file ='Lehman_dimer.out', exist = presence2)

if (.not.presence) then
   write (*,*) '**************************************************'
   write (*,*) 'ERROR:'
   write (*,*) 'Please create an H_QD.input file to run this code.'
   write (*,*) 'This file will contain all the system information.'
   write (*,*) '              Stop.'
   write (*,*) '**************************************************'
elseif (presence2.and.(orb.ne.1))then
   write (*,*) '**************************************************'
   write (*,*) 'USING A TWO LEVEL ORBITAL HAMILTONIAN BY DIEGO'
   write (*,*) 'DIMENSION OF 14 WHEN PREVIOUSLY WE HAD 4'
   write (*,*) 'WE JUMP DIRECTLY TO READ THE LAMBDAS'
   write (*,*) '**************************************************'
   GOTO 111
   elseif (presence2.and.(orb.eq.1))then
   write (*,*) '**************************************************'
   write (*,*) 'USING A TWO LEVEL ORBITAL HAMILTONIAN BY DIEGO'
   write (*,*) 'BUT THE NUMBER OF ORBITAL IS SETTED TO 1,'
   write (*,*) 'CHANGE IT TO 2 OR LARGER. STOP'
   write (*,*) '**************************************************'
   stop
   elseif (orb.ne.1)then
   write (*,*) '**************************************************'
   write (*,*) 'You setted the orbital number to be higher than 1'
   write (*,*) 'but no dimmer file is in the folder. STOP'
   write (*,*) '**************************************************'
   stop
else
!  write (*,*) '**************************************************'
!  write (*,*) 'INFO:'
!  write (*,*) 'Reading H_QD.input'
!  write (*,*) ' '
   open (unit_input, file='H_QD.input', status='old')
   ! read input file
     read (unit_input,*)
     read (unit_input,*)
     read (unit_input,*)
     read (unit_input,*) Nm  ! Number of molecules (sites or spins), including the electron molecule (site)
                     ! for example, no spin, only one electron level is Nm=1
                     ! one single Ti is Nm=1
                     ! one Ti and one Fe is Nm=2
   
   ! we read info of each site:
   ! the electron level is always the first site
      allocate (Spin (Nm))! allocate spin of each site
      allocate (hx (Nm), hy (Nm), hz (Nm)) ! allocate one local magnetic field per site
      allocate (gx (Nm), gy (Nm), gz (Nm)) ! allocate gyromagnetic factors per site
      allocate (nx (Nm), ny (Nm), nz (Nm)) ! allocate one axis per site
      allocate (B20 (Nm), B22 (Nm), B40 (Nm), B44 (Nm)) ! Stephen Operators per site
         read (unit_input,*)
      do i_m = 1, Nm
         read (unit_input,*) Spin (i_m) ! read spin of each site including electron site

         read (unit_input,*) B20 (i_m), B22 (i_m), B40 (i_m), B44 (i_m) ! read local Stephen coefficients
                                                               ! Stephen coefficient units should be meV
         read (unit_input,*) nx (i_m), ny (i_m), nz (i_m) ! axis of the Stephen operators
                                                
         read (unit_input,*) hx (i_m), hy (i_m), hz (i_m) ! local magnetic field per site
                                                 ! in Teslas
         read (unit_input,*) gx (i_m), gy (i_m), gz (i_m) ! gyromagnetic vector
         read (unit_input,*)
      enddo
   ! we read the connection between sites
   ! in this code the connections are just anisotropic dipolar exchange interactions
     read (unit_input, *) Np ! Number of connected pairs including electronic site
      if (Np /=0) then
         allocate (mol1(Np), mol2(Np)) 
         allocate (Jexch(Np, 3)) ! anisotropic exchange: needs three components
                                 ! in GHz
         do i=1, Np
            read (unit_input,*) mol1 (i), mol2 (i) ! Indices of the exchange-connected molecules
            read (unit_input,*) Jexch(i,1), Jexch(i,2), Jexch(i,3) ! Three-component exchange in GHz
         enddo
      endif
            read (unit_input,*)
   ! finally, information on the electronic level:
     if (int(2*Spin(1)+1) /= 2) then
        write (*, *) ' '
        write (*,*) 'ERROR: '
        write (*,*) ' You are using a spin different from 1/2 for the transport electron!'
        write (*,*) ' Stop. '
        write (*,*) '  '
        stop
     endif
     read (unit_input, *) eps_QD  ! electronic level in meV
     read (unit_input, *) U_Hubbard ! Hubbard U (meV) leading to spin polarization SU(2) symmetry is not broken
     read (unit_input,*)
     read (unit_input,*) Name_output
     read (unit_input,*) Nplot ! Number of states to print in Spin_distribution.dat
   close (unit_input)
   write (*,*) 'INFO:'
   write (*,*) 'Reading Hamiltonian parametres finished '
   write (*,*) '**************************************************'
   write (*,*) ' '
   ! reading is finished
endif
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    Transform all quantities into atomic units                                                       !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tesla to a.u.
      hx = gx*hx*BohrMagneton*1000./Hartree; hy = gy*hy*BohrMagneton*1000./Hartree
      hz = gz*hz*BohrMagneton*1000./Hartree
!     print *, hx, hy, hz
! meV to a.u.
      B20 = B20/Hartree; B22=B22/Hartree; B40=B40/Hartree; B44=B44/Hartree
!GHz to a.u.
      if (Np/=0) then
      Jexch=Jexch*GHz/Hartree
      endif
! meV to a.u.
      eps_QD = eps_QD / Hartree; U_Hubbard=U_Hubbard / Hartree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
! Electronic Hamiltonian: first site, dimension N_in (1) = 4                                          !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate (N_in(Nm))
        N_in(1)=4  !int(2*Spin(1)+1)+1+1=4
        allocate (H_el(N_in(1),N_in(1)))
        H_el = zero
        H_el(1,1) = eps_QD; H_el(2,2) = eps_QD ! energy of the one electron state
        H_el(3,3)=0; H_el(4,4) = 2*(eps_QD)+U_Hubbard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    Establish basis set                                                                              !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the basis here is a tensorial product of the first state times the second times...
! each state is given by a spin matrix
! the spin matrix is called Ss (:,:,:,:)
!                  the first entry Ss (i,:,:,:) refers to the site
!                  the second entry Ss (:,i,:,:) refers to the x, y, z component of the spin
!                  the third and fourth entry coincide with the Hamiltonian dimensions Hamiltonian (:,:) 

! Total dimension (Hamiltonian entries) :: N
! each spin is of dimension N_in
! up to a spin i_m the total dimension is the product of previous N_in, stored in N_block
      allocate (N_block(Nm))
          N_block (1) = N_in (1)
! Test to remove after debugging
!         print *, 'N_in',1, N_in(1)
!         print *, 'N_block',1, N_block(1)
       do i_m = 2, Nm
          N_in (i_m)=int(2*Spin(i_m)+1)
          N_block (i_m) = N_block (i_m-1) * N_in (i_m)
!         print *, 'N_in',i_m, N_in(i_m)
!         print *, 'N_block',i_m, N_block(i_m)
       enddo
          N = N_block (Nm) ! full dimension
! MEMORY
      allocate (H (N,N)) ! contains the Hamiltonian before diag and the Eigenstates after
      allocate (W (N)) ! Eigenvalues
      allocate (Identity (N,N))
      allocate (Ss (Nm, 4, N, N)) ! The tensorial spin matrix -basis-
      allocate (Sn (3, N, N)) ! rotated spin matrix
      allocate (SProdS (3, N,N)) ! The following are operated spins 
      allocate (SProdS2 (3, N,N)) 
      allocate (Sp (N,N))
      allocate (Sm (N,N))
      allocate (Sp2 (N,N))
      allocate (Sm2 (N,N))
      allocate (Sp4 (N,N))
      allocate (Sm4 (N,N))

      

      
      
! initializations
     H = zero
     Ss = zero
     Sn = zero
     SProdS = zero
     SProdS2 = zero
     Identity = zero
        do i =1,N
           Identity (i,i) = ur
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the tensorial product of spin matrices (aka the basis set):                          !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call Kronecker_Product (N, Nm, N_in, N_block, Ss, Sx_u, Sy_u, Sz_u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the Hamiltonian                                                                      !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     STEPS:
!
!     ZERO: The electronic contribution to the Hamiltonian
!           we use the same algo as for the spin Ss

       do i2 = 1, N_in (1)
        do i3 = 1, N_in (1)

            do i4 = 1, N/N_block(1)

       H (i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  H_el(i2,i3);

            enddo
         enddo
        enddo



!     FIRST: local magnetic fields

      do i_ = 1, Nm ! Loop on molecules or sites
   
      H (:,:) = H (:,:)+Ss (i_, 1, :,:)*hx(i_)+Ss (i_, 2, :,:)*hy(i_)+Ss (i_, 3, :,:)*hz(i_)

!     SECOND: local Stephen operators
!             We need to compute spin matrices up to 4th power.

      ! We define new directions following the preferential axis given in the input file
      ! in this way we can take the internal "z" axis to be always along the preferential axis
      ! First normalize vector input:
          p_mod=sqrt(nx(i_)**2+ny(i_)**2+nz(i_)**2)
          nx (i_) = nx (i_)/p_mod; ny (i_) = ny (i_)/p_mod; nz (i_) = nz (i_)/p_mod
      ! we make vectors px,py,pz,pxx,... perpendicular to nx,ny,nz
       if (nz (i_) /= 0) then
          if (ny (i_) /= 0 ) then
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=nx (i_); pyy=ny (i_);pzz=-(nx(i_)**2+ny (i_)**2)/nz(i_);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else if (nx (i_) /=0) then
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=nx (i_); pyy=ny (i_);pzz=-(nx(i_)**2+ny (i_)**2)/nz(i_);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else
          pxx = 1.0; pyy = 0.0; pzz=0.0; px = 0.0; py = 1.0; pz=0.0
          endif
       else 
          px=-ny(i_); py=nx(i_); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=0; pyy=0;pzz=1._q
       endif


      ! rotated spins, inside the site loop on i_

       Sn (3,:,:)=Ss(i_,1,:,:)*nx (i_)+Ss(i_,2,:,:)*ny (i_)+Ss(i_,3,:,:)*nz (i_)
       Sn (2,:,:)=Ss(i_,1,:,:)*px+Ss(i_,2,:,:)*py+Ss(i_,3,:,:)*pz
       Sn (1,:,:)=Ss(i_,1,:,:)*pxx+Ss(i_,2,:,:)*pyy+Ss(i_,3,:,:)*pzz

      ! We perform matrix multiplication and use the definitions of Stephens Operators
      ! choosing a few ones up to 4th order -See our paper in J. Phys. Chem. A 124, 2318  (2020)
      ! We add over all molecules the contribution of each local anisotropy, this seems
      ! to work, at least for not extremely coupled spins


        call  MatSquareSpin (N, Sn, SprodS)
        call  MatSquareSpin (N, SprodS, SprodS2)

        H (:,:) = H (:,:)+ B20(i_)*SprodS(3,:,:)  !Longitudinal anisotropy !3 being include in B20 -> 3.*B20(i_)*SprodS(3,:,:) 
        H (:,:) = H (:,:)+ B22(i_)*(SprodS(1,:,:)-SprodS(2,:,:)) ! Transversal anisotropy
!         H (:,:) = H (:,:)+ B22(i_)*(SprodS(1,:,:)+SprodS(2,:,:)) ! Transversal anisotropy
        H (:,:) = H (:,:)+ B40(i_)*35.*SprodS2(3,:,:) ! Fourth order: B40 ain't pretty
        H (:,:) = H (:,:)- B40(i_)*30.*(Spin(i_)*(Spin(i_)+1._q)-2*Spin(i_))*SprodS(3,:,:)
        H (:,:) = H (:,:)+ B40(i_)*(3*(Spin(i_)*(Spin(i_)+1))**2-6*(Spin(i_)*(Spin(i_)+1)))*Identity (:,:)

        ! change to circular spins: S+ and S-
        Sp (:,:) = Sn (1, :,:)+ui*Sn (2, :,:)
        Sm (:,:) = Sn (1, :,:)-ui*Sn (2, :,:)

        call MatSpSm (N, Sp, Sm, Sp2, Sm2)
        call MatSpSm (N, Sp2, Sm2, Sp4, Sm4)

        H = H + B44(i_)*0.5*(Sp4+Sm4) ! Fourth order: B44, sort of fourth-order longitudinal anisotropy
!         H = H + B22(i_)*0.5*(Sp2+Sm2) ! Transversal anisotropy
!         H = H - B44(i_)*0.5*ui*(Sp4-Sm4) ! Fourth order: B44, sort of fourth-order longitudinal anisotropy
   
      enddo
! Test to remove after debugging
!     do i = 1, N
!         write (*, *) (j,H(i,j), j=1, N)
!     enddo


!     THIRD: anisotropic intermolecular exchange interactions
!     we keep the original axis, not the one of the anisotropy

      do i = 1, Np ! loop on pairs

!     if (i ==1) then
!     print *, 'computing first pair!!!'
!     else if (i == 2) then
!     print *, 'computing second pair!!!'
!     else
!     print *, 'computing', i,'th pair!!!'
!     endif

      call MatProdSpin (N, i, mol1, mol2, Ss, SProdS) ! mol1 and mol2 contain the indices of the paired molecules

! Test: remove after debugging
!     print *, 'First Spin'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(Ss(1,3,j1,j)), j=1, 8)
!     enddo
!     print *, 'Second Spin'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(Ss(2,3,j1,j)), j=1, 8)
!     enddo
!     print *, 'Product:'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(SProdS(3,j1,j)), j=1, 8)
!     enddo

      H (:,:) = H (:,:)+Jexch (i, 1)*SProdS (1,:,:)+Jexch (i, 2)*SProdS (2,:,:)+Jexch (i, 3)*SProdS (3,:,:)

      enddo

! Test: remove after debugging
!      print *, 'The Hamiltonian is:'
!    do i = 1,N
!      write (*,'(16g14.4)')  (j,Real(H(i,j)), aimag(H(i,j)), j=1,N)
!    enddo

     if (runs) then
! DO NOT diagonalize if Name_output exists, and use instead the values
! of a previous run
    inquire (file =Name_output, exist = presence)
    else
        presence = .false. !if runs =.false. do not read
    endif
    if (presence) then
      open (unit_spinoutput, file=Name_output, form='unformatted')
      read (unit_spinoutput) H, W
      close (unit_spinoutput)

           print *, ' '
           print *, 'WARNING!!!!!'
           print *, ' '
           print *, ' Reading the Hamiltonian from a previous run!!!'
           print *, ' '
           print *, 'WARNING!!!!!'
           print *, ' '


       ! check the runs are compatible
         if (N /= size (H,1)) then
           print *, ' '
           print *, Name_output, 'IS NOT COMPATIBLE with H_QD.input'
           print *, 'remove', Name_output, 'from the running folder.'
           print *, 'STOP.'
           print *, ' '
           stop
         endif
    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        DIAGONALIZE                                                                                  !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    print *, ' '
!    print *, 'Beginning of diagonalization'
!    print *, ' '
!          print *, 'size H:', size(H,1), N
!     do i =1, N
!     write(126,*) i,(dble(H(i,j)*hartree),aimag(H(i,j)*hartree), j= 1, N)
!     enddo
      call DIAGONALIZE (N,H,W)
   print *, 'Hamiltonian diagonalized'
!    print*,''
!    print*, 'Key energies of the Hamiltonian'
!    print*, 'Zero electron energy:',-Ef*Hartree
!    print*, 'Single level:',(eps_QD-Ef)*Hartree
!    print*, 'Doubly occupied energy:',(2*(eps_QD)+U_Hubbard-Ef)*Hartree
   
   print*,''
   print*, 'Key energies for transport:'
   print*, 'Ionization energy:',eps_QD*Hartree
   print*, 'Coulomb repulsion:',U_Hubbard*Hartree
   print*, 'On site energy plus Ionization:',(U_Hubbard+eps_QD)*Hartree
   print*,''
   print *, 'These are the eigenvalues (GHz)::'
   write (*, '(8g14.4)') ((W(i))*Hartree/GHz, i=1,N)
   print *, 'Eigenvalues (meV)::'
   write (*, '(8g14.4)') ((W(i))*Hartree, i=1,N)
   print *, ' '

! Save it, unformatted so it does not need to be recalculated if it 
! exists on the running folder

      open (unit_spinoutput, file=Name_output, form='unformatted')
      write (unit_spinoutput) H, W
      close (unit_spinoutput)

      
    Ndim=N  
    p=0
    if ((redimension).and.(Nd.gt.N/2).and.(Nd<N)) then
        
!     write (*,*) ' '
    write (*,*) 'Redimensioning from',N,'to New dimension',Nd,'from the input file.'
     Ndim = Nd
!     write (*,*) 'Extra zeros will appear in the eigenvectors due to the old basis being cut'
!
    p=1 ! it will indicate that a redimensionalition was done so we do not do another one
    write (*,*) ' '
    else
    write (*,*) ' '
    write(*,*) 'Input dimension',Nd,' below min dimension',N/2,'above or equal max dimension',Ndim
    write(*,*) 'or redimension is set to be .FALSE.. Redimensioning based on energies differences'
!     redimension=.false.
    write (*,*) ' '
   endif
      
      ! Calculation of output
    allocate (Eigenvalues (N))
    allocate (Delta (N,N))
    allocate (lambda (N,N,2))

        Eigenvalues = W
    do i = 1, N
    do j = 1, N
        Delta (i,j) = W(i)-W(j)
    enddo
    enddo

    
    ! redimension of the problem if the input redimension does not apply:
    ! Nd is a number larger than N or more states are disconnected
    
    if (abs(bias_L).lt.abs(bias_R))then
        bias=abs(bias_R)
    else
        bias=abs(bias_L)
    endif

    if ((redimension).and.(p.eq.0)) then
        call redimensioning (Ndim, Delta, bias, Cutoff)
    elseif (p.eq.0) then
        Ndim=N
        write (*,*) 'Full dimension is taken',N
        write (*,*) ' '
    endif
    
      
      
   print *, 'States and energies written in: Eigenvalues and eigenvectors .dat'
!    print *, ' '
! We go ahead and write the ouput 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        Write output in a useful way                                                                 !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    print *, 'Eigenvalues in meV are written in Hamiltonian_Eigen.dat'
!    print *, ' '
       open (unit_spinoutput, file='Eigenstates.dat')
       open (unit_spinoutput+1, file='Eigenvalues.dat')
       open (unit_spinoutput+2, file='Matrix_transition.dat')

       do i =1, Ndim
          write (unit_spinoutput,*) i,(real(H(j,i)),aimag(H(j,i)), j= 1, N) ! IF REDIMENSION THEN EXTRA ZEROS WILL APPEAR FROM THE OLD BASIS CUT
          
write (unit_spinoutput+1,*) i,hartree*W(i),hartree*W(i)/GHz,hartree*(W(i)-W(1)),hartree*(W(i)-W(1))/GHz
       enddo
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,(hartree*(W(i)-W(j)),j=1,Ndim)
        end do

        write(unit_spinoutput+2,*) ''
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)        
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,((W(i)-W(j))*Hartree/GHz,j=1,Ndim)
        end do
       close (unit_spinoutput)
       close (unit_spinoutput+1)
       close (unit_spinoutput+2)

! Plotting the spin distribution for the first Nplot states 


    if (Nplot > N) then
       Nplot = N
    endif

    allocate (spinX(Nplot,Nplot,Nm),spinY(Nplot,Nplot,Nm),spinZ(Nplot,Nplot,Nm))
    allocate (spin2_T(Nplot,Nplot))
      spinX = zero
      spinY = zero
      spinZ = zero
      open (unit_spinoutput, file='Spin_distribution.dat')
    
    do j = 1, Nm
    do i=1, Nplot
    do ii=1, Nplot
       write (unit_spinoutput,*) 'State=',i
       write (unit_spinoutput,*) '#  Site, Sx, Sy, Sz'

      do j1=1, N
      do j2=1, N
      spinX (i,ii,j)=spinX (i,ii,j)+conjg(H(j1,i))*Ss (j, 1, j1, j2)*H(J2, ii)
      spinY (i,ii,j)=spinY (i,ii,j)+conjg(H(j1,i))*Ss (j, 2, j1, j2)*H(J2, ii)
      spinZ (i,ii,j)=spinZ (i,ii,j)+conjg(H(j1,i))*Ss (j, 3, j1, j2)*H(J2, ii)
      enddo
      enddo

    enddo        
    enddo
    enddo
    
    do i=1,Ndim
        write (unit_spinoutput,*) i, (real(spinX(i,i,j)), real(spinY(i,i,j)), real(spinZ(i,i,j)),j=1,Nm)
    enddo
    
    
    
   print *, 'Spins written in file:  Spin_distribution.dat'
   print *, ' '

        spin2_T = zero  ! spin square

      close (unit_spinoutput)
      
    open (unit_spinoutput, file='Full_energies_and_spin.dat') ! costy writing if total spin is large

      write (unit_spinoutput,*) ' State ', ' Excitation Energy (GHz) ', ' (meV) ', ' Spin^2 ', ' Spin '

    do i = 1, Nplot
    do ii= 1, Nplot
      do j = 1, Nm
       do j4 = 1, Nm
        do j1=1, N
         do j2=1, N
          do j3=1,N
               spin2_T(i,ii)=spin2_T(i,ii)+conjg(H(j1,i))*(Ss (j, 1, j1, j3)*Ss (j4, 1, j3, j2)+  &
      &        Ss (j, 2, j1, j3)*Ss (j4, 2, j3, j2)+ Ss (j, 3, j1, j3)*Ss (j4, 3, j3, j2) &
      &        )*H(J2, ii)
          enddo
         enddo
        enddo
       enddo
      enddo

    enddo
    write (unit_spinoutput,*) i, (W(i))*Hartree/GHz, (W(i))*Hartree, &
        real(spin2_T(i,i)), 0.5*(sqrt(1+4*real(spin2_T(i,i)))-1)
    enddo
    
!     do i=1, Nplot
!         spin_first=real(spinX(i,1))**2+ real(spinY(i,1))**2 +real(spinZ(i,1))**2
!     enddo
    
    
      close (unit_spinoutput)
   print *, 'The Hamiltonian calculation is DONE!'
   write (*,*) '**************************************************'
  endif !end of diverting if Name_output exists


! Calculation of lambda

    lambda = zero

    do i = 1, N
    do j = 1, N
    do i_sigma= 1, 2 !1 spin down as always in this code

! contribution from |0Xsigma|

      i2 = 3 ! 3 is |0>
      i3 = i_sigma ! 1 is down and 2 is up as corresponds to the basis set

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo

! contribution from |\bar{sigm}aX4| where 4 is the doubly occupied state (singlet)

      i2 = (-1)**(i_sigma+1)+i_sigma! if i_sigma=1 then this yields 2
                            ! if i_sigma=2 then this yields 1
      i3 = 4 ! 4 is the doubly occupied state

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo


    enddo
! test lambda
!      write (123,*) i,j, 1, lambda (i,j,1)
!      write (123,*) i,j, 2, lambda (i,j,2)

    enddo
    enddo

    return 
    
      111 write(*,*)''

      ! DIEGO PROGRAM BELOW UNTIL THE END OF THE SUBROUTINE

! Calculation of lambda


    allocate (lambdaR(100000),ind3(100000)) ! large dimension to read any file
    allocate (lambdaI(100000),ind1(100000),ind2(100000))

    lambdaR = zero
    lambdaI = zero
    ind1=zero
    ind2=zero
    ind3=zero

    open(unit=71, file='Lehman_dimer.out', status='old')

    N=0

    outer: DO i=1,100000
        READ(71,*,IOSTAT=Reason) ind1(i),ind2(i),ind3(i),lambdaR(i),lambdaI(i)
        IF (Reason > 0)  THEN
            write(*,*) 'Something wrong with the file, stop'
            stop
        ELSEIF (Reason < 0) THEN
            write(*,*) 'Lambda file read'
            exit outer
        ELSE
            N=N+1 ! increase dimension
        END IF
    ENDDO outer

    print*,N,N/(2*orb),sqrt(real(N/(2*orb))) ! for any number of orbitals, change number 4

    NN=N
    N=N/(2*orb) ! block of 4: site 1, spin up; site 1 spin down; site 2, spin up ...
    N=sqrt(real(N)) ! it should be a interger, N dimension column Hamiltonian
    N=int(N)
    allocate (lambda(N,N,(2*orb))) ! We add the correct dimension of the file
    lambda = zero

    close(71)

    DO j=1,NN
        if ((lambdaR(j)<1E-15).and.(lambdaI(j)<1E-15))then
    lambda(ind1(j),ind2(j),ind3(j))=zero
        else
    lambda(ind1(j),ind2(j),ind3(j))=ur*lambdaR(j)+ui*lambdaI(j)
        endif
    ENDDO


    open(unit=71, file='Energies_dimer.out', status='old')


    allocate (Eigenvalues (N))
    allocate (Delta (N,N))
    allocate (W(N))

    do j=1,N
        READ(71,*) Eigenvalues(j) ! reading energies from Energies_dimer.out
    enddo
    close(71)
    Eigenvalues=Eigenvalues/Hartree
    W=Eigenvalues
    do i = 1, N
        do j = 1, N
            Delta (i,j) = W(i)-W(j)
!             write(*,*)Delta(i,j),i,j
        enddo
    enddo


    write (*,*) ' '
    write (*,*) 'Ajusting the dimension to the dimer one'
     Ndim = N
    write (*,*) 'Dimension',Ndim
    write (*,*) ' '

    ! OUTPUT exactly as we did not have orbital structure

        open (unit_spinoutput, file='Eigenstates.dat')
       open (unit_spinoutput+1, file='Eigenvalues.dat')
       open (unit_spinoutput+2, file='Matrix_transition.dat')

       do i =1, Ndim

write (unit_spinoutput+1,*) i,hartree*W(i),hartree*W(i)/GHz,hartree*(W(i)-W(1)),hartree*(W(i)-W(1))/GHz
       enddo
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form in meV)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,(hartree*(W(i)-W(j)),j=1,Ndim)
        end do

        write(unit_spinoutput+2,*) ''
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form in GHz)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,((W(i)-W(j))*Hartree/GHz,j=1,Ndim)
        end do
       close (unit_spinoutput)
       close (unit_spinoutput+1)
       close (unit_spinoutput+2)

! Plotting the spin distribution for the first Nplot states


    if (Nplot > N) then
       Nplot = N
    endif

        allocate (spinX(Nplot,Nplot,Nm),spinY(Nplot,Nplot,Nm),spinZ(Nplot,Nplot,Nm))
    allocate (spin2(Nplot,Nplot,Nm),spin2_T(Nplot,Nplot))
      spinX = zero
      spinY = zero
      spinZ = zero
      spin2 = zero
      open (unit_spinoutput, file='Spin_distribution.dat')

    do j = 1, Nm
    do i=1, Nplot
    do ii=1, Nplot
!       write (unit_spinoutput,*) 'State=',i
!       write (unit_spinoutput,*) '#  Site, Sx, Sy, Sz'

      do j1=1, N
      do j2=1, N
      spinX (i,ii,j)=spinX (i,ii,j)+conjg(H(j1,i))*Ss (j, 1, j1, j2)*H(J2, ii)
      spinY (i,ii,j)=spinY (i,ii,j)+conjg(H(j1,i))*Ss (j, 2, j1, j2)*H(J2, ii)
      spinZ (i,ii,j)=spinZ (i,ii,j)+conjg(H(j1,i))*Ss (j, 3, j1, j2)*H(J2, ii)
      spin2(i,ii,j)=spin2(i,ii,j)+real(spinX(i,ii,j))**2 +real(spinY(i,ii,j))**2+ real(spinZ(i,ii,j))**2
      enddo
      enddo

    enddo
        write (unit_spinoutput,*) i, j, real(spinX(i,i,j)), real(spinY(i,i,j)), real(spinZ(i,i,j)),&
        real(spin2(i,i,j)),0.5*(sqrt(1+4*real(spin2(i,i,j)))-1)

    enddo
    enddo


   print *, 'Spins written in file:  Spin_distribution.dat'
   print *, ' '

        spin2_T = zero  ! spin square

      close (unit_spinoutput)

    open (unit_spinoutput, file='Full_energies_and_spin.dat')

      write (unit_spinoutput,*) ' State ', ' Excitation Energy (GHz) ', ' (meV) ', ' Spin^2 ', ' Spin '

    do i = 1, Nplot
    do ii= 1, Nplot
      do j = 1, Nm
       do j4 = 1, Nm
        do j1=1, N
         do j2=1, N
          do j3=1,N
               spin2_T(i,ii)=spin2_T(i,ii)+conjg(H(j1,i))*(Ss (j, 1, j1, j3)*Ss (j4, 1, j3, j2)+  &
      &        Ss (j, 2, j1, j3)*Ss (j4, 2, j3, j2)+ Ss (j, 3, j1, j3)*Ss (j4, 3, j3, j2) &
      &        )*H(J2, ii)
          enddo
         enddo
        enddo
       enddo
      enddo

    enddo
    write (unit_spinoutput,*) i, (W(i))*Hartree/GHz, (W(i))*Hartree, &
        real(spin2_T(i,i)), 0.5*(sqrt(1+4*real(spin2_T(i,i)))-1)
    enddo
      close (unit_spinoutput)
   print *, 'The Hamiltonian calculation is DONE!'
   write (*,*) '**************************************************'
    
    
    
    ! END DIEGO EXTENSION
    
    
    
    
    
   end subroutine Hamiltonian
!
! Cutting off dimensions for sizeable calculations
! 
     subroutine redimensioning (Ndim, Delta, bias, Cutoff)
     implicit none
     integer, intent(inout) :: Ndim
     integer :: i, l
     real (q), intent(in) :: bias, Cutoff
     real (q), dimension (:,:), intent(in) :: Delta

       l = 0
        
     do i = 1, Ndim

       if (Delta (i,1) <= bias+10*Cutoff) then
           l = l+1
       endif

     enddo

    write (*,*) 'Redimensioning from',Ndim,'to'
    Ndim = l
    write(*,*) 'New dimension',l,'based on energy differences.'
    write (*,*) ' '


     return
     end subroutine redimensioning


end module H_QD
