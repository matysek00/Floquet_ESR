program Floquet_ESR

! Long-time solution for CW ESR
! based on the theory by Galvez et al
! Phys. Rev. B Physical Review B 104 (24), 245435 (2021)

!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!  Version 1.0 Including Lamb shift with cutoff
!  January 2023

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the main file and calls in:
Use declarations !all variables, constants etc are defined here
Use OpenFiles !all units of files in open statements are here
!Use feed_back_on
Use io !(input/output) the reading and writing routines are here
Use timer
Use H_QD !(Hamiltonian Quantum Dot) contains the Hamiltonian and solution
Use QME_F !(Quantum Master Equation) contains rates and The Matrix
Use Transport_F ! computation of electron current
Use SpinTools
Use Matrix_Coeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          input run's values                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call clock ('STEP 1:: Reading INPUT and Hamiltonian for the QME', 1)

     call reading_input ( gamma_R_0, A_R, gamma_L_0, A_L, phi, orb,&
       NF,NCF,GammaC,Cutoff,WW,gau,redimension,Nd, frequency,&
       bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, B_R, B_L, p_max, write_populations, write_coherences, spinflo,Ef,FermiP,seHa,&
       feedbackon,Iset,tol,VDC,ratio)

!      call clock ('Finished reading INPUT for the QME ', 2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the spin+orbital system                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call clock ('STEP 2:: Diagonalizing Hamiltonian ', 1)

    call Hamiltonian (N,lambda,Delta,Eigenvalues,H,Ss,spinX,spinY,spinZ,&
    &spin2_T,runs,Nm,hx,hy,hz,Ndim,&
    &redimension,orb,bias_R, bias_L,Cutoff)

     call clock ('STEP 2:: Finished diagonalizing Hamiltonian and reading', 2)

    Nmatrix=Ndim*Ndim*NF ! total dimension matrix A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      ALLOCATES                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate (G (Ndim, Ndim, Ndim, Ndim, NF, 2))
     allocate (GC (Ndim, Ndim, Ndim, Ndim, NF, 2))
     allocate (GA (Ndim, Ndim, Ndim, Ndim, NF), GCA (Ndim, Ndim, Ndim, Ndim, NF))
     allocate (A (Nmatrix, Nmatrix),A_old(Nmatrix,Nmatrix))
     allocate (B (Nmatrix))
     allocate (X (Nmatrix),XX(Nmatrix))
     allocate (rho(Ndim,Ndim))
     allocate (WORK (Nmatrix))
     allocate (RWORK (Nmatrix))
     allocate (SWORK (Nmatrix*(Nmatrix+1)))
     allocate (IPIV (Nmatrix))
     allocate (curr(NF))

! loop on driving freq
!uencies. Parallelized with coarrays

! we open the outfile here too
        k=0
        open (unit_curr, file='Current_0.dat')
        open (unit_error, file='strange_populations.dat')
        open (unit_output+k, file='POPULATIONS_0.dat')
        k=k+1
        open (unit_output+k, file='POPULATIONS_1.dat')
        k=k+1
        open (unit_output+k, file='POPULATIONS_n1.dat')
        k=k+1
        open (unit_output+k, file='COHERENCES_0.dat')
        k=k+1
        open (unit_output+k, file='COHERENCES_1.dat')
        k=k+1
        open (unit_output+k, file='COHERENCES_n1.dat')
        k=k+1
        open (unit_output+k, file='SpinFloquet_0.dat')
        k=k+1
        open (unit_output+k, file='SpinFloquet_1.dat')
        k=k+1
        open (unit_output+k, file='SpinFloquet_n1.dat')


! add more floquet numbers (\pm 2) if the driving is large enough to make them important

A=zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    Solve QME in Floquet basis                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call clock ('STEP 3:: Calculating rates off resonance freq and feedback loop', 1)
! the off resonance frequency should be too large

! First evaluate rates in the H_QD basis and Floquet indices
!      call clock ('STEP 3:: Computing rates ', 1)
  
  !     right electrode
  call rates (Ndim, orb, frequency, gamma_R_0, lambda, Spin_polarization_R,&
            NCF, p_max, A_R, phi, B_R, GammaC, bias_R, Delta, Cutoff, Temperature,& 
            seHa, WW, gau, N_int, GA, GCA)
  G(:,:,:,:,:,1)  = GA(:,:,:,:,:)
  GC(:,:,:,:,:,1) = GCA(:,:,:,:,:)
  
  !     left electrode
  call rates (Ndim, orb, frequency, gamma_L_0, lambda, Spin_polarization_L,&
    NCF, p_max, A_L, phi, B_L, GammaC, bias_L, Delta, Cutoff, Temperature,& 
    seHa, WW, gau, N_int, GA, GCA)
  G(:,:,:,:,:,2)  = GA(:,:,:,:,:)
  GC(:,:,:,:,:,2) = GCA(:,:,:,:,:)
  
  call coeff_matrix (Ndim, frequency, NF, NCF, G, Nmatrix, A, B)
  
  open (unit_rates, file='lala.dat')
  do i=1,Nmatrix
  do j=1,Nmatrix
  if (A(i,j).ne.zero) then
          write (unit_rates,*) i,j, dble(A(i,j)), imag(A(i,j))
  endif     
  enddo
  enddo
  close(unit_rates)
!     call clock ('STEP 3:: Finished computing rates ', 2)
    ! Solve the linear equations A*x=B
! using LAPACK zcgesv for square matrices

!      call clock ('STEP 4:: Solving linear system of QME ', 1)
     NRHS = 1 !one column for B
!

!       Solve the equations A*X = B. X retuned in B 
!       A  is n x m matrix
!       B, X are  l x k matrix
!       zgsev (n, l, A, m, IPIV, B, l, INFO)
    call zgesv(Nmatrix,1,A,Nmatrix,IPIV,B,Nmatrix,INFO)
    ! better solver and similar performance
    X = B
    XX=zero
    do i=1,NF
        XX(1+ndim*ndim*(NF-i):ndim*ndim*(NF+1-i))=X(1+ndim*ndim*(i-1):ndim*ndim*i)
    enddo

           if (INFO > 0) then
               print *, 'The linear problem is singular and cannot be solved, check input: dimension and bias. INFO=', INFO
!                i, U(i,i) computed in DOUBLE PRECISION is
!                 exactly zero.  The factorization has been completed,
!                 but the factor U is exactly singular, so the solution
!                 could not be computed.
           else if (INFO < 0) then
               print *, 'if INFO = -i, the i-th argument had an illegal value. INFO=', INFO
           endif

! Compute the DC electron current
     call Current (Ndim, NF, NCF, X, XX, GC, curr, Electrode)

if (feedbackon) then
  print *, 'Feedback not implemented. Sorry :('
endif
!if (feedbackon) then
!  open (unit_feedbackon, file='feedbackon_data.dat')
!  print*,'Feedback on. We will increase/decrease polarized coupling (moving tip...)'
!  print*,'Initial values for the current and couplings',curr(NCF+1)*pA, gamma_L_0*hartree,gamma_R_0*hartree
!  print*,''
!  call feedback(Spin_polarization_L,Spin_polarization_R,gamma_L_0,gamma_R_0,&
!   &curr,Iset,tol,bias_L,bias_R,feedbackon,GC,G,A_fast_left,A_fast_right,X,XX,ratio)
!   close(unit_feedbackon)
!    call clock ('Feedback on loop ended. New rates computed', 2)
!else
!   call clock ('Tip remains static. Feedback off. Rates computed.', 2)
!endif

  ! Here we will write all rates of zero floquet to be analyzed later on
     open (unit_rates, file='rates_floquet.dat')
     open (unit_rates+1, file='rates_floquet0.dat')

    do l=1,Ndim
    do j=1,Ndim
        do u=1,Ndim
        do v=1,Ndim
        if ((G (l,j,u,v,NCF+1,1)+G (l,j,u,v,NCF+1,2)).ne.zero)then
write(unit_rates,*)hartree*dble(G(l,j,u,v,NCF+1-2:NCF+1+2,2)),hartree*dimag(G(l,j,u,v,NCF+1-2:NCF+1+2,2)),&
    &hartree*dble(G (l,j,u,v,NCF+1-2:NCF+1+2,1)), hartree*dimag(G (l,j,u,v,NCF+1-2:NCF+1+2,1)), l,j,u,v, &
    &frequency/(2*pi_d*time_unit)
    ! WE ONLY WRITE 5 FLOQUET NUMBERS: -2,-1,0,1,2.
    write(unit_rates+1,*)hartree*dble(G(l,j,u,v,NCF+1,2)),hartree*dimag(G(l,j,u,v,NCF+1,2)),&
    &hartree*dble(G (l,j,u,v,NCF+1,1)), hartree*dimag(G (l,j,u,v,NCF+1,1)), l,j,u,v, &
    &frequency/(2*pi_d*time_unit)
    ! simplified version of writing the rates, only Floquet 0
        endif
        enddo
        enddo
    enddo
    enddo

     close(unit_rates)
     close(unit_rates+1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    output the results                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Current

!      call clock ('STEP 6:: Writing output ', 1)
write (unit_curr, *) frequency/(2*pi_d*time_unit), curr(NCF+1)*pA,curr(NCF+2)*pA,curr(NCF+3)*pA ! 0,1,2

           do l=1,Ndim
            i = 1 + ndim*ndim*(NCF) + (l-1)*(Ndim+1)
              rho (l,l) = X(i)
              sum_rule  = 0
              do i = 1, Ndim
                sum_rule  = sum_rule  + dble (rho(i,i))
              enddo

              if((dble(rho(l,l))<-1.E-8).or.(dble(rho(l,l))>1.000001_q).or.(dabs(dimag(rho(l,l)))>1E-6)&
              &.or.(sum_rule>1.000001_q))then
                write(unit_error,*) rho(l,l),l,sum_rule,Ndim,bias_R*hartree,bias_L*hartree,frequency/(2*pi_d*time_unit)
                write(*,*) 'WEIRD RESULTS! -> CHECK strange_populations.dat but we continue...'
!                 stop
              end if
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit), (dble(rho (l,l)), l= 1, Ndim)
            k=k+1


! Floquet n = 1
           do l=1,Ndim
            i = 1 + ndim*ndim*(NCF+1) + (l-1)*(Ndim+1)
              rho (l,l) = X(i)
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit), (dble(rho (l,l)), l= 1, Ndim)
            k=k+1
! Floquet n = -1
           do l=1,Ndim
            i = 1 + ndim*ndim*(NCF-1) + (l-1)*(Ndim+1)
              rho (l,l) = X(i)
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit), (dble(rho (l,l)), l= 1, Ndim)
            k=k+1

      if (write_coherences) then
        i1=1
        ! Floquet n = 0
          do i=1+ndim*ndim*(NCF),ndim*ndim*(NCF+1)
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit),dble(rho (:,:)),dimag(rho (:,:))
              ! the order is alter by the way fortran print this array so we have for Ndim=4
              ! t Re 11 Re 21 Re 31 Re 41 Re 12 Re 22 Re 32  Re 42 ...
              ! Im 11 Im 21 Im 31 Im 41 Im 12 Im 22 Im 32 Im 42 ...

            k=k+1
        i1=1
        ! Floquet n = 1
          do i=1+ndim*ndim*(NCF+1),ndim*ndim*(NCF+2)
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit),dble(rho (:,:)),dimag(rho (:,:))
            k=k+1
        i1=1
        ! Floquet n = -1
           do i=1+ndim*ndim*(NCF-1),ndim*ndim*(NCF)
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo
              write (unit_output+k,*) frequency/(2*pi_d*time_unit),dble(rho (:,:)),dimag(rho (:,:))
            k=k+1
      endif


      if (spinflo) then
! if (i_omega==1)then
!       write(*,*)omega, spinX
!       endif
        i1=1
          do i=1+ndim*ndim*(NCF),ndim*ndim*(NCF+1) ! Floquet n = 0
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo

    call SpinFloquet (Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, rho, hx, hy, hz,&
    &Sx,Sy,Sz,Sh,spin2_ave)

            write (unit_output+k,*) frequency/(2*pi_d*time_unit),&
            &(real(Sx(i)), real(Sy(i)), real(Sz(i)), real(Sh(i)),  i=1, Nm),&
            & real(spin2_ave),real(sqrt(1+4*(spin2_ave))-1)*0.5
            k=k+1
        i1=1

          do i=1+ndim*ndim*(NCF+1),ndim*ndim*(NCF+2) ! Floquet n = 1
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo
    call SpinFloquet (Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, rho, hx, hy, hz,&
    &Sx,Sy,Sz,Sh,spin2_ave)

            write (unit_output+k,*) frequency/(2*pi_d*time_unit),&
            &(real(Sx(i)), real(Sy(i)), real(Sz(i)), real(Sh(i)),  i=1, Nm),&
            & real(spin2_ave),real(sqrt(1+4*(spin2_ave))-1)*0.5
            k=k+1

        i1=1
           do i=1+ndim*ndim*(NCF-1),ndim*ndim*(NCF) ! Floquet n = -1
            l=1+(i1-1)/ndim
            j=i1-(l-1)*ndim
            i1=i1+1
              rho (l,j) = X(i)
           enddo
    call SpinFloquet (Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, rho, hx, hy, hz,&
    &Sx,Sy,Sz,Sh,spin2_ave)

            write (unit_output+k,*) frequency/(2*pi_d*time_unit),&
            &(real(Sx(i)), real(Sy(i)), real(Sz(i)), real(Sh(i)),  i=1, Nm),&
            & real(spin2_ave),real(sqrt(1+4*(spin2_ave))-1)*0.5
            k=k+1
      endif

close (unit_curr)
close (unit_error)
do i=unit_output,unit_output+k
    close (i)
enddo
    write(*,*) ''
    write(*,*) '************************************************************************************************'
    write(*,*) 'Calculation done and output written!!!'
    write(*,*) ''
    write(*,*) 'Every output is divided by Floquet numbers 0,-1,1 (n indicates the minus sign)'
    write(*,*) ''

        if (write_populations) then
            write(*,*) 'Populations are written in POPULATIONS_'
        endif
    write(*,*) ''
        if (write_coherences) then
            write(*,*) 'Coherences are written in COHERENCES_'
        endif
    write(*,*) ''
        if (spinflo) then
            write(*,*) 'Spin Sx,Sy,Sz,Sh per site are written in SpinFloquet_'
        endif

    write(*,*) ''
     write(*,*) 'No zero real and imag left and right rates for Floquet numbers -2,-1,0,1,2 written in rates_floquet.dat'
     write(*,*) 'while rates_floquet0.dat contains only the Floquet zero'
    write(*,*) ''
    write(*,*) 'DC current wrote in Current_0.dat'

     call clock ('Final STEP:: Everything done', 2)

end program Floquet_ESR
