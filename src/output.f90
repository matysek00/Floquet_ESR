module output
use declarations
use OpenFiles !all units of files in open statements are here
use SpinTools
CONTAINS

    subroutine write_rate_out(G, NCF, Ndim, frequency)
    implicit none
    ! Declare all arguments and local variables
    complex(q), intent(in) :: G(:,:,:,:,:,:)
    real(q), intent(in) :: frequency
    integer, intent(in) :: NCF, Ndim

    integer :: l, j, u, v
    complex(qc) :: G_temp(NCF, 2)

    open (unit_rates, file='rates_floquet.dat')
    open (unit_rates+1, file='rates_floquet0.dat')
!   TODO: add a header to the file

    do l=1,Ndim
    do j=1,Ndim
    do u=1,Ndim
    do v=1,Ndim
        G_temp = G(l,j,u,v,:,:)*hartree ! convert to Hartree

        if ((G_temp(NCF+1,1)+G_temp(NCF+1,2)) == (0.0_q,0.0_q)) then
!           if the rates are zero, we do not write them
            cycle
        end if

        write (unit_rates,*) frequency/(2*pi_d*time_unit), l,j,u,v,& 
            dble(G_temp(:,2)), dimag(G_temp(:,2)),&
            dble(G_temp(:,1)), dimag(G_temp(:,1))

!       simplified version of writing the rates, only Floquet 0  
        write (unit_rates+1,*) frequency/(2*pi_d*time_unit), l,j,u,v,& 
            dble(G_temp(NCF+1,2)), dimag(G_temp(NCF+1,2)),&
            dble(G_temp(NCF+1,1)), dimag(G_temp(NCF+1,1))
    end do
    end do
    end do
    end do

    close(unit_rates)
    close(unit_rates+1)

    return
end subroutine write_rate_out


subroutine check_populations(Rho, NCF, Ndim, bias_R, bias_L, frequency)
    implicit none
    real(q), intent(in) :: bias_R, bias_L, frequency
    complex(qc), intent(in) :: Rho(:,:,:)
    integer, intent(in) :: NCF, Ndim
    integer :: l
    real(q) :: sum_rule

    open (unit_error, file='strange_populations.dat')
    do l=1,Ndim

        sum_rule  = sum_rule  + dble (Rho (l,l,NCF))
!       TODO: This is not working I need to fix it
        if((dble(Rho(l,l,NCF))<-1.E-8).or.(dble(Rho(l,l,NCF))>1.000001_q).or.(dabs(dimag(Rho(l,l,NCF)))>1E-6)&
                &.or.(sum_rule>1.000001_q))then

            write(unit_error,*) Rho(l,l,NCF),l,sum_rule,Ndim,bias_R*hartree,bias_L*hartree,frequency/(2*pi_d*time_unit)
            write(*,*) 'WEIRD RESULTS! -> CHECK strange_populations.dat but we continue...'
        end if
    enddo
    close (unit_error)
    
    return
end subroutine check_populations


subroutine write_cur_out(current, frequency)
    implicit none
    real(q), intent(in) :: current(:)
    real(q), intent(in) :: frequency
    
    open (unit_curr, file='Current_0.dat')
    write (unit_curr, *) 'Frequency (GHz) / Current (pA)', (i, i=1, NF)
    write (unit_curr, *) frequency/(2*pi_d*time_unit), (curr(i)*pA, i=NCF+1, NF)
    close (unit_curr)

    return
end subroutine write_cur_out


subroutine write_pop_out(Rho, NCF, Ndim, frequency)
    implicit none
    integer, intent(in) :: NCF, Ndim
    real(q), intent(in) :: frequency
    complex (qc), intent(in) :: Rho(:,:,:)

    integer :: i, l
    character(len=50) :: filename
    
    write(*,*) 'Populations are written in POPULATIONS_'
    write(*,*) ''
    
    do i = -NCF, NCF
        write (filename, '(A13, I0.3, A4)') 'POPULATIONS_n', i, '.dat'
        open  (unit_pop+i, file=filename)
        write (unit_pop+i,*) frequency/(2*pi_d*time_unit), (dble(Rho (l,l,i+NCF+1)), l= 1, Ndim)
        close (unit_pop+i)
    enddo

    return
end subroutine write_pop_out


subroutine write_coh_out(Rho, NCF, frequency)
    implicit none
    integer, intent(in) :: NCF
    real(q), intent(in) :: frequency
    complex(qc), intent(in) :: Rho(:,:,:)
    integer :: i
    character(len=30) :: filename
    
    write(*,*) 'Coherences are written in COHERENCES_'
    write(*,*) ''

    do i = -NCF, NCF
        write (filename, '(A12, I0.3, A4)') 'COHERENCES_n', i, '.dat'
        open  (unit_coh+i, file=filename)
        write (unit_coh+i,*) frequency/(2*pi_d*time_unit), dble(Rho (:,:,i+NCF+1)), dimag(Rho (:,:,i+NCF+1))
        close (unit_coh+i)
    enddo

    return
end subroutine write_coh_out


subroutine write_spin_out (Rho, NCF, Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, hx, hy, hz,&
                           Sx, Sy, Sz, Sh, spin2_ave)
    implicit none
    integer, intent(in) :: NCF, Nm, Ndim
    real(q), intent(in) ::  hx(:), hy(:), hz(:)
    complex(q), intent(in) :: Ss(:,:,:,:), spinX(:,:,:), spinY(:,:,:), spinZ(:,:,:), spin2_T(:,:)
    complex(q), intent(in) :: Rho(:,:,:), H(:,:)
    complex(q), intent(inout) :: spin2_ave
    complex(q), intent(inout), allocatable :: Sx(:), Sy(:), Sz(:), Sh(:)
    real(q) :: frequency
    character(len=30) :: filename
    
    ! Declare local variables

    integer ::  i, l
    
    write(*,*) 'Spin Sx,Sy,Sz,Sh per site are written in SpinFloquet_'
    write(*,*) ''
    
    do i = -NCF, NCF
        
        write(filename, '(A13, I0.3, A4)') 'SpinFloquet_n', i, '.dat' 
        open (unit_floq+i, file=filename)

        call SpinFloquet (Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, Rho(:,:,i+NCF+1), hx, hy, hz,&
                  Sx,Sy,Sz,Sh,spin2_ave)
                
        write (unit_floq+i,*) frequency/(2*pi_d*time_unit),&
            (real(Sx(l)), real(Sy(l)), real(Sz(l)), real(Sh(l)),  l=1, Nm),&
            real(spin2_ave),real(sqrt(1+4*(spin2_ave))-1)*0.5
        close (unit_floq+i)
    enddo

    return
end subroutine write_spin_out

end module output