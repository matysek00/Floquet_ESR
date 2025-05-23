module Matrix_Coeff
    Use declarations
    CONTAINS

subroutine coeff_matrix2 (Ndim, NF, NCF, GA, QME_Tensor_A)
    implicit none
    integer, intent (in) :: Ndim, NF, NCF
    complex (qc), intent (in) :: GA (Ndim, Ndim, Ndim, Ndim, NF)

    complex (qc), intent(out), allocatable :: QME_Tensor_A (:,:,:,:,:,:)
    integer :: l, j, v, u, n, m 
            
    allocate (QME_Tensor_A (Ndim, Ndim, NF, Ndim, Ndim, NF))
    !allocate (QME_Tensor_A (NF, Ndim, Ndim, NF, Ndim, Ndim))
    QME_Tensor_A = zero
    
    fourier_n: do n= 1, NF
    fourier_m: do m = max(1,n-NCF),min(NF,n+NCF)
    level_l: do l=1, Ndim
    level_j: do j=1, Ndim
    level_v: do v=1, Ndim
    level_u: do u=1, Ndim
!               Gamma_{j,v,v,u:n-m} rho_{l,u,m}
            QME_Tensor_A (j,l,n, u,l,m) = QME_Tensor_A (j,l,n, u,l,m) + ui*GA (j,v,v,u,n-m+NCF+1)

!               Gamma_{l,v,v,u:m-n} rho_{u,j,m}
            QME_Tensor_A (j,l,n, j,u,m) = QME_Tensor_A (j,l,n, j,u,m) + ui*CONJG(GA (l,v,v,u,m-n+NCF+1))

!               Gamma_{u,j,l,v:m-n} rho_{v,u,m}
            QME_Tensor_A (j,l,n, u,v,m) = QME_Tensor_A (j,l,n, u,v,m) - ui*CONJG(GA (u,j,l,v,m-n+NCF+1))
            !
!               !Gamma_{v,l,j,u:n-m} rho_{v,u,m}
            QME_Tensor_A (j,l,n, u,v,m) = QME_Tensor_A (j,l,n, u,v,m) - ui*(GA (v,l,j,u,n-m+NCF+1))
            
    enddo level_u
    enddo level_v
    enddo level_j
    enddo level_l
    enddo fourier_m
    enddo fourier_n
    return 
end subroutine coeff_matrix2


subroutine coeff_matrix (Ndim, frequency, NF, NCF, G, Nmatrix, A, B)
    implicit none
    integer,intent (in) :: Ndim, NF, NCF
    integer,intent (in) :: Nmatrix
    real (q), intent(in) :: frequency
    complex (qc), intent (in) :: G (Ndim, Ndim, Ndim, Ndim, NF, 2)
    complex (qc), intent (inout) :: A (Nmatrix, Nmatrix)
    complex (qc), intent (out) :: B (Nmatrix)
    
    complex (qc), allocatable :: QME_Tensor_R (:,:,:,:,:,:), QME_Tensor_L (:,:,:,:,:,:)
    complex (qc), allocatable :: QME_Tensor(:,:,:,:,:,:), QME_Tensor_D(:,:,:,:,:,:)

    integer :: l, j, n, i2

!       Return A that solves rho_{l,j,n} A_{l,j,n u,v,m}  = delta(v,1) * delta(u,1) * delta (n,0)
!       B is the right hand side of the equation 
!       B = delta(v,1) * delta(u,1) * delta (n,0)

    allocate (QME_Tensor (Ndim, Ndim, NF, Ndim, Ndim, NF))
    allocate (QME_Tensor_D (Ndim, Ndim, NF, Ndim, Ndim, NF))
    !allocate (QME_Tensor (NF, Ndim, Ndim, NF, Ndim, Ndim))
    !allocate (QME_Tensor_D (NF, Ndim, Ndim, NF, Ndim, Ndim))
!       In dot zeroth-order dynamics 

    QME_Tensor_D = zero
    QME_Tensor = zero
    do n=1,NF
    do l = 1, Ndim
    do j = 1, Ndim
            QME_Tensor_D (j,l,n, j,l,n) = Delta (l,j) + frequency*(n-NCF-1)
    enddo
    enddo
    enddo

!       Right and left electrode hopping
    call coeff_matrix2 (Ndim, NF, NCF, G(:,:,:,:,:,1), QME_Tensor_R) 
    call coeff_matrix2 (Ndim, NF, NCF, G(:,:,:,:,:,2), QME_Tensor_L) 
    
    QME_Tensor = QME_Tensor_D + QME_Tensor_L + QME_Tensor_R      

!       I don't understand why we do this
    do n=1,NF
            QME_Tensor(1,1,n, :,:,:) = zero
    do l = 1, Ndim
            QME_Tensor(1,1,n, l,l,n) = one
    enddo  
    enddo
    
    A = reshape(QME_Tensor, (/Nmatrix, Nmatrix/))
    i2 = 0
    
!       B definition entails detailed balance
    B = zero
    B (1+ndim*ndim*(NCF)) = one
    
    return
end subroutine coeff_matrix 
end module Matrix_Coeff