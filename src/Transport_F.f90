module Transport_F
Use declarations
CONTAINS


subroutine Current (Ndim, NF, NCF, Rho, GC, curr)
    implicit none
    integer, intent (in) :: Ndim, NF, NCF
    integer :: l, u, j, pn,pm, i, Electrode
    real (q), intent(out) :: curr (:)
    complex (qc), intent (in) ::  GC (Ndim, Ndim, Ndim,Ndim, NF), Rho (Ndim, Ndim, NF)

    curr = 0._q
    
    
    do l = 1, Ndim
    do u = 1, Ndim
    do j = 1, Ndim
    do pn = -NCF, NCF
!       TODO: these could be combined possibly but maybe this is more readable
!       rho(u,l,-m-n) * GC(l,j,j,u,m) 
        do pm = max(-NCF-pn,-NCF), min(NCF-pn, NCF)
            curr (pn+NCF+1) = curr (pn+(NCF+1)) + Rho (l,u, -pn-pm+(NCF+1)) * GC (l,j,j,u,pm+(NCF+1)) 
        enddo
        
        do pm = max(-NCF+pn,-NCF), min(NCF+pn,NCF)
!           rho*(u,l,n-m) * GC*(l,j,j,u,m) 
            curr (pn+NCF+1) = curr (pn+NCF+1) + conjg(Rho (l,u, pn-pm+(NCF+1)) * GC (l,j,j,u,pm+NCF+1))
        enddo

    enddo
    enddo
    enddo
    enddo
    
    return
end subroutine Current
end module Transport_F
