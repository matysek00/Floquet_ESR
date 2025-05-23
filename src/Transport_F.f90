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
    
    do pn = -NCF, NCF
    do l = 1, Ndim
    do u = 1, Ndim
    do j = 1, Ndim

!       TODO: these could be combined
        do pm = max(pn-NCF,1), min(pn+NCF, NCF)
            curr (pn+NCF+1) = Rho (l,u, NCF-pn-pm+1) * GC (l,j,j,u,pm+NCF+1) 
        enddo
        
        do pm = max(pn+NCF,1), min(pn-NCF,1)
            curr (pn+NCF+1) = curr (pn+NCF+1) + conjg(Rho (l,u, NCF+pn-pm+1) * GC (l,j,j,u,pm+NCF+1))
        enddo

    enddo
    enddo
    enddo
    enddo
    
     return
     end subroutine Current
end module Transport_F
