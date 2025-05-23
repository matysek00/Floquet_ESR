module Transport_F
Use declarations
CONTAINS
     subroutine Current (Ndim, NF, NCF, X, XX, GC, curr, Electrode)
     implicit none
     integer, intent (in) :: Ndim, NF, NCF
     integer :: l, u, j, pn,pnn, i, Electrode
     real (q), intent(out) :: curr(:)
     complex (qc), intent (in) ::  GC (:,:,:,:,:,:), X (:),XX(:)

    curr = 0._q
     
    i=1
    ! Floquet n = 0
    do pn = NF,1,-1

    do l = 1, Ndim
    do u = 1, Ndim
    do j = 1, Ndim
        
        curr(NCF+1) = curr(NCF+1) + 2._q*dble( (GC ( l, j, j, u, pn,1)*(Electrode)&
        & + GC ( l, j, j, u, pn,2)*(1-Electrode) )* X (i))
    enddo
    
        i=i+1
    enddo
    enddo
    enddo
    
    ! Floquet n = 1 for the current, so rho goes from -1 to 3 and conjg from 1 to -3
    
    i=1+ndim*ndim*(NCF-1)

    do pn = NCF+1+2,NCF+1-2,-1
    
    if (i>ndim*ndim*NF)then
        exit
    endif
    
    pnn= NF - pn + 1
    
    do l = 1, Ndim
    do u = 1, Ndim
    do j = 1, Ndim
        curr(NCF+2) = curr(NCF+2) + (GC(l,j,j,u,pn,1)*X(i)+conjg(GC(l,j,j,u,pnn,1)*XX(i)))*Electrode&
        & + (GC(l,j,j,u,pn,2)*X(i)+conjg(GC(l,j,j,u,pnn,2)*XX(i)))*(1-Electrode)
    enddo
        i=i+1
    enddo
    enddo
    enddo
    
    
    ! Floquet n = 2 for the current, so rho goes from 0 to 4 and conjg from -4 to 0

    i=1+ndim*ndim*(NCF)

    do pn = NCF+1+2,NCF+1-2,-1
    
    if (i>ndim*ndim*NF)then
        exit
    endif
    
    pnn= NF - pn + 1
    
    do l = 1, Ndim
    do u = 1, Ndim
    do j = 1, Ndim
        curr(NCF+3) = curr(NCF+3) + (GC(l,j,j,u,pn,1)*X(i)+conjg(GC(l,j,j,u,pnn,1)*XX(i)))*Electrode&
        & + (GC(l,j,j,u,pnn,2)*X(i)+conjg(GC(l,j,j,u,pnn,2)*XX(i)))*(1-Electrode)
    enddo
        i=i+1
    enddo
    enddo
    enddo

    
     return
     end subroutine Current
end module Transport_F
