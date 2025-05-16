module feed_back_on
use declarations
Use OpenFiles !all units of files in open statements are here
Use QME_F !(Quantum Master Equation) contains rates and The Matrix
Use Transport_F ! computation of electron current
CONTAINS
   ! feedbackon Subroutine
   Subroutine feedback(Spin_polarization_L,Spin_polarization_R,gamma_L_0,gamma_R_0,&
   &curr,Iset,tol,bias_L,bias_R,feedbackon,GC,G,A_fast_left,A_fast_right,X,XX,ratio)
   integer :: i_feed
   real(q), intent (in) :: Spin_polarization_L,Spin_polarization_R,Iset,tol,ratio
   real(q), intent (inout) :: gamma_L_0,gamma_R_0
   real(q) :: bias_L,bias_R!,ratio
   logical, intent (inout) :: feedbackon
   complex (qc), intent (inout) :: GC(:,:,:,:,:,:),G(:,:,:,:,:,:),X(:),XX(:)
   real(q), intent (inout) :: curr(:)
   complex (qc), intent (inout) :: A_fast_left(:,:),A_fast_right(:,:)

   i_feed=1
!    ratio=10._q ! max ratio between gammas for feedbackon to stop

do while (feedbackon)
  if ((Spin_polarization_L.ne.0).and.(gamma_L_0.gt.ratio*gamma_R_0))then
   print*,'Tip crashed since',gamma_L_0*hartree,'larger than',gamma_R_0*hartree,'Feedback on stop'
   print*,'Current',curr(NCF+1)*pA,'Iset',Iset

!    do while (feedbackon)
  write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
!    feedbackon=.false.
!    ! We remove the coupling dependence on the rates to use later the update couplings R,L
!   GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)/gamma_R_0
!   GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)/gamma_L_0
!   G (:,:,:,:,:,1)=G (:,:,:,:,:,1)/gamma_R_0
!   G (:,:,:,:,:,2)=G (:,:,:,:,:,2)/gamma_L_0
!   A_fast_left=A_fast_left/gamma_L_0
!   A_fast_right=A_fast_right/gamma_R_0
!
!    ! So, if the ratio condition is achieved, we increase the gammaR. Physically the tip push the molecule more
!    ! and it couples more to the R electrode, unpolarized
!
! gamma_R_0=gamma_R_0+gamma_R_0*(1.0-(1.0-abs(Spin_polarization_L)/(abs(Spin_polarization_L)+1E-10)))*&
! &(tol*10)*sign(1._q,Iset-abs(curr(NCF+1)*pA))*0.5*&
! &(((Iset-tol/10-abs(curr(NCF+1)*pA))/(abs(curr(NCF+1)*pA)+Iset))*(1+(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA))))))+&
! &((abs(curr(NCF+1)*pA)-Iset-tol/10)/(abs(curr(NCF+1)*pA)+Iset))*(1-(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA)))))))
!
!   GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)*gamma_R_0
!   GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)*gamma_L_0
!   G (:,:,:,:,:,1)=G (:,:,:,:,:,1)*gamma_R_0
!   G (:,:,:,:,:,2)=G (:,:,:,:,:,2)*gamma_L_0
!   A_fast_left=A_fast_left*gamma_L_0
!   A_fast_right=A_fast_right*gamma_R_0
!
!   call coeff_matrix (Ndim, omega, NF, NCF, G, Nmatrix, A, B, A_fast_left,A_fast_right, faster)
!
!      UPLO = 'U'
!      NRHS = 1 !one column for B
!      LDA = Nmatrix
!      LDB = Nmatrix
!      LDX = Nmatrix
!
!     call zgesv(Nmatrix,NRHS,A,LDA,IPIV,B,LDB,INFO)
!     ! better solver and similar performance
!     X=B
!     XX=zero
!     do i=1,NF
!         XX(1+ndim*ndim*(NF-i):ndim*ndim*(NF+1-i))=X(1+ndim*ndim*(i-1):ndim*ndim*i)
!     enddo
!
!   ! Compute the DC electron current
!       call Current (Ndim, NF, NCF, X, XX, GC, curr, Electrode)
!   !      call clock ('STEP 5:: Finished DC current ', 2)
!   i_feed=i_feed+1
! !   print*,curr(NCF+1)*pA, gamma_L_0*hartree,gamma_R_0*hartree,i_feed
!     if ((Spin_polarization_L.ne.0).and.(gamma_R_0.gt.gamma_L_0))then
!         feedbackon=.false.
!     print*,'Cannot reach Iset even after increasing unpolarized coupling',&
!     &gamma_L_0*hartree,'equals to',gamma_R_0*hartree,'Feedback on stop'
!     print*,'Current',curr(NCF+1)*pA,'Iset',Iset
!     elseif (((abs(curr(NCF+1)*pA).gt.(Iset-tol)).and.((abs(curr(NCF+1)*pA)).lt.(Iset+tol))).and.(i_seed.lt.1E5))then
!       feedbackon=.false.
!       print*,'Condition for Iset sucessfully achieved'
!       print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
!       write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
!     elseif (i_feed.gt.1E5)then
!       print*,'Max number of iterations reached, please check if current',curr(NCF+1)*pA,'is close to Iset',Iset
!       print*,'otherwise, increase/decrease the polarized coupling to improve convergence.'
!       print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
!       write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
      feedbackon=.false.
!     endif
!     enddo

  elseif ((Spin_polarization_R.ne.0).and.(gamma_R_0.gt.ratio*gamma_L_0))then
   print*,'Tip crashed since',gamma_R_0*hartree,'larger than',gamma_L_0*hartree,'Feedback on stop'
   print*,'Current',curr(NCF+1)*pA,'Iset',Iset

!    do while (feedbackon)
  write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
   !    feedbackon=.false.
!    ! We remove the coupling dependence on the rates to use later the update couplings R,L
!   GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)/gamma_R_0
!   GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)/gamma_L_0
!   G (:,:,:,:,:,1)=G (:,:,:,:,:,1)/gamma_R_0
!   G (:,:,:,:,:,2)=G (:,:,:,:,:,2)/gamma_L_0
!   A_fast_left=A_fast_left/gamma_L_0
!   A_fast_right=A_fast_right/gamma_R_0
!
!    ! So, if the ratio condition is achieved, we increase the gammaR. Physically the tip push the molecule more
!    ! and it couples more to the R electrode, unpolarized
!
! gamma_L_0=gamma_L_0+gamma_L_0*(1.0-(1.0-abs(Spin_polarization_R)/(abs(Spin_polarization_R)+1E-10)))*&
! &(tol*10)*sign(1._q,Iset-abs(curr(NCF+1)*pA))*0.5*&
! &(((Iset-tol/10-abs(curr(NCF+1)*pA))/(abs(curr(NCF+1)*pA)+Iset))*(1+(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA))))))+&
! &((abs(curr(NCF+1)*pA)-Iset-tol/10)/(abs(curr(NCF+1)*pA)+Iset))*(1-(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA)))))))
!
!   GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)*gamma_R_0
!   GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)*gamma_L_0
!   G (:,:,:,:,:,1)=G (:,:,:,:,:,1)*gamma_R_0
!   G (:,:,:,:,:,2)=G (:,:,:,:,:,2)*gamma_L_0
!   A_fast_left=A_fast_left*gamma_L_0
!   A_fast_right=A_fast_right*gamma_R_0
!
!   call coeff_matrix (Ndim, omega, NF, NCF, G, Nmatrix, A, B, A_fast_left,A_fast_right, faster)
!
!      UPLO = 'U'
!      NRHS = 1 !one column for B
!      LDA = Nmatrix
!      LDB = Nmatrix
!      LDX = Nmatrix
!
!     call zgesv(Nmatrix,NRHS,A,LDA,IPIV,B,LDB,INFO)
!     ! better solver and similar performance
!     X=B
!     XX=zero
!     do i=1,NF
!         XX(1+ndim*ndim*(NF-i):ndim*ndim*(NF+1-i))=X(1+ndim*ndim*(i-1):ndim*ndim*i)
!     enddo
!
!   ! Compute the DC electron current
!       call Current (Ndim, NF, NCF, X, XX, GC, curr, Electrode)
!   !      call clock ('STEP 5:: Finished DC current ', 2)
!   i_feed=i_feed+1
! !   print*,curr(NCF+1)*pA, gamma_L_0*hartree,gamma_R_0*hartree,i_feed
!     if ((Spin_polarization_R.ne.0).and.(gamma_L_0.gt.gamma_R_0))then
!         feedbackon=.false.
!     print*,'Cannot reach Iset even after increasing unpolarized coupling',&
!     &gamma_R_0*hartree,'equals to',gamma_L_0*hartree,'Feedback on stop'
!     print*,'Current',curr(NCF+1)*pA,'Iset',Iset
!     elseif (((abs(curr(NCF+1)*pA).gt.(Iset-tol)).and.((abs(curr(NCF+1)*pA)).lt.(Iset+tol))).and.(i_seed.lt.1E5))then
!       feedbackon=.false.
!       print*,'Condition for Iset sucessfully achieved'
!       print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
!       write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
!     elseif (i_feed.gt.1E5)then
!     print*,'Max number of iterations reached, please check if current',curr(NCF+1)*pA,'is close to Iset',Iset
!     print*,'otherwise, increase/decrease the polarized coupling to improve convergence.'
!     print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
!   write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
    feedbackon=.false.
!     endif
!     enddo

  elseif (((abs(curr(NCF+1)*pA).lt.(Iset-tol)).or.((abs(curr(NCF+1)*pA)).gt.(Iset+tol))).and.(i_feed.le.1E5))then
  !.and.(VDC.ne.0._q))then
  write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree

  ! We remove the coupling dependence on the rates to use later the update couplings R,L
  GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)/gamma_R_0
  GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)/gamma_L_0
  G (:,:,:,:,:,1)=G (:,:,:,:,:,1)/gamma_R_0
  G (:,:,:,:,:,2)=G (:,:,:,:,:,2)/gamma_L_0
  A_fast_left=A_fast_left/gamma_L_0
  A_fast_right=A_fast_right/gamma_R_0

gamma_L_0=gamma_L_0+gamma_L_0*(1.0-(1.0-abs(Spin_polarization_L)/(abs(Spin_polarization_L)+1E-10)))*&
&(tol*10)*sign(1._q,Iset-abs(curr(NCF+1)*pA))*0.5*&
&(((Iset-tol/10-abs(curr(NCF+1)*pA))/(abs(curr(NCF+1)*pA)+Iset))*(1+(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA))))))+&
&((abs(curr(NCF+1)*pA)-Iset-tol/10)/(abs(curr(NCF+1)*pA)+Iset))*(1-(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA)))))))

gamma_R_0=gamma_R_0+gamma_R_0*(1.0-(1.0-abs(Spin_polarization_R)/(abs(Spin_polarization_R)+1E-10)))*&
&(tol*10)*sign(1._q,Iset-abs(curr(NCF+1)*pA))*0.5*&
&(((Iset-tol/10-abs(curr(NCF+1)*pA))/(abs(curr(NCF+1)*pA)+Iset))*(1+(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA))))))+&
&((abs(curr(NCF+1)*pA)-Iset-tol/10)/(abs(curr(NCF+1)*pA)+Iset))*(1-(-1)**(dble(int(1.-0.1*sign(1._q,Iset-abs(curr(NCF+1)*pA)))))))

  GC (:,:,:,:,:,1)=GC (:,:,:,:,:,1)*gamma_R_0
  GC (:,:,:,:,:,2)=GC (:,:,:,:,:,2)*gamma_L_0
  G (:,:,:,:,:,1)=G (:,:,:,:,:,1)*gamma_R_0
  G (:,:,:,:,:,2)=G (:,:,:,:,:,2)*gamma_L_0
  A_fast_left=A_fast_left*gamma_L_0
  A_fast_right=A_fast_right*gamma_R_0

  call coeff_matrix (Ndim, omega, NF, NCF, G, Nmatrix, A, B, A_fast_left,A_fast_right, faster)

     UPLO = 'U'
     NRHS = 1 !one column for B
     LDA = Nmatrix
     LDB = Nmatrix
     LDX = Nmatrix

    call zgesv(Nmatrix,NRHS,A,LDA,IPIV,B,LDB,INFO)
    ! better solver and similar performance
    X=B
    XX=zero
    do i=1,NF
        XX(1+ndim*ndim*(NF-i):ndim*ndim*(NF+1-i))=X(1+ndim*ndim*(i-1):ndim*ndim*i)
    enddo

  ! Compute the DC electron current
      call Current (Ndim, NF, NCF, X, XX, GC, curr, Electrode)
  !      call clock ('STEP 5:: Finished DC current ', 2)
  i_feed=i_feed+1
!   print*,curr(NCF+1)*pA, gamma_L_0*hartree,gamma_R_0*hartree,i_feed

  elseif (i_feed.gt.1E5)then
    print*,'Max number of iterations reached, please check if current',curr(NCF+1)*pA,'is close to Iset',Iset
    print*,'otherwise, increase/decrease the polarized coupling to improve convergence.'
    print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
  write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
    feedbackon=.false.
  else
   print*,'Condition for Iset sucessfully achieved'
   print*,'New coupling:', gamma_L_0*hartree,gamma_R_0*hartree
  write(unit_feedbackon,*) curr(NCF+1)*pA,bias_L*hartree,bias_R*hartree, gamma_L_0*hartree,gamma_R_0*hartree
   feedbackon=.false.
  endif
enddo

   end Subroutine feedback
end module  feed_back_on
