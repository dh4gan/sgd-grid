! **************************************************************
!
!  Interpolates opacity table.
!
! **************************************************************

       subroutine get_kappa(dens,kappa,T)

	use eosdata

       real :: dens, kappa
       real :: upart1, upart2,T
       real :: m, d
       
       integer ::  i, j       
 

       if (dens<1.0e-24) dens = 1.0e-24

       i = 1
       do while((eostable(i,1,1)<=dens).and.(i<nrhopoints))
         i = i + 1
       enddo

       If (T<1.0) T = 1.0

       j = 1
       do while ((eostable(i-1,j,2)<=T).and.(j<nUpoints))
         j = j + 1
       enddo

       m = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
          (eostable(i-1,j-1,2) - eostable(i-1,j,2))
       d = eostable(i-1,j,6) - m*eostable(i-1,j,2)

       upart1 = m*T + d
 
       j = 1
       do while ((eostable(i,j,2)<=T).and.(j<1000))
         j = j + 1
       enddo

       m = (eostable(i,j-1,6) - eostable(i,j,6))/ &
          (eostable(i,j-1,2) - eostable(i,j,2))
       d = eostable(i,j,6) - m*eostable(i,j,2)
   
       upart2 = m*T + d

       m = (upart2 - upart1)/  &
          (eostable(i,1,1)-eostable(i-1,1,1))
       d = upart2 - m*eostable(i,1,1)
 
       kappa = m*dens + d
			       
       return

       end subroutine get_kappa
