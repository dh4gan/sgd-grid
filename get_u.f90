! **************************************************************
!
!  Interpolates opacity table.
!
! **************************************************************

       subroutine get_u(dens,intenerg,T)

	use sphNGdata
	use serendata
	use eosdata

       real :: dens, intenerg
       real :: upart1, upart2,T
       real :: m, d
       
       integer ::  i, j       
 

       if (dens<1.0e-24) dens = 1.0e-24

       i = 1
       do while((eostable(i,1,1)<=dens).and.(i<260))
         i = i + 1
       enddo

       If (T<1.0) T = 1.0

       j = 1
       do while ((eostable(i-1,j,2)<=T).and.(j<1000))
         j = j + 1
       enddo

       m = (eostable(i-1,j-1,3) - eostable(i-1,j,3))/ &
          (eostable(i-1,j-1,2) - eostable(i-1,j,2))
       d = eostable(i-1,j,3) - m*eostable(i-1,j,2)

       upart1 = m*T + d
 
       j = 1
       do while ((eostable(i,j,2)<=T).and.(j<1000))
         j = j + 1
       enddo

       m = (eostable(i,j-1,3) - eostable(i,j,3))/ &
          (eostable(i,j-1,2) - eostable(i,j,2))
       d = eostable(i,j,3) - m*eostable(i,j,2)
   
       upart2 = m*T + d

       m = (upart2 - upart1)/  &
          (eostable(i,1,1)-eostable(i-1,1,1))
       d = upart2 - m*eostable(i,1,1)
 
       intenerg = m*dens + d
			       
       return

       end subroutine get_u
