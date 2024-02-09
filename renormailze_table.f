CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Generic subroutine that renormalizes the scattering table
C     for a given valley
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine renormalize_table(n_lev)
      implicit real*8 (a-h,o-z)

      integer n_lev
      integer ii
      real*8 tau(3,10)

       common
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/table/scatt_table(3,10,50,4000)
     &/counter/i_count(3,10)
     &/index/k_sub(3)
  
      do ival = 1,3
      do isub = 1,k_sub(ival)

      ii = i_count(ival,isub)

         if(ii.gt.1)then
         
            do i = 2, ii
            do k = 1, n_lev
            scatt_table(ival,isub,i,k)=scatt_table(ival,isub,i-1,k)
     1                            + scatt_table(ival,isub,i,k)
            enddo
            enddo
            
         endif
      
         if(ival.eq.1.and.isub.eq.1)then
         do k = 1,n_lev
            open(unit=1,file='scatt_table_unprimed',status='unknown')
            write(1,*)k,scatt_table(ival,isub,ii,k)
         enddo
         close(1)
         endif

          tau(ival,isub) = 0.D0
          do i = 1,n_lev
             if(scatt_table(ival,isub,ii,i).gt.tau(ival,isub))
     1               tau(ival,isub) = scatt_table(ival,isub,ii,i)   
          enddo

          do i = 1, ii
          do k = 1, n_lev
            scatt_table(ival,isub,i,k) = 
     1                  scatt_table(ival,isub,i,k)/tau(ival,isub)
         enddo
         enddo

         tau_max(ival,isub) = 1.D0/tau(ival,isub)
c         print*,'Maximum gama'
c         print*,ival,isub,tau_max(ival,isub)
      
      enddo
      enddo

      return
      end
