CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      MAIN PROGRAM FOR MONTE CARLO SOLUTION OF THE BOLTZMANN
C      TRANSPORT EQUATION FOR A THREE-VALLEY SEMICONDUCTOR
C     
       parameter(n_lev=1000, nele=10000)
C
C      n_lev  => # of energy levels in the scattering table       
C      nele   => # total number of electrons that are simulated
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEVER NEVER NEVER USE implicit
       implicit real*8 (a-h, o-z)
C       implicit none      

       character*40 file_name
       integer nsim
C       real*8 time,flag_write
C       integer j
       
       common
     &/time_1/dt,dtau,tot_time

C     Read parameters from input file
       call readin(n_lev)

C     Calculate overlap integrals
      call read_wavefunctions()

C     Calculate scattering table
       call sc_table(n_lev)

C     Renormalize table
      call renormalize_table(n_lev)

C     Initialize carriers
       nsim = nele
       call init(nsim)
       file_name = 'initial_distribution'
       call histograms(nsim,file_name)     

C     Start the Monte carlo procedure       
       time = 0.D0
       j = 0
       flag_write = 0.D0

       do while(time.le.tot_time)
        
          j = j + 1          
          time=dt*dfloat(j)
          print*,time                  
          call free_flight_scatter(n_lev,nsim)
          file_name = 'current_distribution'
          call histograms(nsim,file_name)
          call write(nsim,j,time,flag_write)
          flag_write = 1.D0
                                   
       enddo	! End of the time loop
       
       end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SAVE HISTOGRAMS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine histograms(nsim,file_name)
       implicit real*8 (a-h, o-z)
C       implicit none

       integer nsim
       real*8 kx,ky,e
       character*30 file_name
       
       common
     &/variables/p(20000,3),ip(20000,2),energy(20000)


       open(unit=6,file=file_name,status='unknown')
       write(6,*)'kx   ky   energy'

       do i = 1, nsim
          kx = p(i,1)
	  ky = p(i,2)
	  e = energy(i)	      		 
          write(6,*)kx,ky,e
       enddo
       close(6)

       return
       end 
       
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      WRITE AVERAGES IN FILES
C      (these averages correspond to one time step)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine write(nsim,iter,time,flag_write)
       implicit real*8 (a-h,o-z)
       
       integer nsim,iter,iv,is
       real*8 kb,time,kx,ky
       character*40 file_1,file_2
       character*40 file_3,file_4,file_5
       integer i_nele(3,10)
       real*8 velocity_sumx(3)
       real*8 velocity_avx(3)
       real*8 velocity_sumy(3)
       real*8 velocity_avy(3)

       integer i_valley(3)
       real*8 energy_sum(3)
       real*8 energy_av(3)
       
       common
     &/fund_const/a0,q,h,hbar,kb,am0,eps_0
     &/variables/p(20000,3),ip(20000,2),energy(20000)
     &/mas_con/r_mc(3)
     &/mas_dos/r_md(3)
     &/index/k_sub(3)
     

C     Open output files

      file_1 = 'vx_time_averages'
      file_2 = 'vy_time_averages'
      file_3 = 'energy_time_averages'
      file_4 = 'unprimed_valley_occupation'
      file_5 = 'primed_valley_occupation'
      open(unit=1,file=file_1,status='unknown')
      open(unit=2,file=file_2,status='unknown')
      open(unit=3,file=file_3,status='unknown')   
      open(unit=4,file=file_4,status='unknown')
      open(unit=5,file=file_5,status='unknown')
      
       do iv = 1,3
          i_valley(iv) = 0
          velocity_sumx(iv) = 0.D0
          velocity_sumy(iv) = 0.D0
          energy_sum(iv) = 0.D0
          
          do isub = 1,k_sub(iv)
             i_nele(iv,isub) = 0
          enddo
          
       enddo

       do i = 1,nsim
      
         iv = ip(i,1)
         is = ip(i,2)
         ee = energy(i)
         kx = p(i,1)
         ky = p(i,2)
         i_nele(iv,is) = i_nele(iv,is) + 1
         i_valley(iv) = i_valley(iv) + 1
         velocity_sumx(iv) = h/r_mc(iv)*kx + velocity_sumx(iv)
         velocity_sumy(iv) = h/r_mc(iv)*ky + velocity_sumy(iv)
c         velocity_sumx(iv) = kx*kx + velocity_sumx(iv)
c         velocity_sumy(iv) = ky*ky + velocity_sumy(iv)         
c         energy_sum(iv) = energy_sum(iv) + 2.D0*r_mc(iv)*ee*q/h/h
         energy_sum(iv) = energy_sum(iv) + ee

       write(10,*)kx
       write(11,*)ky
       
       enddo
       
       close(10)
       close(11) 
       
       do iv = 1,3
           velocity_avx(iv) = velocity_sumx(iv)/dfloat(i_valley(iv))
           velocity_avy(iv) = velocity_sumy(iv)/dfloat(i_valley(iv))
           energy_av(iv) = energy_sum(iv)/dfloat(i_valley(iv))
       enddo
       
       velocity_finalx = 0.D0
       velocity_finaly = 0.D0
       energy_final  = 0.D0
       do iv = 1,3
          velocity_finalx = velocity_finalx + velocity_avx(iv)/3.D0
          velocity_finaly = velocity_finaly + velocity_avy(iv)/3.D0
          energy_final = energy_final + energy_av(iv)/3.D0
       enddo
       
       write(1,*)time, velocity_finalx
       write(2,*)time, velocity_finaly
       write(3,*)time, energy_final
       write(4,*)time, i_nele(1,1)
       write(5,*)time, i_nele(2,1)
     
c       do iv = 1,3
c       do isub = 1,k_sub(iv)
c          print*,iv,isub,i_nele(iv,isub)
c       enddo
c       enddo   

       return
       end
                    		

