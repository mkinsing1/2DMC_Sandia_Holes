CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       PERFORM THE FREE-FLIGHT AND SCATTER 
C       PART WITHIN ONE TIME INTERVAL    
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine free_flight_scatter(n_lev,nsim)
       implicit real*8(a-h, o-z)
       
       integer n_lev,nsim,iv,is
       real*8 kx,ky,e
        
       common
     &/ran_var/iso
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/time_1/dt,dtau,tot_time
     &/variables/p(20000,3),ip(20000,2),energy(20000)
     &/particle_atr/kx,ky,iv,is,e
       
      do i = 1, nsim	! loop for all carriers
      !write(*,*) 'nsim',i,nsim

C     Inverse mapping of particle atributes

      kx = p(i,1)
      ky = p(i,2)
      dtau = p(i,3)
      iv = ip(i,1)
      is = ip(i,2)
      e = energy(i)

C     Initial free-flight of the carriers

      dte = dtau
      if(dte.ge.dt)then
         dt2=dt
      else
         dt2=dte
      endif
      
      call drift(dt2)
           
      if(dte.gt.dt)goto 401

      !write(*,*) 'dte.gt.dt',dte,dt
C     Free-flight and scatter part

402   dte2=dte

      call scatter_carrier(n_lev)
                	
219   rr=ran(iso)
      if(rr.le.1e-6) go to 219
      dt3=-(dlog(rr))*tau_max(iv,is)
      dtp = dt - dte2	! remaining time to scatter in dt-interval
      if(dt3.le.dtp)then
	 dt2 = dt3
      else
	 dt2 = dtp
      endif

      call drift(dt2)
	   
c     Update times

      dte2 = dte2 + dt3
      dte = dte2
      if(dte.lt.dt)goto 402
           
401   dte = dte - dt
      dtau = dte

C     Map particle atributes

      p(i,1) = kx
      p(i,2) = ky
      p(i,3) = dtau
      ip(i,1) = iv
      ip(i,2) = is
      energy(i) = e	             

      enddo

      return	
      end
