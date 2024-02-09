CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       INITIALIZE CARRIER ENERGY AND WAVEVECTOR ACCORDING TO THE
C       MAXWELL-BOLTZMANN STATISTICS
C
C       Assumption:  All carriers are initially in valley 1 for GaAs
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
       subroutine init(nsim)
       implicit real*8 (a-h, o-z)

       real*8 k,kx,ky,e,sume,rr
       integer iv,is,nsim
               
       common
     &/ran_var/iso
     &/pi/pi,two_pi  
     &/temp/tem,Vt
     &/mass/smh(3),hhm(3)
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/variables/p(20000,3),ip(20000,2),energy(20000)
     
       sume = 0.D0

       do i = 1, nsim
        
       e = -Vt*log(ran(iso))
       sume = sume + e

c      Initial valley index and region
	
       rr = 3.D0*ran(iso)	
       if(rr.le.1.)then
          iv = 1
          is = 1
       elseif(rr.le.2.)then  
	  iv = 2
	  is = 1
       elseif(rr.le.3.)then
          iv = 3
          is = 1
       endif

C      Initial wavevector

       k=smh(iv)*sqrt(e)
       fai=two_pi*ran(iso)
       kx=k*cos(fai)
       ky=k*sin(fai)

C      Initial free-flight

103    rr = ran(iso)
       if(rr.le.1.e-6)go to 103
       tc = -(dlog(rr))*tau_max(iv,is)

C      Map particle atributes
       
       p(i,1) = kx             
       p(i,2) = ky
       p(i,3) = tc
 
       ip(i,1) = iv
       ip(i,2) = is
       energy(i) = e

       enddo

       sume = sume/dfloat(nsim)
       print*,'Average carrier energy - initial distribution'
       print*,'Energy = ',sume
       print*,'   '
	   
       return
       end
