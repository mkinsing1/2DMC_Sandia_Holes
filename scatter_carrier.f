CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       SELECT SCATTERING MECHANISM AND PERFORM 
C       THE SCATTERING PART THAT MODIFIES PARTICLE ATRIBUTES
C       (kx, ky, kz, iv, energy, i_region)    
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine scatter_carrier(n_lev)
       implicit real*8 (a-h, o-z)

       integer n_lev,iv,is,loc
       real*8 kx,ky,e

       common
     &/ran_var/iso
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par2/flag_mech(3,10,50)
     &/table/scatt_table(3,10,50,4000)
     &/particle_atr/kx,ky,iv,is,e
     &/counter/i_count(3,10)

C    Calculate index to the scattering table

      loc = int(e/de)
      if(loc.eq.0)loc=1
      if(loc.gt.n_lev)loc=n_lev
	   
C    Select scattering mechanism

      i_top = i_count(iv,is)
      rr = ran(iso)
      if(rr.ge.scatt_table(iv,is,i_top,loc))then
         goto 222  ! self-scattering
      endif
      if(rr.lt.scatt_table(iv,is,1,loc))then
         i_fix = 1
         goto 111
      endif
      if(i_top.gt.1)then
        do i=1,i_top-1
           bound_lower = scatt_table(iv,is,i,loc)
           bound_upper = scatt_table(iv,is,i+1,loc)
           if(rr.ge.bound_lower.and.rr.lt.bound_upper)then
              i_fix = i + 1
              goto 111
           endif
        enddo
      endif

111   continue

C    Perform scattering (change energy and randomize momentum)    

      select_mech = flag_mech(iv,is,i_fix)
      if(select_mech.eq.1)then
         call isotropic(i_fix)
      elseif(select_mech.eq.2)then
         call isotropic(i_fix)
      elseif(select_mech.eq.3)then
         call roughness_angle(i_fix)
      endif

222   continue

      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       SCATTERING SUBROUTINES
C       (CHANGE ENERGY AND WAVEVECTORS OF PARTICLES)
C
C       In the definition of the scattering rates, a variable called
C       'flag_mech' has been defined.  The values assigned to this
C       variable correspond to:
C       
C       1  ==>  Isotropic scattering (acoustic, intervalley g-phonons)
C       2  ==>  Isotropic scattering (intervalley f-phonons)
C       3  ==>  Coulomb scattering ==> small-angle scattering
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       ISOTROPIC SCATTERING PROCESS
C       uniform probability density for scattering in all directions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine isotropic(i_fix)
       implicit real*8 (a-h,o-z)
       
       integer iv,is,iv1,is1,i_fix
       real*8 kx,ky,e
       
       common
     &/ran_var/iso
     &/pi/pi,two_pi 
     &/mass/smh(3),hhm(3)
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par3/subband(3,10,50),valley(3,10,50)
     &/particle_atr/kx,ky,iv,is,e    
      
C     Update carrier energy
      e = e + w(iv,is,i_fix)
      if(e.lt.0)then
         e = 0.D0
      endif

C     Update carrier wavevector
      rknew = smh(iv)*dsqrt(e)
      fi = two_pi*ran(iso)
      kx = rknew*dcos(fi)
      ky = rknew*dsin(fi)
      iv1 = int(valley(iv,is,i_fix))
      is1 = int(subband(iv,is,i_fix))
c      print*,'Phonons'
c      print*,kx,ky,e,rknew
      iv = iv1
      is = is1

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SURFACE-ROUGHNESS SCATTERING MECHANISM
C       non-uniform probability density function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine roughness_angle(i_fix)
       implicit real*8 (a-h,o-z)
       
       integer iv,is,iv1,is1,i_fix
       real*8 kx,ky,e,ef,kb
       
       common
     &/ran_var/iso
     &/pi/pi,two_pi
     &/fund_const/q,h,kb,am0,eps_0 
     &/mass/smh(3),hhm(3)
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par3/subband(3,10,50),valley(3,10,50)
     &/particle_atr/kx,ky,iv,is,e
     &/mas_dos/r_md(3)
     &/surface_roughness/delta,corr_length    
      
C     Update carrier energy
      ef = e + w(iv,is,i_fix)
      if(ef.lt.0)then
         ef = 0.D0
      endif

      final_mass = r_md(iv)
      factor = final_mass/h*corr_length*corr_length/h*q*
     1         dsqrt(e*ef)     

11    fi1 = ran(iso)
      fi2 = two_pi*ran(iso)      
      C = dexp(factor)
      f1 = fi1*C
      f2 = dexp(-factor*dcos(fi2))      
      if(f1.gt.f2)goto 11
c      write(12,*)fi2,C
      
C     Update carrier wavevector
      rknew = smh(iv)*dsqrt(ef)
      if(kx.gt.0.D0)then
         fi_0 = datan(ky/kx)
      else
         fi_0 = pi + datan(ky/kx)
      endif
      fi_prime = fi_0 - fi2            
      kx = rknew*dcos(fi_prime)
      ky = rknew*dsin(fi_prime)
      
      iv1 = int(valley(iv,is,i_fix))
      is1 = int(subband(iv,is,i_fix))
      iv = iv1
      is = is1
      e = ef

      return
      end
      
