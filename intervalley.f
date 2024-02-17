CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Generic subroutine for the calculation of 
C     INTERVALLEY PHONONS scattering rate
c     (absorption + emission)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine intervalley(n_lev,w0,
     1            init_valley,ifin_valley)
       implicit real*8 (a-h,o-z)

      integer n_lev
      real*8 kb, w0,rnq
           
       common
     &/temp/tem,Vt
     &/fund_const/a0,q,h,hbar,kb,am0,eps_0 
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par2/flag_mech(3,10,50)
     &/scatt_par3/subband(3,10,50),valley(3,10,50)
     &/table/scatt_table(3,10,50,4000)
     &/density/density,sound_velocity
     &/intervalley1/coupling_constant
     &/intervalley2/final_valleys,i_mech
     &/wavefunc/wavefunc(10000,10,3),sub_energy(10,3)
     &/overlap_intervalley/overlap_intervalley(3,10,3,10)
     &/index/k_sub(3)
     &/counter/i_count(3,10)
     &/mas_dos/r_md(3)
     
C    Calculate constants

      rnq = 1.D0/(exp(w0/Vt)-1.D0)
      final_mass = r_md(ifin_valley)
      const = final_mass/2.D0/h*(coupling_constant**2)
     1        /w0*final_valleys*q/density
      
C    (a) Scattering rate - absorption

      do isub = 1, k_sub(init_valley)
      do jsub = 1, k_sub(ifin_valley)
      
         i_count(init_valley,isub) = 
     1           i_count(init_valley,isub) + 1
         ab = rnq*const*
     1        overlap_intervalley(init_valley,isub,
     2        ifin_valley,jsub) 
         do i = 1, n_lev
         ee = de*dfloat(i)
         ef = ee + w0 + sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,ifin_valley)
         if(ef.le.0.D0)then
            absorption = 0.D0
         else
            absorption = ab
         endif
         ii = i_count(init_valley,isub)   
         scatt_table(init_valley,isub,ii,i) = absorption 
         enddo
         flag_mech(init_valley,isub,ii) = i_mech
         w(init_valley,isub,ii) = w0 + 
     1                  sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,ifin_valley)
         subband(init_valley,isub,ii) = jsub
         valley(init_valley,isub,ii)  = ifin_valley
         
      enddo
      enddo

C    (b) Scattering rate - emission

      do isub = 1, k_sub(init_valley)
      do jsub = 1, k_sub(ifin_valley)
         i_count(init_valley,isub) = 
     1           i_count(init_valley,isub) + 1
         em = (rnq+1)*const*
     1        overlap_intervalley(init_valley,isub,
     2        ifin_valley,jsub) 
         do i = 1, n_lev
         ee = de*dfloat(i)
         ef = ee - w0 + sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,ifin_valley)
         if(ef.le.0)then
            emission = 0.D0
         else
            emission = em
         endif
         ii = i_count(init_valley,isub)  
         scatt_table(init_valley,isub,ii,i) = emission 
         enddo
         flag_mech(init_valley,isub,ii) = i_mech
         w(init_valley,isub,ii) = - w0 + 
     1                  sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,ifin_valley)
         subband(init_valley,isub,ii) = jsub
         valley(init_valley,isub,ii)  = ifin_valley
         
      enddo
      enddo
      
      return
      end
