CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Generic subroutine for the calculation of 
C     INTERVALLEY PHONONS scattering rate
c     (absorption + emission)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine acoustic_scatt(n_lev,init_valley,ifin_valley)      
       implicit real*8 (a-h,o-z)

      integer n_lev
      real*8 kb, const, de, final_valleys
           
       common 
     &/temp/tem,Vt
     &/fund_const/q,h,kb,am0,eps_0 
     &/nonp/af,af2,af4
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par2/flag_mech(3,10,50)
     &/scatt_par3/subband(3,10,50)
     &/scatt_par3/valley(3,10,50)
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

      final_mass = r_md(ifin_valley)
      const = final_mass*tem*kb/h*(coupling_constant**2)
     1        *final_valleys*q/h*q/h/density/sound_velocity/
     1        sound_velocity
      
      do isub = 1, k_sub(init_valley)
      do jsub = 1, k_sub(ifin_valley)
         i_count(init_valley,isub) = 
     1           i_count(init_valley,isub) + 1
         ac   = const*overlap_intervalley(init_valley,isub,
     1        ifin_valley,jsub)
      
         do i = 1, n_lev
         
         ee = de*dfloat(i)
         ef = ee  + sub_energy(isub,init_valley)-
     1              sub_energy(jsub,ifin_valley)
 
         if(ef.le.0)then
            acoustic_rate = 0.D0
         else
            acoustic_rate = ac
         endif
         
         ii = i_count(init_valley,isub)   
         scatt_table(init_valley,isub,ii,i) = acoustic_rate
         if (init_valley.eq.1) then
           if (isub.eq.1) then
             write(85,*) scatt_table(1,1,ii,i)
           endif
         endif
         
         enddo
         
         flag_mech(init_valley,isub,ii) = i_mech
         w(init_valley,isub,ii) =  
     1                  sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,ifin_valley)
         subband(init_valley,isub,ii) = jsub
         valley(init_valley,isub,ii)  = ifin_valley
         
      enddo
      enddo

      return
      end
