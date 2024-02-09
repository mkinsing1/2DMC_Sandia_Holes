CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Generic subroutine for the calculation of 
C     INTERVALLEY PHONONS scattering rate
c     (absorption + emission)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine roughness(n_lev,init_valley,ifin_valley)
       implicit real*8 (a-h,o-z)

      integer n_lev
      real*8 kb
           
       common
     &/fund_const/q,h,kb,am0,eps_0 
     &/pi/pi,two_pi
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)
     &/scatt_par2/flag_mech(3,10,50)
     &/scatt_par3/subband(3,10,50),valley(3,10,50)
     &/table/scatt_table(3,10,50,4000)
     &/intervalley2/final_valleys,i_mech
     &/wavefunc/wavefunc(10000,10,3),sub_energy(10,3)
     &/overlap_sr/overlap_sr(3,10,10)
     &/surface_roughness/delta,corr_length
     &/index/k_sub(3)
     &/counter/i_count(3,10)
     &/mas_dos/r_md(3)
     &/sr_flag/sr_flag
     
C    Calculate constants

      final_mass = r_md(ifin_valley)
      const = 0.5D0*final_mass*delta*delta/h*corr_length*
     1        corr_length/h/h*q*q
      
C    (a) Scattering rate - absorption

      do isub = 1, k_sub(init_valley)
      do jsub = 1, k_sub(init_valley)
      
         i_count(init_valley,isub) = 
     1           i_count(init_valley,isub) + 1
     
         sr = const
         
         do i = 1, n_lev
         
         ee = de*dfloat(i)
         
         ef = ee + sub_energy(isub,init_valley) -
     1                  sub_energy(jsub,init_valley)
     
         if(ef.le.0.D0)then         
            rough = 0.D0           
         else         
            sum = 0.D0
            d_fi = two_pi/20.D0
            
            do j = 1,20
               fi = dfloat(j)*d_fi
               q_vec2 = 2.D0*final_mass/h*q/h*
     1            (ee + ef - 2.D0*dsqrt(ee*ef)*dcos(fi))
               screened_sr = overlap_sr(init_valley,isub,jsub)   
               sum = sum + dexp(-q_vec2*corr_length*
     1               corr_length/4.D0)*d_fi*(screened_sr**2)              
            enddo               
            rough = sum*sr
         endif
         
         ii = i_count(init_valley,isub)   
         scatt_table(init_valley,isub,ii,i) = rough
          
         enddo
         
         flag_mech(init_valley,isub,ii) = i_mech
         w(init_valley,isub,ii) = sub_energy(isub,init_valley)-
     1                  sub_energy(jsub,init_valley)
         subband(init_valley,isub,ii) = jsub
         valley(init_valley,isub,ii)  = init_valley
         
      enddo
      enddo
      
      return
      end
