CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       SUBROUTINE THAT CREATES THE SCATTERING TABLE
C       flag_mech = 1 ==> isotropic scattering process
C       flag_mech = 2 ==> polar optical phonons
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine sc_table(n_lev)
       implicit real*8 (a-h, o-z)

       integer n_lev, i_mech
       
       common    
     &/select_acouctic/acoustic_scattering
C     &/select_intervalley_1/intervalley_zero_g
     &/select_intervalley_2/intervalley_zero_f
     &/select_sr/surface_roughness
     &/select_alloy/alloy_disorder
     &/select_coulomb/coulomb_scattering
     &/sigma_acoustic/sigma_acoustic
     &/pot_alloy/pot_alloy
     &/coulomb/doping_density(5),Energy_debye(5)
     &/Def_pot_1/DefPot_zero_g
     &/Def_pot_2/DefPot_zero_f
     &/interval_phonons_1/phonon_zero_g
     &/interval_phonons_2/phonon_zero_f
     &/intervalley1/coupling_constant
     &/intervalley2/final_valleys,i_mech

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CREATE THE SCATTERING TABLE FOR EACH REGION
C     Scattering mechanism:  - Acoustic phonons
C                            - Intervalley g-phonons
C                            - Intervalley f-phonons    
C                            - Coulomb scattering
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    CONVERT THIS CODE FROM ELECTRONS IN SI TO HOLES IN SI
C     ToDo Changes:          - 3 Valleys are HH,LH,SO bands
C                              - HH=1, LH=2, SO=3     
C                            - Degeneracies == 1 (change final valleys to 1)
C                            - Change Phonon energies
C                            - Add Alloy Disorder scattering    
C                            - Add Coulomb scattering
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

C    Acoustic phonons scattering rate

      if(acoustic_scattering.eq.1)then
         coupling_constant = sigma_acoustic
         i_mech = 1
         
         init_valley = 1
         ifin_valley = 1
         final_valleys = 1.D0
         call acoustic_scatt(n_lev,init_valley,ifin_valley)

         init_valley = 2
         ifin_valley = 2
         final_valleys = 1.D0
         call acoustic_scatt(n_lev,init_valley,ifin_valley)
         
         init_valley = 3
         ifin_valley = 3
         final_valleys = 1.D0
         call acoustic_scatt(n_lev,init_valley,ifin_valley)
                  
      endif
      
C    Surface-roughness scattering rate

      if(surface_roughness.eq.1)then
      
         i_mech = 3
         
         init_valley = 1
         ifin_valley = 1
         final_valleys = 1.D0
         call roughness(n_lev,init_valley,ifin_valley)

         init_valley = 2
         ifin_valley = 2
         final_valleys = 1.D0
         call roughness(n_lev,init_valley,ifin_valley)
         
         init_valley = 3
         ifin_valley = 3
         final_valleys = 1.D0
         call roughness(n_lev,init_valley,ifin_valley)
                  
      endif             
      
C    Alloy disorder scattering rate

      if(alloy_disorder.eq.1)then
         coupling_constant = pot_alloy
         i_mech = 1
         
         init_valley = 1
         ifin_valley = 1
         final_valleys = 1.D0
         call alloy_scatt(n_lev,init_valley,ifin_valley)

         init_valley = 2
         ifin_valley = 2
         final_valleys = 1.D0
         call alloy_scatt(n_lev,init_valley,ifin_valley)
         
         init_valley = 3
         ifin_valley = 3
         final_valleys = 1.D0
         call alloy_scatt(n_lev,init_valley,ifin_valley)
                  
      endif
      
C    Intervalley scattering: zero-order interaction (g-process)
C
C      if(intervalley_zero_g.eq.1)then
C         w0 = phonon_zero_g
C         coupling_constant = DefPot_zero_g
C         final_valleys = 1.D0
C         i_mech = 1
C        
C         init_valley = 1
C         ifin_valley = 1    
C         call intervalley(n_lev,w0,
C     1                    init_valley, ifin_valley)
C         init_valley = 2
C         ifin_valley = 2
C         call intervalley(n_lev,w0,
C     1                    init_valley, ifin_valley)
C     
C         init_valley = 3
C         ifin_valley = 3 
C         call intervalley(n_lev,w0,
C     1                    init_valley, ifin_valley)        
C         
C      endif

C    Intervalley scattering: zero-order interaction (f-process)

      if(intervalley_zero_f.eq.1)then
         w0 = phonon_zero_f
         coupling_constant = DefPot_zero_f
C         final_valleys = 2.D0
         final_valleys = 1.D0
         i_mech = 2
C HH to LH band
         init_valley = 1
         ifin_valley = 2    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley)
C HH to SO band
         init_valley = 1
         ifin_valley = 3    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley)   
C LH to HH band
         init_valley = 2
         ifin_valley = 1    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley)
C LH to SO band 
         init_valley = 2
         ifin_valley = 3    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley) 
C SO to HH band     
         init_valley = 3
         ifin_valley = 1    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley)
C SO to LH band     
         init_valley = 3
         ifin_valley = 2    
         call intervalley(n_lev,w0,
     1                    init_valley, ifin_valley)   
      endif
        
      return
      end
