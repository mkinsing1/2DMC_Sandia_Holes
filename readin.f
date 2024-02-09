CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      READ INPUT PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine readin(n_lev)
       implicit real*8 (a-h, o-z)
       
       real*8 kb, am0
       real*8 fx,fy

       common
     &/ran_var/iso
     &/pi/pi,two_pi
     &/fund_const/q,h,kb,am0,eps_0 
     &/dri/qh
     &/temp/tem,Vt
     &/mass/smh(3),hhm(3)
     &/scatt_par/emax,de,w(3,10,50),tau_max(3,10)  
     &/density/density,sound_velocity
     &/dielec_func/eps_high
     &/dielec_func_oxide/eps_ox
     &/time_1/dt,dtau,tot_time
     &/force/fx,fy
     &/select_acouctic/acoustic_scattering
C     &/select_intervalley_1/intervalley_zero_g
     &/select_intervalley_2/intervalley_zero_f
     &/select_sr/surface_roughness
     &/select_alloy/alloy_disorder
     &/select_coulomb/coulomb_scattering
     &/sigma_acoustic/sigma_acoustic
     &/Def_pot_1/DefPot_zero_g
     &/Def_pot_2/DefPot_zero_f
     &/interval_phonons_1/phonon_zero_g
     &/interval_phonons_2/phonon_zero_f
     &/surface_roughness/delta,corr_length
     &/mas_dos/r_md(3)
     &/mas_con/r_mc(3)

C      Define fundamental constants and general parameters

       iso=1345
       am0=9.11D-31
       h=1.05459D-34
       q=1.60219D-19
       qh=q/h  
       eps_0=8.85419D-12
       kb=1.38066D-23
       pi=3.161592654D0
       two_pi = 2.D0*pi

C      Define time step and maximum simulation time

       dt=1D-13
       tot_time = 20.D-12

C      Set temperature and doping density in various regions

       tem=300.D0
       Vt=kb*tem/q

C      Set the electric field

       fx = 10.D5       
       fy = 0.D0
               
C      Define high-frequency dielectric constant

       eps_high = 10.92D0
       eps_high = eps_high*eps_0
       eps_ox = 3.9D0
       eps_ox = eps_ox*eps_0

C      Define crystal density and sound velocity

       density = 2329.D0
       sound_velocity = 9040.D0

C      Define parameters for the scattering table

       emax=1.D0
       de=emax/dfloat(n_lev)

C      Select scattering mechanisms

       acoustic_scattering = 1
C       intervalley_zero_g = 1  (Not needed becuase no degeneracies for HH,LL,SO hole bands)
       intervalley_zero_f = 1
       surface_roughness = 1
C      alloy_disorder = 1  (Need to add)
C      coulomb_scattering = 1 (Need to add)

C      Define coupling constants

       sigma_acoustic = 9.D0  !  [eV]
       DefPot_zero_g = 5.23D10   ! [eV/m]
       DefPot_zero_f = 5.23D10   ! [eV/m]
   
       phonon_zero_g = 0.06D0  !  [eV]
       phonon_zero_f = 0.06D0  !  [eV]

C      Define Longitudinal and transverse masses
 
       am_l = 0.91D0
       am_t = 0.19D0
       r_md(1) = am_t*am0
       r_md(2) = sqrt(am_l*am_t)*am0
       r_md(3) = sqrt(am_l*am_t)*am0

       r_mc(1) = 2.D0*am0/(1.D0/am_t+1.D0/am_t)
       r_mc(2) = 2.D0*am0/(1.D0/am_l+1.D0/am_t)
       r_mc(3) = 2.D0*am0/(1.D0/am_l+1.D0/am_t)
       
       do i = 1,3
          smh(i) = dsqrt(2.D0*r_mc(i))*dsqrt(q)/h
          hhm(i) = h/r_mc(i)/q*h/2.D0
       enddo
       
C      define parameters for surface roughness
       delta = 2.84D-10
       corr_length = 24.2D-10       

       return
       end  
       
