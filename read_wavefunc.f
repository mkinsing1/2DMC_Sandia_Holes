       subroutine read_wavefunctions() 
       implicit real*8(a-h,o-z)

       real*8 Ec(10000)
       real*8 y_vec(10000)
       real*8 Na(10000)
       integer nsub_1, nsub_2
       real*8 sum, term1,term2
       real*8 kb
       
       common 
     &/wavefunc/wavefunc(10000,10,3),sub_energy(10,3)
     &/sub_pop/rn_sub(10,3)
     &/x_vector/x_vec(10000),j_max
     &/overlap_intervalley/overlap_intervalley(3,10,3,10)
     &/overlap_sr/overlap_sr(3,10,10)
     &/overlap_srq/overlap_srq(3,10,10,1000)
     &/screening_terms/q_vector(1000),q_max,dq_max,nq_max     
     &/screened_sr_element/screened_sr(3,10,10,1000)     
     &/index/k_sub(3)
     &/fund_const/a0,q,h,hbar,kb,am0,eps_0
     &/pi/pi,two_pi
     &/dielec_func/eps_high
     &/dielec_func_oxide/eps_ox
     &/sr_flag/sr_flag

          nsub_1 = 2
          nsub_2 = 1
          
          k_sub(1) = nsub_1
          k_sub(2) = nsub_2
          k_sub(3) = nsub_2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C         Read the wavefunctions from the file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          j = 0
          open(unit = 1, file ='wavefunctions.dat',status='old')
          x0 = 0.D0
66        read(1,*,END=88)x,Ec(j),
     1                    (wavefunc(j,i,1),i=1,nsub_1),
     2                    (wavefunc(j,i,2),i=1,nsub_2)         
          x = x*1.D-9
          if(j.ne.0)then
             x_vec(j-1) = x - x0
          endif
          j = j + 1
          x0 = x
          goto 66
88        j_fix = j-1
          j_max = j_fix
          close(1)

          do j = 1, j_fix
          do i = 1, nsub_2
             wavefunc(j,i,3)=wavefunc(j,i,2)
          enddo
          enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C         Read the subband energies from a file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          open(unit=3, file ='sub_energy.dat', status = 'old')
          read(3,*)Vg, RNs, Ef,
     1            (sub_energy(i,1),i=1,nsub_1),
     2            (sub_energy(i,2),i=1,nsub_2)
          close(3)
          
          do i = 1, nsub_2
             sub_energy(i,3) = sub_energy(i,2)
          enddo
          do ival = 1,3
          do isub = 1,k_sub(ival)
             sub_energy(isub,ival)=1.D-3*sub_energy(isub,ival)
          enddo
          enddo  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C         Read the subband population from a file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          open(unit=4, file ='sub_pop.dat', status = 'old')
          read(4,*)Vg,RNs,delta_2,delta_4,
     1            (rn_sub(i,1),i=1,nsub_1),
     2            (rn_sub(i,2),i=1,nsub_2)
          close(4)
          
          do i = 1, nsub_1
             rn_sub(i,1) = RNs*rn_sub(i,1)*1.D2
          enddo
          do i = 1, nsub_2
             rn_sub(i,2) = RNs*rn_sub(i,2)*1.D2/2.D0
             rn_sub(i,3) = rn_sub(i,2)
          enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C         Calculate the depletion charge density and the effective field
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          j = 0
          open(unit = 1, file ='cond_band.dat',status='old')
          y0 = 0.D0
77        read(1,*,END=99)y,potential,e_density,charge_density,
     1                    doping_density
          Na(j) = charge_density
          Na(j) = Na(j)*1.D6        
          y = y*1.D-9
          if(j.ne.0)then
             y_vec(j-1) = y - y0
          endif
          j = j + 1
          y0 = y
          goto 77
99        j1_fix = j-1
          close(1)

          RN_depl = 0.D0
          do j = 0, j1_fix - 1
             Rn_depl = Rn_depl + Na(j)*y_vec(j)
          enddo

          Rn_depl = dabs(Rn_depl)
          RNs = RNs*1.D4
          E_eff = q/eps_high*(-0.5D0*RNs + Rn_depl)
          print*,'E_eff [V/cm] = ',E_eff/100.D0
          print*,'Sheet electron density = ',RNs
          print*,'Depletion charge density = ',Rn_depl
          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         Calculate the overlap  for intervalley scattering
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          do in_val = 1,3
          do ifin_val = 1,3
          
          do in_sub = 1,k_sub(in_val)
          do ifin_sub = 1,k_sub(ifin_val)         

                sum = 0
                do j = 0, j_fix
                   term1 = wavefunc(j,in_sub,in_val)**2
                   term2 = wavefunc(j,ifin_sub,ifin_val)**2
                   sum = sum + term1*term2 *x_vec(j)
                enddo  
                overlap_intervalley(in_val,in_sub,ifin_val,ifin_sub)=sum  
                
          enddo
          enddo
          enddo
          enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         Calculate the overlap for surface roughness scattering
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C         Calculate the bare matrix elements

          do in_val = 1,3        
          do in_sub = 1,k_sub(in_val)
          do ifin_sub = 1,k_sub(in_val)         

                sum = 0
                do j = 1, j_max-1
                   e_field = (Ec(j+1)-Ec(j))/x_vec(j)
                   term1 = wavefunc(j,in_sub,in_val)*e_field*
     1                     wavefunc(j,ifin_sub,in_val)
                   term2 = sub_energy(ifin_sub,in_val)*
     1                     (wavefunc(j+1,in_sub,in_val)-
     2                      wavefunc(j,in_sub,in_val))*
     3                      wavefunc(j,ifin_sub,in_val)/x_vec(j)
                   term3 = sub_energy(in_sub,in_val)*
     1                     (wavefunc(j+1,ifin_sub,in_val)-
     2                      wavefunc(j,ifin_sub,in_val))*
     3                      wavefunc(j,in_sub,in_val)/x_vec(j)
                   term = term1 - term2 + term3
                   term = term1
                   sum = sum + term*x_vec(j)
                enddo     
                overlap_sr(in_val,in_sub,ifin_sub)=sum
            
          enddo
          enddo
          enddo
         
          return                     
          end
