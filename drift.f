CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PERFORM THE K-SPACE AND REAL-SPACE MOTION OF THE CARRIERS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine drift(tau)
       implicit real*8(a-h, o-z)
 
       real*8 kx,ky,e,tau
       real*8 qh1,dkx,dky
       real*8 fx,fy
       integer iv
       
       common
     &/dri/qh
     &/mass/smh(3),hhm(3)
     &/force/fx,fy
     &/particle_atr/kx,ky,iv,is,e
                          
       qh1 = -qh*tau
       dkx=qh1*fx
       dky=qh1*fy
       
       kx = kx+dkx
       ky = ky+dky

       skx = kx*kx
       sky = ky*ky
       sk = skx+sky
       gk = hhm(iv)*sk
       e = gk

       return
       end
               		
