       subroutine screen(q_vec,init_valley1,isub1,jsub1,
     1                   screened_element)
       implicit real*8 (a-h, o-z)
       
       real*8 kb
       real*8 wavefunc_1(10000,30)
       real*8 rn_sub_1(30)
       real*8 overlap_factor_1(30)
       real*8 overlap_fac(3,10,10)
       real*8 screen_vector(30)
       real*8 r_md_1(30)
       real*8 form_factor(30,30)
       real*8 matrix_element(30,30)
       real*8 eps_inverse(30,30)
       real*8 screened_factor(30)
       integer indx(30)
       
        
       common
     &/fund_const/a0,q,h,hbar,kb,am0,eps_0
     &/temp/tem,Vt
     &/mas_dos/r_md(3)
     &/dielec_func/eps_high
     &/dielec_func_oxide/eps_ox                                    
     &/index/k_sub(3)        
     &/wavefunc/wavefunc(10000,10,3),sub_energy(10,3)
     &/x_vector/x_vec(10000),j_max         
     &/sub_pop/rn_sub(10,3)
     &/sr_flag/sr_flag
     &/overlap_sr/overlap_sr(3,10,10)
        
C      Mapping the appropriate matrix element that needs to be screened

       if(sr_flag.eq.1)then
          do ival = 1,3
          do isub = 1, k_sub(ival)
          do jsub = 1, k_sub(ival)
             overlap_fac(ival,isub,jsub) = overlap_sr(ival,isub,jsub)   
          enddo
          enddo
          enddo
       endif

C      Creating a 1D array for the calculation of the screened diagonal matrix
C      elements
              
       kk = 1
       do ival = 1,3
       do isub = 1,k_sub(ival)
          do j = 1, j_max
             wavefunc_1(j,kk)=wavefunc(j,isub,ival)
          enddo
          overlap_factor_1(kk) = overlap_fac(ival,isub,isub)
          rn_sub_1(kk) = rn_sub(isub,ival)
          r_md_1(kk) = r_md(ival)
          if(init_valley1.eq.ival.and.isub1.eq.isub.
     1       and.jsub1.eq.isub)then
             kk_fix = kk
          endif
          kk = kk + 1      
       enddo
       enddo
       kk_tot = k_sub(1) + k_sub(2) + k_sub(3)
       
C      Calculation of the screening wavevector

       do kk = 1, kk_tot
       
          factor_1 = q*q/2.D0/eps_high/q_vec*rn_sub_1(kk)/kb/tem
          factor_2 = dsqrt(2.D0*r_md_1(kk)*kb/h*tem/h)
          argument = 0.5D0*q_vec*dsqrt(h/2.D0/r_md_1(kk)*h/kb/tem)
          dy = argument/20.D0
          sum = 0.D0
          do j = 1,20
             y = dfloat(j)*dy
             sum = sum + dexp(y*y)*dy
          enddo
          factor_3 = 2.D0*dexp(-argument**2)*sum
          screen_vector(kk) = factor_1*factor_2*factor_3
          
       enddo
              
C      Calculation of the form_factors

       dielec_factor = (eps_high - eps_ox)/(eps_high + eps_ox)
       do ii = 1, kk_tot
       do jj = 1, kk_tot
       
          sum = 0.D0
          zz = 0.D0
          
          do kz = 1,j_max
             zp = 0.D0
             do kzp = 1, j_max
             
                term1 = dexp(-q_vec*dabs(zz-zp))
                term2 = dielec_factor*dexp(-q_vec*(zz+zp))
                term_tot = term1 + term2
                wavefunc1 = wavefunc_1(kz,ii)**2
                wavefunc2 = wavefunc_1(kzp,jj)**2
                sum = sum + 
     1                wavefunc1*term_tot*wavefunc2*x_vec(kz)*x_vec(kzp)
                zp = zp + x_vec(kzp)
                
             enddo
             zz = zz + x_vec(kz)
          enddo
          
          form_factor(ii,jj) = sum
                 
       enddo
       enddo          
       
C      Evaluate the elements of the screening matrix

       do ii = 1, kk_tot
       do jj = 1, kk_tot
       
          if(ii.eq.jj)then
             matrix_element(ii,jj) = 1.D0 + form_factor(ii,jj)*
     1                 screen_vector(jj)/q_vec
c             matrix_element(ii,jj) = 1.D0          
          else         
             matrix_element(ii,jj) = form_factor(ii,jj)*
     1                 screen_vector(jj)/q_vec
c             matrix_element(ii,jj) = 0.D0                   
          endif
       
       enddo
       enddo

C      Evaluate the inverse of the dielectric function
       
       do ii = 1, kk_tot
          do jj = 1, kk_tot
             eps_inverse(ii,jj) = 0.D0
          enddo
          eps_inverse(ii,ii) = 1.D0
       enddo
       
       call LUDCMP(matrix_element,kk_tot,kk_tot,INDX,D)
       do jj = 1, kk_tot          
          call LUBKSB(matrix_element,kk_tot,kk_tot,INDX,
     1                eps_inverse(1,jj))
       enddo
       
c       do ii = 1, kk_tot
c       do jj = 1, kk_tot
c          print*,ii,jj,eps_inverse(ii,jj)
c       enddo
c       enddo
       
C      Calculate the screened matrix elements

       do ii = 1,kk_tot
          sum = 0.D0
          do jj = 1, kk_tot
             sum = sum + eps_inverse(ii,jj)*overlap_factor_1(jj)
          enddo
          screened_factor(ii) = sum      
       enddo
       
       screened_element = screened_factor(kk_fix)
       
       print*,'screening subroutine'
       print*,'Bare element = ',q_vec,overlap_factor_1(kk_fix)
       print*,'Screened element = ',screened_element                
        
       return
       end
       

C     Subroutines used for the matrix inversion

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit real*8 (a-h,o-z)
      
      real*8 A(30,30),B(30)
      integer INDX(30)
      
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END



      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      implicit real*8 (a-h,o-z)
      
      PARAMETER (NMAX=100,TINY=1.0D-20)
      real*8 A(30,30),VV(NMAX)
      integer INDX(30)
      
      D=1.D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END


       
