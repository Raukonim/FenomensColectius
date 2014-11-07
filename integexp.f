c23456789012345678901234567890123456789012345678901234567890123456789012
c     ******************************************************************
c     *        INTEG MC                                                *
c     *                                                                *
c     *                                                                *
c     ******************************************************************
c
c     * declaracio de variables ****************************************
c
      implicit none
      integer llav
      integer j
c
c     * aquest es el vector que contindra els nombres a l'atzar
c     * la dimensio de rvec ha de ser N+24 com a minim
      real*4 rvec(50024)
c
      real*4 sx,sx2,ex
      integer N
c  
c     * llavor del generador
c
      llav=139177
c
c     * nombre d'elements de cada serie
c
      N=10000
c
c     Generem N uniformes 
c
      call rcarin(llav,rvec,N)
      call rcarry(rvec,N)
c
c     * variables pel calcul dels promitjos
c
      sx=0.0d0
      sx2=0.0d0

      open(13,file='integexp.out')
      
      do j=1,N
         sx=sx+exp(rvec(j))
         sx2=sx2+exp(rvec(j))*exp(rvec(j))
         ex=sx2/float(j)-(sx*sx)/(float(j)*float(j))
         ex=1.96*sqrt(ex)/sqrt(float(j))
         write(13,*) j, sx/float(j),ex
      enddo
c
c     * normalitzacio dels promitjos
c
      sx=sx/dfloat(N)
      sx2=sx2/dfloat(N)
c
c        sortida de resultats per pantalla
c
      write(*,*) ' llavor = ', llav
      write(*,*) ' <I>    = ',sx,' ',exp(1.0)-1.0
      write(*,*) ' V(I)   = ',sx2-sx*sx
      write(*,*) ' error  = ',1.96*sqrt((sx2-sx*sx)/N)
       
      close(13)  
      stop
      end

c
c     ******************************************************************
c     *                     SUBRUTINA RCARIN                           *
c     ******************************************************************
c

      SUBROUTINE RCARIN(IJKL,RVEC,LENV)

C----------------------------------------------------------------------
C Inicializa valores antes de llamar a la subrutina RCARRY.
C IJKL debe estar en el rango 0<IJKL<900 000 000.
C Para conseguir los valores standar usados por Marsaglia y Zaman en su
C articulo poner IJKL = 54217137 (I=12, J=34, K=56, L=78)
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
      COMMON /RAN1/ CARRY
      DIMENSION RVEC(LENV+24)

      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177,177) + 2
      J = MOD(IJ,177)     + 2
      K = MOD(KL/169,178) + 1
      L = MOD(KL,169)

      DO 2 II=24,1,-1
        S = 0.0
        T = 0.5
        DO 3 JJ=1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF (MOD(L*M,64).GE.32) S = S+T
          T = 0.5*T
3       CONTINUE
        RVEC(II) = S
2     CONTINUE

      CARRY = 0.0

      RETURN
      END
c
c     ******************************************************************
c     *                     SUBRUTINA RCARRY                           *
c     ******************************************************************
c
      SUBROUTINE RCARRY(RVEC,LENV)
C----------------------------------------------------------------------
C Generador de numeros pseudo-aleatorios. Algoritmo de G. Marsaglia y
C A. Zaman. Genera numeros reales de 32-bits con mantisas de 24 bits,
C comprendidos entre 0 y 1 (1, explicitamente excluido).
C Periodo aproximado : 10**171.
C Admite la generacion de subsecuencias disjuntas.
C                   F. James, 1989
C Version modificada (mas rapida que el original). (2/9/91)
C----------------------------------------------------------------------
      DIMENSION RVEC(LENV+24)
      COMMON /RAN1/ CARRY
      PARAMETER (TWOM24=1.0/16777216.0)
C
      DO 100 IVEC=25,LENV+24
        UNI = RVEC(IVEC-24) - RVEC(IVEC-10) - CARRY
        IF (UNI.LT.0.) THEN
          UNI = UNI + 1.0
          CARRY = TWOM24
        ELSE
          CARRY = 0.0
        ENDIF

        IF(UNI.EQ.0.)THEN
          UNI=RVEC(IVEC-24)*TWOM24
            in48=-48
          IF(UNI.EQ.0.)UNI=2**(in48)
        ENDIF

        RVEC(IVEC) = UNI
100   CONTINUE

      DO 200 I=1,24
200   RVEC(I)=RVEC(LENV+I)

      RETURN
      END
     
