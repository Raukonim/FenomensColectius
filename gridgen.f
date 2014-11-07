c23456789012345678901234567890123456789012345678901234567890123456789012
c     ******************************************************************
c     *       Grid generator for an isin simulation                    *
c     *                                                                *
c     *                                                                *
c     ******************************************************************
c
c
c     Spin states +-1 in a 2D grid-> 2 coordinates: i, j from 0 to 1000?
c     matrix [[rvec][rvec]....[rvec]] if srvec(i,j).ge.0.5 sij=1 else -1
c
c     * variable declarations ****************************************
c
      implicit none
      integer llav, i ,j, N, sij
      
      
      integer grid(10024,10024)
      real*4 rvec(10024)
      
c     * generator seed

      llav=139177
      
c     * series elements
      
      N=1000
      
c     * random numbers generation
      
      call rcarin(llav,rvec,N)
      call rcarry(rvec,N)
      
c     * Open file

      open(13,file='grid.out')
      
c     * do loops for matrix inicialization     
      do i=1,N
            do j=1,N
                  if (rvec(j).ge.0.5) then
                        sij=1
                  else
                        sij=-1
                  endif
                  grid(i,j) = sij
            enddo
            write(13,*)( grid(i,j) ,j=1,N)
            call rcarry(rvec,N)
      enddo
      
      stop
      end
      
      
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
     
