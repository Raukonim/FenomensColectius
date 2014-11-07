c23456789012345678901234567890123456789012345678901234567890123456789012
c     ******************************************************************
c     *        SIMULACIO D'UNA CADENA DE MARKOV DISCRETA               *
c     *                                                                *
c     *        GENERADOR RCARRY                                        *
c     *                                                                *
c     *                                                                *
c     *        N     - nombre d'estats  (<20)                         *
c     *        P0(1:N)  - probabilitat inicial                         *
c     *        W(1:N, 1:N)  - matriu de probabilitats de transicio     *
c     *        NLLAV - nombre de llavors                               *
c     *        NPAS  - nombre de passes de la cadena (<1000)           *
c     *                                                                *
c     ******************************************************************
c
c     * declaracio de variables ****************************************
c
      implicit none
      integer idau
      integer N, NPAS, I, J, K, ISTAT
      REAL*4 P0(1:20), W(1:20,1:20), PROB(1:20)
      real*4 suma
      real*4 xhist(1:20, 0:1000)
c
c     * aquest es el vector que contindra els nombres a l'atzar
c     * la dimensio de rvec ha de ser NPAS+24 com a minim
c
      real*4 rvec(1024),rand
      integer ivec
      integer llav, nllav, llav0
c
c     dades inicials
c

      N = 6
      llav0 = 1347
      NPAS=50
      nllav = 1000000
c
c     llegim P0 i comprovem la suma
c
      suma=0.0d0
      open (unit=9,file="Borratxo-p0.dat")
     
      write(*,*) 'Llegim P0s'

      do i=1,N
         read(9,*) P0(i)
         suma=suma+P0(i)
      enddo
      
      write (*,*) 'Suma probabilitats inicials = ', suma

      close (9)
c
c     Llegim W(i->j) i comprovem la suma
c
      open(unit=10, file="Borratxo-W.dat")

      write(*,*) 'Llegim Ws'

      do i=1,N
         suma=0.0d0
         do j=1,N
            read (10,*) W(i,j)
            suma=suma+W(i,j)
         enddo
         write (*,*) 'Suma prob de transicio des de ', i, ' = ', suma
      enddo

      close(10)


c
c     Inicialitzem el generador de nombres a l'atzar
c
      call rcarin(llav,rvec,NPAS+1)
c
c     Bucle de llavors
c
      do llav=llav0, llav0+nllav-1, 1

c
c        Generem NPAS+1 nombres a l'atzar
c

         call rcarry(rvec,NPAS+1)
c
c        Amb el primer nombre generem l estat incial
c
         ivec=1
         rand=rvec(ivec)
         ivec=ivec+1
         ISTAT=IDAU(rand,N,P0)

         xhist(istat,0) = xhist(istat,0)+1.0
c
c        Passes de la cadena
c
         do k=1,NPAS
            
            do j=1,N
               PROB(j)=W(ISTAT,J)               
            enddo

            rand=rvec(ivec)
            ivec=ivec+1

            ISTAT = IDAU(rand,N,PROB)

            xhist(istat,k)=xhist(istat,k)+1.0

         enddo      

      enddo

c
c     Normalitzacio de xhist
c
      do i=1,N
        do k=0,NPAS
           xhist(i,k) = xhist(i,k)/REAL(NLLAV)
        enddo
      enddo
c
c
c
      open (unit=16,file='Markov.out')
   
      do k=0,NPAS
        write(16,1699) k,(xhist(i,k), i=1,N)
 1699   format (i4, 20(1X, f7.5))
      enddo
    
      close (16)

      stop
      end

c     ******************************************************************
c     *                    SUBRUTINA IDAU                               *
c     *****************************************************************
      INTEGER FUNCTION IDAU(RAND,N,PROB)
      real*4 RAND,XLIM
      integer N,K,I0
      real*4 PROB(1:20)
      XLIM=0.0D0
      do k=1,N
         XLIM=XLIM+PROB(k)
         IF (RAND.LT.XLIM) then
            IDAU=K
            RETURN
         ENDIF
      enddo
      END
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
     
