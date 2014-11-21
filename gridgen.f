c23456789012345678901234567890123456789012345678901234567890123456789012
c     ******************************************************************
c     *       Grid generator for an isin simulation                    *
c     *                                                                *
c     *                                                                *
c     ******************************************************************
c
c
c     Spin states +-1 in a 2D grid-> 2 coordinates: i, j from 1 to 256?
c     matrix [[rvec][rvec]....[rvec]] if rvec(i,j).ge.0.5 sij=1 else -1
c
c     * variable declarations ****************************************
c
      Implicit None
      Real*8 Magne, Ener                 ! Declarem les funcions
      Real*8 magin, enerin
      Integer*4 llav, i ,j, N, L, sij, k
      Integer*4 Nrand, ivec, PBC(0:10+1), montec
      
      
      Integer*2 S(1:256,1:256)
c     * Dimensió rvec = 256*256 + 24
      Real*4 rvec(65536+24)

      L=10
      N=L*L
      
      PBC(0)=L
      PBC(L+1)=1

      Do k=1,L
        PBC(k)=k
      EndDo
      
c     * generator seed

      llav=17258
      
c     * series elements

      Nrand=N
c     * random numbers generation
      
      Call rcarin(llav,rvec,Nrand)
      Call rcarry(rvec,Nrand)
      
c     * Open file

      Open(13,File='grid.out')
      
      ivec=1
      
c     * do loops for matrix inicialization     
      Do i=1,L
        Do j=1,L
          If ((rvec(ivec)).ge.0.5)Then
            S(i,j) = 1
          Else
            S(i,j) = -1
          EndIf
          ivec=ivec+1
          Write(*,*) i,' ',j,' ', S(i,j)
        EndDo
        Write(13,*)( S(i,j) ,j=1,N)
      EndDo
c     cridem la funció Magne que calcula la magnetització de la xarxa
      magin=Magne(S,L)
      Write(*,*) "La magnetitzacció inicial és: ", magin
c     cridem la funcií Ener que calcula l'energia de la xar
      enerin=Ener(S,L, PBC)
      Write(*,*) "L'Enrgia inicial és:  ", enerin
      Stop
      End
      
      
c      
c     Funció Magne que calcula la magnetització per spin 

      Real*8 Function Magne(S,L)
      
      Integer*2 S(1:256,1:256)
      Integer*4 N, suma, L, i, j
      
      suma=0
      N=L*L
      
      Do i=1,L
        Do j=1,L
          suma=suma+S(i,j)
        Enddo
      Enddo
      
      Magne= Real(suma)/Real(N)

c     Do de Montecarlo
      Do i=1,montec
c     Do N cops
         Do j=1,N
            ! codi metropolis
         EndDo
      EndDo

      Return
      End

c      
c     Funció Ener que calcula l'energia per spin

      Real*8 Function Ener(S,L,PBC)
      
      Integer*2 S(1:256,1:256)
      Integer*4 N, suma, L, k, PBC(0:10+1)
      
      suma=0
      
      Do i=1,L
        Do j=1,L
          suma=suma-S(i,j)*(S(PBC(i+1),j)+S(i,PBC(j+1)))
        EndDo
      EndDo
      
      Ener=Real(suma)/Real(L*L)
      
      Return
      End


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
     
