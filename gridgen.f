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
      !Implicit None
      Real*8 Magne, Ener, Boltz               ! Declarem les funcions
      Real*8 magin, enerin, prob
      Real*8 eners, magnes, enerpas, magnepas
      Real*8 sumam, sume, sumam2, sume2, sumabsm, suma0
      Real*8 start, finish, vprob(1:8)
      Real*4 A, B, C, T
      Integer*4 illav, illav0, nllav, i ,j, N, L, sij, k
      Integer*4 IMC, MCtotal, MCini, MCstep, DeltaEner, enerdif
      Integer*4 Nrand, ivec, ipas, PBC(0:257), montec
      Integer*4 num_args, ix
      Character(Len=12) args(4)
      
      
      Integer*2 S(1:256,1:256)
c     * Dimensió rvec = 3*(256*256) + 24
      Real*4 rvec(196632)


      Call cpu_time(start)
 
c     Variables que poden ser modificades per obtenir corbes més smooth
c     les llegim com arguments de la linia de comandes
      !REAL(8)    :: A,B


      num_args = command_argument_count()
      !Allocate(args(num_args))  ! I've omitted checking the return status of the allocation 

      Do ix = 1, 4
        Call get_command_argument(ix,args(ix))
      EndDo
      read (args(1),*) L !10
      read (args(2),*) MCtotal !10000
      read (args(3),*) T !1.0d0
      read (args(4),*) nllav !30
      
      MCini=1000
      MCstep=10
      
c     * generator seed

c      illav=17258
      illav0=114

      
      PBC(0)=L
      PBC(L+1)=1

      Do k=1,L
        PBC(k)=k
      EndDo
      

c     * series elements
      N=L*L
      Nrand=3*N
      
      
c     Probs

      vprob(1)=0d0
      vprob(2)=0d0
      vprob(3)=0d0
      vprob(4)=dexp(-4.0d0/T)
      vprob(5)=0d0
      vprob(6)=0d0
      vprob(7)=0d0
      vprob(8)=dexp(-8.0d0/T)

      
      eners=0.0d0
      enerpas=0.0d0
      magnes=0.0d0
      magnepas=0.0d0
      sumam=0.0d0
      sume=0.0d0
      sumam2=0.0d0
      sume2=0.0d0
      sumabsm=0.0d0
      suma0=0.0d0
      
                

c     * Open file

c      Open(13,File='grid.out')
      Open(14,File='dades.out',position="append")
        
c     bucle de temperatures
c      Do T=0,5.0d0,0.1d0
c     bucle de llavors illav a nllav
      Do illav=illav0,illav0+nllav,1
c     * random numbers generation
      
        Call rcarin(illav,rvec,Nrand)

c     generació de la matriu inicial
        Call rcarry(rvec,Nrand)
      
      
      ivec=0
      
c     * do loops for matrix generation    
        Do i=1,L
          Do j=1,L
            If ((rvec(ivec)).ge.0.5)Then
              S(i,j) = 1
            Else
              S(i,j) = -1
            EndIf
            ivec=ivec+1
c          Write(*,*) i,' ',j,' ', S(i,j)
          EndDo
c          Write(13,*)( S(i,j) ,j=1,L)
        EndDo
c     cridem la funció Magne que calcula la magnetització de la xarxa
c        magin=Magne(S,L)
c        Write(*,*) "La magnetitzacció inicial és: ", magin
c     cridem la funció Ener que calcula l'energia de la xarxa
c        enerin=Ener(S,L, PBC)
c        Write(*,*) "L'Enrgia inicial és:  ", enerin


c     Do de Montecarlo
        Do IMC=1,MCtotal
          Call rcarry(rvec, Nrand)
          ivec=0
c        Do N cops
          Do ipas=1,N
            ! codi metropolis 3 nums aleatoris
            i=INT(L*rvec(ivec))+1
            ivec=ivec+1
            j=INT(L*rvec(ivec))+1
            ivec=ivec+1
c           Imprimim ij per comprovar que els index siguin del tipus
c           que desitgem
c            write(*,*) i, j

c           Calculem delta H
            enerdif=DeltaEner(i, j, S, PBC)
c           if per comprovar si l'increment energètic és negatiu
c           o positiu
            If (enerdif.le.0.0d0) Then
c              write(6,*)'acceptem el canvi'
               S(i,j)=-S(i,j)
c               calcul energia
c               eners=Ener(S, L, PBC)
c               magnes=Magne(S, L)
            Else
c             solament acceptem el canvi si l'exponencial de -deltaH/T 
c             és major que un nombre a l'atzar
c              write(*,*) enerdif
              C=rvec(ivec)
              ivec=ivec+1
c              prob=Boltz(enerdif, T)
c               Write(*,*) prob
              If (C.le.vprob(enerdif)) Then
c                             write(6,*)'acceptem el canvi'
                S(i,j)=-S(i,j)
c                calcul energia
c                eners=Ener(S, L, PBC)
c                magnes=Magne(S, L)
              EndIf
            EndIf
          EndDo
c        calcul dels promitjos(real*8) energia  imantacio sumam sume sumam2
c        sume2 sumabsm per mc prou grans mcini=1000 i multiples de mcstep=10
c        IF(mc.gt.mcini)AND(imc.eq.mcstep*(imc/imcstep)) 
c        definim fora bucle montecarlo sum* on hi afegim Sum0 que
c        compta el nombre de condicions complerrtes  de mcstep
c        Write ene1 i ene2 ca
          IF((IMC.gt.MCini).AND.((Mod(IMC,MCstep)).eq.0)) Then
            enerpas = Ener(S, L, PBC)
            magnepas = Magne(S,L)
c            Write(*,*) IMC, eners, enerpas
c            Write(*,*) IMC, magnes, magnepas
            sumam = sumam + magnepas
            sume = sume + enerpas
            sumam2 = sumam2 + (magnepas**2)
            sume2 = sume2 + (enerpas**2)
            sumabsm = sumabsm + Abs(magnepas)
            suma0 = suma0+1.0d0
c            Write(14,*) suma0, enerpas, sume, sume2,
c     +      magnepas, sumam, sumam2, sumabsm
          EndIf
        EndDo
      EndDo
c     normalitzem divindint els sumatoris per suma0 que compta
c     la quantitat de passos que hem sumat
      
      sume=sume/suma0
      sume2=sume2/suma0
      sumam=+sumam/suma0
      sumam2=sumam2/suma0
      sumabsm=sumabsm/suma0
      
      Write(14,*) L, T, nllav, MCtotal, sume, sume2,
     +       sumam, sumam2, sumabsm, suma0
      
c      Close(13)
      Close(14)
      
      Call cpu_time(finish)
      print '("Time = ",f12.3," seconds.")',finish-start

      Stop
      End
      
c     ******************************************************************
c     *                       Funció Magne                             *
c     *                                                                *
c     *               Calcula la magnetització per spin                *
c     *                                                                *
c     ******************************************************************
c

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


      Return
      End

c     ******************************************************************
c     *                        Funció Ener                             *
c     *                                                                *
c     *                  Calcula l'energia per spin                    *
c     *                                                                *
c     ******************************************************************
c

      Real*8 Function Ener(S,L,PBC)
      
      Integer*2 S(1:256,1:256)
      Integer*4 N, suma, L, k, PBC(0:257)
      
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
c     *                        Funció DeltaEner                        *
c     *                                                                *
c     *                  calcula l'increment de l'energia              * 
c     *                            i, j, s(ij)                         *
c     *                                                                *
c     ******************************************************************
c


      Integer*4 Function DeltaEner(i, j, S,PBC)
      
      Integer*2 S(1:256,1:256)
      Integer*4 sumveins, PBC(0:257), i, j
      
      sumveins=0
      
      sumveins=(S(i,PBC(j+1))+S(i,PBC(j-1)))
      sumveins=sumveins+(S(PBC(i+1),j)+S(PBC(i-1),j))
      
      DeltaEner=2*S(i,j)*sumveins
      
c      Write(*,*) S(i,j), S(i,PBC(j+1)), S(i,PBC(j-1))
c     +    ,S(PBC(i+1),j), S(PBC(i-1),j), sumveins, DeltaEner
      
      Return
      End
      
c     ******************************************************************
c     *                        Funció Boltz                            *
c     *                                                                *
c     *                  Calcula l'exponencial de boltzmann            *
c     *                                                                *
c     ******************************************************************
c

      Real*8 Function Boltz(enerdif, T)
      
      Real*8 enerdif
      Real*4 T
      
  
      Boltz=dexp((-enerdif)/T)
      
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
     
