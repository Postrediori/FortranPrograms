      SUBROUTINE PDDATA(ltwo)
C----------------------------------------------------------------------
C    BõáOB èPOñEÑìPõ èOÑÉOTOBKà àCXOÑHõX ÑAHHõX à èOÑÉOTOBKA HAÉPìáKà
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(1100),RP2(1100),RP(1100)
     */BL7/XY(8,2),AKT(16,16)
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PP(2,90),PT(2,90)
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1,ANU
C
      X(1)=-0.57735027
      X(3)=0.0
      X(2)=-X(1)
      C(1)=1.0
      C(2)=1.0
      C(3)=0.0
C
C !!!  XDATA - èPOñEÑìPA PAáÅàBKà OÅãACTà HA KOHEóHõE ùãEMEHTõ  !!!
C
      CALL XADATA
C
      write(6,200)
  200 FORMAT(//25X,'*** áAÑAóA èãOCKOâ ÑEîOPMAñàà ***')
C
      DO 300 I=1,NZ
 300  AMZ(I)=0.0
C
C      îOPMàPOBAHàE  HAÉPìáKà HA èEPBOM ìóACTKE
C
      NN=2*NNOD
C
      IF(NEP.LE.0)  GOTO 601
      DO 560 I=1,NN
  560 A1(I)=0.0
      CALL PDFOR(NN,1)
      DO 562 I=1,NN
  562 RP1(I)=A1(I)
C
C    îOPMàPOBAHàE  HAÉPìáKà  HA  BTOPOM  ìóACTKE
C
  601 IF(NEP1.LE.0)  GOTO 600
      DO 1 I=1,NEP1
      LTEL(I)=LTEL(NEP+I)
      LST(I)=LST(NEP+I)
      PP(1,I)=PP(1,NEP+I)
      PP(2,I)=PP(2,NEP+I)
      PT(1,I)=PT(1,NEP+I)
    1 PT(2,I)=PT(2,NEP+I)
      if(ltwo.eq.1) then
         do i=(nep1/2+1),nep1,2
            pp(1,i)=-pp(1,i)
            pp(2,i)=-pp(2,i)
         end do
      end if
      NEP=NEP1
      DO 582 I=1,NN
  582 A1(I)=0.0
      CALL PDFOR(NN,1)
      DO 564 I=1,NN
  564 RP2(I)=A1(I)
  600 CONTINUE
      RETURN
      END
C*****
      SUBROUTINE PDFOR(NN,IFF)
C----------------------------------------------------------------------
C
C               îOPMàPOBAHàE MATPàñõ A(NN*LH)
C               à èPABOâ óACTà A1(NN) CàCTEMõ
C
C----------------------------------------------------------------------
      COMMON
     */BL1/A(400000)
     */BL1A/LH,A1(1100)
     */BL2A/NEL,NC(8,110)
     */BL7/XY(8,2),AK(16,16)
     */BL15/LEC(16),JIN(8)
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PP(2,90),PT(2,90)
     */BLHT/HPT(16)
C
      JJ=NN*LH
      IF(IFF) 5,5,7
    5 DO 10 I=1,JJ
   10 A(I)=0.0
      GO TO 13
    7 DO 6 I=1,NN
    6 A1(I)=0.0
C
  13  IP=1
C
      DO 50 JJ=1,NEL
C
      DO 49 I1=1,16
   49 HPT(I1)=0.0
C
      CALL PDELM(JJ,IG)
      IF(IFF.EQ.0) CALL PDSTIF(JJ,IG)
C
      IF(IFF) 44,44,48
   48 IF(NEP) 44,44,41
C
C       PACèPEÑEãEHHAü HAÉPìáKA
C
   41 DO 43 J=1,5
      IF(LTEL(IP)-JJ) 44,42,44
   42 CALL PDPRES(LST(IP),PP(1,IP),PP(2,IP),PT(1,IP),PT(2,IP),HPT)
   43 IP=IP+1
   44 DO 50 J=1,16
      M=LEC(J)
      A1(M)=A1(M)+HPT(J)
C
      IF(IFF) 55,55,50
   55 DO 45 K=1,16
      IF(LEC(K).LT.LEC(J)) GOTO 45
      M=(LEC(J)-1)*LH+LEC(K)-LEC(J)+1
      A(M)=A(M)+AK(J,K)
   45 CONTINUE
   50 CONTINUE
C
   81 DO 90 I=1,NZ
      IF(IFF.GT.0)  GOTO 85
C
C        ÜECTKOE áAKPEèãEHàE
C
      IF(IFF.GT.0)  GOTO 85
      M=(MZD(I)-1)*LH+1
      J=M+LH-1
      DO 80 K=M,J
   80 A(K)=0.0
      A(M)=1.0
   85 M=MZD(I)
      A1(M)=0.0
   90 CONTINUE
      IF(IFF.GT.0) GOTO 110
      DO 100 I=1,NN
      DO 95 J=1,NZ
      IF(I.LE.MZD(J)-LH.OR.I.GE.MZD(J)) GOTO 95
      M=(I-1)*LH+MZD(J)-I+1
      A(M)=0.0
   95 CONTINUE
  100 CONTINUE
  110 RETURN
      END
C*****
      SUBROUTINE PDFF(N)
      COMMON
     */BL1A/LH,AA1(1100)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(1100),RP2(1100),RP(1100)
     */BL15/LEC(16),JIN(8)
     */BLS/SUM(10000),FP(10000)
     */BLTEOR/ITEOR
     */blita/m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
     */bltcom/tcom(1500000)
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
      DIMENSION DUL(16),DS(4),DR(3),de1(3),UL(16),DE(4)
C
      NN=NNOD*2
C
      DO 10 I=1,N
   10 FP(I)=0.0
C
      numiter=3
      LL=0
      L2=0
   20 LL=LL+1
C
      DO 30 I=1,NN
   30 AA1(I)=RP(I)
           CALL PDFOR(NN,0)
           CALL PDSTRT(1,NN,NN,-1)
           CALL PDBSUB(1,NN,0)
      DO 40 II=1,NEL
           CALL PDELM(II,IG)
      DO 50 I=1,16
      J=LEC(I)
      UL(I)=SUM(J)
   50 DUL(I)=AA1(J)
      KK=1
      DO 40 M=1,2
      DO 40 N2=1,2
           CALL PDCOB(X(M),X(N2),IG)
           CALL PDDEF(DUL,DE)
           CALL PDDDD(II,KK)
           CALL PDDES(DE,  DS,DR,de1,DK)
      N1=NN+(II-1)*32+(KK-1)*8
      DO 80 J=1,4
   80 FP(N1+J)=DS(J)
      DO 90 J=1,3
   90 FP(N1+4+J)=DR(J)
      FP(N1+8)=DK
      IF(LL.NE.numiter) GO TO 45
      sxx=sum(n1+1)
      syy=sum(n1+2)
      sxy=sum(n1+3)
      dsxx=ds(1)
      dsyy=ds(2)
      dsxy=ds(3)
      dexx=de(1)
      deyy=de(2)
      dexy=de(3)
      www(kk,ii)=www(kk,ii)+((sxx+dsxx)*dexx+
     +           2.0*(sxy+dsxy)*dexy+(syy+deyy)*deyy)/2.0
      alambd(kk,ii)=alambd(kk,ii)+sqrt(de1(1)*de1(1)+de1(2)*de1(2)+
     +           de1(3)*de1(3))
           CALL PDST3(DR,DK)
      N3=(II-1)*4+KK
      DO I=1,LP
           TCOM(LP*(N3-1)+I)=T(I)
      END DO
   45 KK=KK+1
      L2=L2+NTA
   40 CONTINUE
C
      IF(LL.EQ.1) GO TO 20
      IF(L2.GT.0.AND.LL.LT.numiter) GOTO 20
C
      DO 110 I=1,NN
  110 FP(I)=AA1(I)
C
      PRINT 120,L2
  120 FORMAT(20X,'óàCãO ÄäíàÇçõï çÄèêÄÇãÖçàâ PABHO  ',I6/)
C
      RETURN
      END
      SUBROUTINE PDSE(K,DA,LA,K1,A0)
C----------------------------------------------------------------------
C     A,DA(3)-BEKTOP à EÉO èPàPAôEHàE B TOóKE-K;
C     LA(K1)-HOMEPA ìÉãOBõX TOóEK TPAEKTOPàà;
C     K1-HOMEP áBEHA,OTBEóAûôEÉO TOóKE-K;
C     A0(3,K1)-KOOPÑàHATõ K1-Oâ ìÉãOBOâ TOóKà
C     (èPEÑBAPàTEãúHO HEOÅXOÑàMO èOãOÜàTú K1=1)
C----------------------------------------------------------------------
      DIMENSION DA(3),LA(10),A0(3,10)
C
      N1=LA(K1)
      N2=LA(K1+1)
      DO 10 I=1,3
   10 DA(I)=(A0(I,K1+1)-A0(I,K1))/(N2-N1)
      IF(K.EQ.N2-1) K1=K1+1
      RETURN
      END
C
      SUBROUTINE PDSOL(N)
C----------------------------------------------------------------------
C    PEòEHàÖ CàCTEMõ Ñàî.ìPABHEHàâ METOÑOM ùâãEPA
C----------------------------------------------------------------------
      COMMON
     */BLS/SUM(10000),FP(10000)
C
      CALL PDFF(N)
      DO 1 I=1,N
    1 SUM(I)=SUM(I)+FP(I)
      RETURN
      END
C
      FUNCTION PDRR(X,XI,ET)
      DIMENSION AN(8),X(8)
      X1=1.0-XI
      X2=1.0+XI
      E1=1.0-ET
      E2=1.0+ET
      AN(1)=-X1*E1*(X2+ET)/4.0
      AN(3)=X2*E1*(XI-E2)/4.0
      AN(5)=X2*E2*(ET-X1)/4.0
      AN(7)=-X1*E2*(E1+XI)/4.0
      AN(2)=X1*X2*E1/2.0
      AN(4)=E1*E2*X2/2.0
      AN(6)=X1*X2*E2/2.0
      AN(8)=E1*E2*X1/2.0
      PDRR=0.0
      DO 1 I=1,8
    1 PDRR=PDRR+X(I)*AN(I)
      RETURN
      END
C
      SUBROUTINE XBLOCK(NELX,ITEORX)
C----------------------------------------------------------------------
C
C     èEPEÑAóA ÑAHHõX óEPEá ÅãOKà COMMON
C
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2/XYT(2,550)
     */BL4/X(3),C(3),NNOD
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PP(2,90),PT(2,90)
     */BL2A/NEL,NC(8,110)
     */BLTEOR/ ITEOR
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */btop/ itope(8),itopu
     */AL1A/LHX
     */AL2/XYTX(2,550)
     */AL4/XX(3),CX(3),ISTRX,PRECX,NNODX
     */AL16/MZDX(110),AMZX(110),NZX,NZ1X
     */AL21/NEPX,NEP1X,LTELX(90),LSTX(90),PPX(2,90),PTX(2,90)
     */AL2A/NCX(8,110)
     */ajint/jnkontx,jnkelx,jkelx(50),jksx(50),jkfx(50)
c
      jnkont=jnkontx
      jnkel=jnkelx
      do i=1,jnkel
         jkel(i)=jkelx(i)
         jks(i)=jksx(i)
         jkf(i)=jkfx(i)
      end do
C
      ITEOR=ITEORX
C
      LH=LHX
      NEL=NELX
      NNOD=NNODX
C
      DO 1 I=1,NNOD
      DO 1 J=1,2
    1 XYT(J,I)=XYTX(J,I)
C
      DO 2 I=1,NEL
      DO 2 J=1,8
    2 NC(J,I)=NCX(J,I)
C
      NZ=NZX
C
      DO 3 I=1,NZ
    3 MZD(I)=MZDX(I)
C
      NEP=NEPX
      NEP1=NEP1X
C
      N0=NEP+NEP1
      DO 4 I=1,N0
      LTEL(I)=LTELX(I)
      LST(I)=LSTX(I)
      DO 4 J=1,2
      PT(J,I)=PTX(J,I)
      PP(J,I)=PPX(J,I)
    4 CONTINUE
C
      RETURN
      END
