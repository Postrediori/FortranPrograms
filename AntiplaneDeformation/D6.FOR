      SUBROUTINE APDAI(A1,A2,A3)
C----------------------------------------------------------------------
C     BõóàCãEHàE ìHàBEPCAãúHõX îìHKñàâ MATEPàAãA
C----------------------------------------------------------------------
      COMMON/MDF/ T0,ES,AA,AN,R1,AK0,ETA,AN1
           A1=AA
           A2=A1*R1
           A3=A1*AK0
      RETURN
      END
      SUBROUTINE APDCOB(XI,ET,IG)
C----------------------------------------------------------------------
C     îOPMàPOBAHàE MATPàñõ ÑàîîEPEHñàPOBAHàü
C----------------------------------------------------------------------
      COMMON
     */BL4/X(3),C(3),NNOD
     */BL7/XY(8,2),AKT(8,8)
     */BL8/BT(2,8),DET,R,Z,yak(2,2)
C
      DIMENSION DNE(2),DNY(2),AJ(2,2),AN(8),XEN(2,8)
C
      X1=1-XI
      X2=1+XI
      E1=1-ET
      E2=1+ET
C
      AN(1)=-X1*E1*(X2+ET)
      AN(3)=X2*E1*(XI-E2)
      AN(5)=X2*E2*(ET-X1)
      AN(7)=-X1*E2*(E1+XI)
      AN(2)=X1*X2*E1
      AN(4)=E1*E2*X2
      AN(6)=X1*X2*E2
      AN(8)=E1*E2*X1
C
      XEN(1,1)=E1*(2*XI+ET)
      XEN(2,1)=X1*(2*ET+XI)
      XEN(1,3)=E1*(2*XI-ET)
      XEN(2,3)=X2*(2*ET-XI)
      XEN(1,5)=E2*(2*XI+ET)
      XEN(2,5)=X2*(2*ET+XI)
      XEN(1,7)=E2*(2*XI-ET)
      XEN(2,7)=X1*(2*ET-XI)
      XEN(1,2)=-XI*E1
      XEN(2,2)=-0.5*X1*X2
      XEN(1,4)=0.5*E1*E2
      XEN(2,4)=-ET*X2
      XEN(1,6)=-XI*E2
      XEN(2,6)=-XEN(2,2)
      XEN(1,8)=-XEN(1,4)
      XEN(2,8)=-ET*X1
C
      DO 10 I=1,7,2
      AN(I)=0.25*AN(I)
      AN(I+1)=0.5*AN(I+1)
      XEN(1,I)=0.25*XEN(1,I)
   10 XEN(2,I)=0.25*XEN(2,I)
C
      IF(IG) 14,14,11
C
   11 R=0.125*X1*X2*E1*E2
      Z=-0.25*XI*E1*E2
      T=-0.25*X1*X2*ET
      AN(IG-1)=AN(IG-1)+R
      AN(IG)=AN(IG)-2*R
      J=IG+1
      IF(J.GT.8)J=J-8
      AN(J)=AN(J)+R
      XEN(1,IG-1)=XEN(1,IG-1)+Z
      XEN(2,IG-1)=XEN(2,IG-1)+T
      XEN(1,IG)=XEN(1,IG)-2*Z
      XEN(2,IG)=XEN(2,IG)-2*T
      XEN(1,J)=XEN(1,J)+Z
      XEN(2,J)=XEN(2,J)+T
C
   14 R=0.0
      Z=0.0
      DO 15 I=1,8
      R=R+AN(I)*XY(I,1)
      Z=Z+AN(I)*XY(I,2)
   15 CONTINUE
C
      DO 20 I=1,2
      DO 20 J=1,2
      AJ(I,J)=0.0
      DO 20 K=1,8
   20 AJ(I,J)=AJ(I,J)+XEN(I,K)*XY(K,J)
C
      do i=1,2
        do j=1,2
           yak(i,j)=aj(i,j)
        end do
      end do
      XS=AJ(1,1)
      YS=AJ(1,2)
      DET=AJ(1,1)*AJ(2,2)-AJ(2,1)*AJ(1,2)
      FK=AJ(1,1)
      AJ(1,1)=AJ(2,2)
      AJ(2,2)=FK
      AJ(1,2)=-AJ(1,2)
      AJ(2,1)=-AJ(2,1)
      DO 25 I=1,2
      DO 25 J=1,2
   25 AJ(I,J)=AJ(I,J)/DET
C
      DO 30 I=1,2
      DO 30 J=1,8
   30 BT(I,J)=0.0
C
      DO 55 K=1,8
      DNE(1)=XEN(1,K)
      DNE(2)=XEN(2,K)
           DO 35 I=1,2
           DNY(I)=0.0
           DO 35 J=1,2
   35      DNY(I)=DNY(I)+AJ(I,J)*DNE(J)
      BT(1,K)=DNY(1)
      BT(2,K)=DNY(2)
   55 CONTINUE
C
      RETURN
      END
      SUBROUTINE APDDAT
C----------------------------------------------------------------------
C    BõáOB èPOñEÑìPõ èOÑÉOTOBKà àCXOÑHõX ÑAHHõX à èOÑÉOTOBKA HAÉPìáKà
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(500),RP2(500),RP(500)
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PT(2,90)
C
      X(1)=-0.57735027
      X(3)=0.0
      X(2)=-X(1)
      C(1)=1.0
      C(2)=1.0
      C(3)=0.0
C
C !!!  XADATA - èPOñEÑìPA PAáÅàBKà OÅãACTà HA KOHEóHõE ùãEMEHTõ  !!!
C
      CALL XADATA
C
      WRITE(6,1)
    1 FORMAT(//22X,'*** áAÑAóA AHTàèãOCKOâ ÑEîOPMAñàà ***')
C
      DO 2 I=1,NZ
    2 AMZ(I)=0.0
C
C      îOPMàPOBAHàE  HAÉPìáKà HA èEPBOM ìóACTKE
C
      IF(NEP.LE.0)  GOTO 3
      DO 5 I=1,NNOD
    5 A1(I)=0.0
           CALL APDFOR(NNOD,0)
      DO 6 I=1,NNOD
    6 RP1(I)=A1(I)
C
C      îOPMàPOBAHàE  HAÉPìáKà HA BTOPOM ìóACTKE
C
    3 IF(NEP1.LE.0)  GOTO 4
      DO 7 I=1,NEP1
      LTEL(I)=LTEL(NEP+I)
      LST(I)=LST(NEP+I)
      PT(1,I)=PT(1,NEP+I)
    7 PT(2,I)=PT(2,NEP+I)
      NEP=NEP1
      DO 8 I=1,NNOD
    8 A1(I)=0.0
           CALL APDFOR(NNOD,0)
      DO 9 I=1,NNOD
    9 RP2(I)=A1(I)
    4 CONTINUE
C
      RETURN
      END
      SUBROUTINE APDDP(II,KK,r,dr,AK,dk)
C----------------------------------------------------------------------
C  èOCTPOEHàE MATPàñõ KACATEãúHOâ ÜECTKOCTà èO TEOPàà MàKPOÑEîOPMAñàà
C               B ìÉãOBOâ TOóKE TPAEKTOPàà HAÉPìÜEHàü
C                     Ñãü AHTàèãOCKOâ ÑEîOPMAñàà
C----------------------------------------------------------------------
      COMMON
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1
     */BL9/ DP(2,2),DD(2,2)
     */BLF/ F(2),A2,AM1
     */BLITA/ M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLTEOR/ ITEOR
     */BLALFA/ ALFA(4,110),AL
C
      DIMENSION G(2,2),G1(2,2),r(2),dr(2)
C
      NTA=0
      IF(ITEOR.EQ.1) THEN
          CALL APIN2(r,dr,dk,    F,G,OM)
          IF(NTA.GT.0) GOTO 10
      END IF
      IF(ITEOR.EQ.2.OR.ITEOR.EQ.3) THEN
          ALMAX=ALFA(KK,II)
          CALL APDIN2M(r,dr,ak,ALMAX,   AL,F,G,OM)
          IF(AL.GT.0) GOTO 10
      END IF
C----------
C     ìèPìÉOE ÑEîOPMàPOBAHàE
C
      DO 20 I=1,2
      DO 20 J=1,2
   20 DP(I,J)=0.0
C
      GOTO 30
C----------
C     èãACTàóECKOE ÑEîOPMàPOBAHàE
C
   10 A1=1.0
      A2=0.0
      A3=0.0
      CALL APDAI(A1,A2,A3)
      AM1=A3/(A1+A3*OM)
C
      DO 40 I=1,2
      DO 40 J=1,2
      G(I,J)=G(I,J)-AM1*F(I)*F(J)
      Z=0.0
      IF(I.EQ.J) Z=1.0
   40 G1(I,J)=Z+A2*G(I,J)/A1
C
      CALL APMATR(G1)
C
      DO 50 I=1,2
      DO 50 J=1,2
      DP(I,J)=0.0
      DO 50 L=1,2
   50 DP(I,J)=DP(I,J)+G1(I,L)*G(L,J)/A1
C
C----------
C     îOPMàPOBAHàE MATPàñõ CBüáà MEÜÑì TEHáOPOM HAèPüÜEHàü à TEHáOPOM
C         ÑEîOPMAñàà Ñãü ìèPìÉOÉO à èãACTàóECKOÉO ÑEîOPMàPOBAHàü
   30 DO 60 I=1,2
      DO 60 J=1,2
      Z=0.0
      IF(I.EQ.J) Z=1.0
      DD(I,J)=T0/ES*(Z-DP(I,J))
   60 CONTINUE
      RETURN
      END
      SUBROUTINE APDDEF(U,E)
C----------------------------------------------------------------------
C    BõóàCãEHàE ÑEîOPMAñàâ B ìáãAX àHTEÉPàPOBAHàü
C----------------------------------------------------------------------
      COMMON /BL8/B(2,8),DET,R,Z,yak(2,2)
C
      DIMENSION U(8),E(2)
C
      DO 1 I=1,2
      E(I)=0.0
      DO 1 J=1,8
    1 E(I)=E(I)+B(I,J)*U(J)
C
      RETURN
      END
      SUBROUTINE APDES(ii,kk,DE,  DS,DR,de1,DK)
C----------------------------------------------------------------------
C èPOñEÑìPA èOCTPOEHàü CKOPOCTEâ àáMEHEHàü HAèPüÜEHàâ-DS(2)
C    èAPAMETPA KAèèA-DK à AKTàBHõX HAèPüÜEHàâ -DR(2)
C        èO áAÑAHHõM CKOPOCTüM ÑEîOPMAñàà-DE(2)
C----------------------------------------------------------------------
      DIMENSION DS(2),DE(2),DR(2),DE1(2)
C
      COMMON
     */BL9/DP(2,2),DD(2,2)
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1
     */BLF/F(2),A2,AM1
     */BLITA/ M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLALFA/ ALFA(4,110),AL
     */BLTEOR/ ITEOR
C
C      ìèPìÉOE ÑEîOPMàPOBAHàE
      IF(ITEOR.EQ.1.AND.NTA.GT.0) GOTO 10
      IF((ITEOR.EQ.2.OR.ITEOR.EQ.3).AND.AL.GE.alfa(kk,ii)) GOTO 10
      DO 20 I=1,2
      DS(I)=DE(I)*T0/ES
      DR(I)=DE(I)
      DE1(I)=0.0
   20 CONTINUE
      DK=0.0
      RETURN
C
C     èãACTàóECKOE ÑEîOPMàPOBAHàE
C
   10 DK=0.0
      DO 30 I=1,2
      DE1(I)=0.0
      DO 40 J=1,2
   40 DE1(I)=DE1(I)+DP(I,J)*DE(J)
      DR(I)=DE(I)-A2*DE1(I)
      DS(I)=(DE(I)-DE1(I))*T0/ES
   30 DK=DK+F(I)*DR(I)*AM1
      RETURN
      END
      SUBROUTINE APDDGN(IG,JIN)
C----------------------------------------------------------------------
C èPEOÅPAáOBAHàE  óETõPEXìÉOãúHõX  KOHEóHõX  ùãEMEHTOB  B  TPEìÉOãúHõE
C----------------------------------------------------------------------
      DIMENSION JIN(8)
C
      IG=0
C
      DO 10 I=1,7,2
      IF(JIN(I)-JIN(I+1))10,20,10
   10 CONTINUE
C
      RETURN
C
   20 IG=I+5
      IF(IG.GT.8) IG=IG-8
C
      RETURN
      END
      SUBROUTINE APDELM(II,IG)
C----------------------------------------------------------------------
C     II - HOMEP ùãEMEHTA
C     IG - èPàáHAK BõPOÜÑEHHOCTà ùãEMEHTA
C----------------------------------------------------------------------
      COMMON
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL7/XY(8,2),AKT(8,8)
     */BL15/JIN(8)
C
      DO 1 K=1,8
   1  JIN(K)=NC(K,II)
C
      CALL APDDGN(IG,JIN)
C
      DO 2 K=1,8
      J=JIN(K)
      XY(K,1)=XYT(1,J)
   2  XY(K,2)=XYT(2,J)
C
      RETURN
      END
      SUBROUTINE APDFOR(LL,IFF)
C----------------------------------------------------------------------
C    îOPMàPOBAHàE MATPàñõ CàCTEMõ à EE èPABOâ óACTà
C      ECãà IFF.EQ.0 , TO BõóàCãüETCü èPABAü óACTú
C      ECãà IFF.EQ.1 , TO îOPMàPìETCü MATPàñA ÜECTKOCTà CàCTEMõ
C----------------------------------------------------------------------
      COMMON
     */BL1/A(60000)
     */BL1A/LH,A1(1100)
     */BL2A/NEL,NC(8,110)
     */BL7/XY(8,2),AK(8,8)
     */BL15/JIN(8)
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PT(2,90)
     */BLHT/HPT(8)
C
      IF(IFF.EQ.0) GOTO 10
C
      JJJ=LL*LH
      DO 20 I=1,JJJ
   20 A(I)=0.
C
   10 IP=1
      DO 100 JJ=1,NEL
C
      DO 30 K=1,8
   30 HPT(K)=0.0
C
      CALL APDELM(JJ,IG)
C
      IF(IFF.EQ.0) GOTO 40
C
      CALL APDST(IG,JJ)
C
C     îOPMàPOBAHàE MATPàñõ ÜECTKOCTà CàCTEMõ
C
      DO 50 I=1,8
      DO 50 J=1,8
      I1=JIN(I)
      J1=JIN(J)
      IF(J1.LT.I1) GOTO 50
      M=(I1-1)*LH+J1-I1+1
      A(M)=A(M)+AK(I,J)
   50 CONTINUE
      GOTO 100
C
C     îOPMàPOBAHàE èPABOâ óACTà OT PACèPEÑEãEHHOâ HAÉPìáKà
C
   40 IF(NEP.LE.0)  GOTO 100
      DO 60 J=1,5
      IF(LTEL(IP)-JJ) 70,80,70
   80 CALL APDPR(LST(IP),PT(1,IP),PT(2,IP),HPT)
   60 IP=IP+1
C
   70 DO 90 J=1,8
      M=JIN(J)
   90 A1(M)=A1(M)+HPT(J)
C
  100 CONTINUE
C
C        ÜECTKOE áAKPEèãEHàE
C
      DO 95 I=1,NZ
      IF(IFF.LE.0)  GOTO 85
      M=(MZD(I)-1)*LH+1
      J=M+LH-1
      DO 86 K=M,J
   86 A(K)=0.0
      A(M)=1.0
   85 M=MZD(I)
      A1(M)=0.0
   95 CONTINUE
      IF(IFF.LE.0) GOTO 115
      DO 116 I=1,LL
      DO 96 J=1,NZ
      IF(I.LE.MZD(J)-LH.OR.I.GE.MZD(J)) GOTO 96
      M=(I-1)*LH+MZD(J)-I+1
      A(M)=0.0
   96 CONTINUE
  116 CONTINUE
  115 RETURN
      END
C
      SUBROUTINE APDNOD(W,NNOD)
C----------------------------------------------------------------------
C    èPOñEÑìPA ùKCTPAèOãüñàà èOãEâ áAÑAHHõX B ìáãAX àHTEÉPàPOBAHàü
C                   HA ìáãõ KOHEóHOÉO ùãEMEHTA
C----------------------------------------------------------------------
      COMMON
     */BL2A/NEL,NC(8,110)
C
      DIMENSION
     *W(550),SN(550),LK(550),NPT(8),T(4,4)
C
      DATA T/1.8660254,-0.5,0.1339746,-0.5,
     *       -0.5,0.1339746,-0.5,1.8660254,
     *       -0.5,1.8660254,-0.5,0.1339746,
     *       0.1339746,-0.5,1.8660254,-0.5/
C
      DO 100 I=1,NNOD
      SN(I)=0.0
  100 LK(I)=0
C
      DO 300 K=1,NEL
        DO 1 I=1,8
    1   NPT(I)=NC(I,K)
      DO 300 I=1,4
      N=NPT(2*I-1)
      C=0.0
        DO 200 J=1,4
        II=4*(K-1)+J
  200   C=C+T(I,J)*W(II)
      SN(N)=SN(N)+C
  300 LK(N)=LK(N)+1
C
      DO 600 K=1,NEL
        DO 7 I=1,8
    7   NPT(I)=NC(I,K)
      DO 600 I=1,4
      N=NPT(2*I-1)
      J=2*I+1
      IF(J.GT.8)J=1
      N2=NPT(J)
      N1=NPT(2*I)
      SN(N1)=(SN(N)/LK(N)+SN(N2)/LK(N2))/2.
  600 LK(N1)=1
C
      DO 800 I=1,NNOD
  800 W(I)=SN(I)/LK(I)
C
      RETURN
      END
      SUBROUTINE APDPR(JS,P1,P2,AH)
C----------------------------------------------------------------------
C       BõóàCãEHàE BEKTPA BHEòHEâ HAÉPìáKà ùãEMEHTA AH(8)
C           (JS - HOMEP ìáãA B ãOKAãúHOâ HìMEPAñàà)
C----------------------------------------------------------------------
      COMMON /BL7/XY(8,2) ,AKT(8,8)
C
      DIMENSION XP(3,2),AN(3),AH(8),XIN(3),X(3),C(3)
C
      DO 1 I=1,8
    1 AH(I)=0.0
C
      X(1)=-0.77459667
      X(2)=-X(1)
      X(3)=0.0
      C(1)=0.55555556
      C(2)=0.55555556
      C(3)=0.88888889
C
      DO 10 I=1,3
      J=I+JS-1
      IF(J.GT.8) J=1
      XP(I,1)=XY(J,1)
   10 XP(I,2)=XY(J,2)
C
      DO 90 M=1,3
      XI=X(M)
      X1=1.0-XI
      X2=1.0+XI
      AN(1)=-0.5*X1*XI
      AN(2)=X1*X2
      AN(3)=0.5*X2*XI
      XIN(1)=0.5*(2.0*XI-1.0)
      XIN(2)=-2.0*XI
      XIN(3)=0.5*(2.0*XI+1.0)
      P=P1+0.5*(P2-P1)*(1.0+XI)
      R=0.0
      YS=0.0
      XS=0.0
           DO 20 I=1,3
           R=R+AN(I)*XP(I,1)
           YS=YS+XIN(I)*XP(I,2)
   20      XS=XS+XIN(I)*XP(I,1)
                DS=SQRT(XS*XS+YS*YS)
                DO 80 I=1,3
                J=(I+JS-1)
                IF(J.GT.8) J=1
                ANI=AN(I)*DS*P*C(M)
  80            AH(J)=AH(J)+ANI
  90  CONTINUE
C
      RETURN
      END
      FUNCTION APDRR(X,XI,ET)
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
                         APDRR=0.0
                         DO 1 I=1,8
    1                    APDRR=APDRR+X(I)*AN(I)
      RETURN
      END
      SUBROUTINE APDSE(K,DA,LA,K1,A0)
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
C
      RETURN
      END
      SUBROUTINE APDSOL(N)
C----------------------------------------------------------------------
C  èPOñEÑìPA PEòEHàü CàCTEMõ Ñàî.ìPABHEHàâ METOÑOM ùâãEPA
C----------------------------------------------------------------------
      COMMON
     */BLS/SUM(5500),FP(5500)
C
      CALL APDFF(N)
C
      DO 1 I=1,N
    1 SUM(I)=SUM(I)+FP(I)
C
      RETURN
      END
      SUBROUTINE APDST(IG,II)
C----------------------------------------------------------------------
C     îOPMàPOBAHàE MATPàñõ ÜECTKOCTà ùãEMEHTA AKT(8,8)
C----------------------------------------------------------------------
      COMMON
     */BL4/X(3),C(3),NNOD
     */BL7/XY(8,2),AKT(8,8)
     */BL8/B(2,8),DET,R,Z,yak(2,2)
     */BL9/DPO(2,2),DD(2,2)
C
      DIMENSION CC(8,2)
C
      DO 10 I=1,8
      DO 10 J=1,8
   10 AKT(I,J)=0.0
C
      KK=1
C
      DO 20 M=1,2
      DO 20 N=1,2
C
      CALL APDCOB(X(M),X(N),IG)
      CALL APDDDD(II,KK)
C
      DO 30 I=1,8
      DO 30 J=1,2
      CC(I,J)=0.0
      DO 30 K=1,2
   30 CC(I,J)=CC(I,J)+B(K,I)*DD(K,J)
C
      DO 40 I=1,8
      DO 40 J=1,8
      DO 40 K=1,2
   40 AKT(I,J)=AKT(I,J)+CC(I,K)*B(K,J)*C(M)*C(N)*DET
C
      KK=KK+1
C
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE APDSTR(K)
C----------------------------------------------------------------------
C   BõóàCãEHàE HAèPüÜEHàâ. K- HOMEP ìáãA àHTEÉPàPOBAHàü
C----------------------------------------------------------------------
      COMMON
     */BL9/DPO(2,2),DD(2,2)
     */BL14/EPS(2,4),SIG(2,4),UL(8)
C
      DO 1 I=1,2
      SIG(I,K)=0.0
      DO 1 J=1,2
    1 SIG(I,K)=SIG(I,K)+EPS(J,K)*DD(I,J)
C
      RETURN
      END
      SUBROUTINE APDSUB(L1,N1,I1)
C
      COMMON
     */BL1/A(60000)
     */BL1A/LH,A1(1100)
C
      IDK=(N1-1)*LH+1
      ID1=IDK-I1*LH
      M1=I1+1
      N2=L1+N1
      MM=LH-1
C
         DO 40 I=M1,N1
            S=0.0
            IR=ID1
            DO 10 M=1,MM
               IR=IR+LH
               IF(IR-IDK)5,5,20
    5          II=ID1+M
               JJ=N2-I+M
               S=S+A(II)*A1(JJ)
   10       CONTINUE
   20       II=N2-I
            A1(II)=(A1(II)-S)/A(ID1)
   30       ID1=ID1-LH
   40    CONTINUE
C
      RETURN
      END
      SUBROUTINE APDTRT(L1,N1,N,IB)
C
      COMMON/BL1/A(60000)
     *      /BL1A/LH,A1(1100)
C
      IDK=(N1-1)*LH+1
      ID=1
      IH1=LH-1
      L2=L1+N-1
C
      DO 40 I=L1,L2
         ID1=ID
         DO 20 M1=1,IH1
            ID1=ID1+LH
            IF(ID1-IDK) 5,5,30
    5       IF(IB) 7,7,15
    7       JJ=ID+M1
            IF(A(JJ).EQ.0.0)GO TO 15
            MM=IH1-M1+1
            DO 10 M=1,MM
               II=ID1+M-1
               KK=JJ+M-1
               A(II)=A(II)-A(JJ)/A(ID)*A(KK)
   10       CONTINUE
   15       II=I+M1
            JJ=ID+M1
            A1(II)=A1(II)-A(JJ)/A(ID)*A1(I)
   20    CONTINUE
   30 ID=ID+LH
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE XABLOC(NELX,ITEORX)
C----------------------------------------------------------------------
C     èEPEÑAóA ÑAHHõX óEPEá ÅãOKà COMMON
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2/XYT(2,550)
     */BL4/X(3),C(3),NNOD
     */BL16/MZD(110),AMZ(110),NZ,NZ1
     */BL21/NEP,NEP1,LTEL(90),LST(90),PT(2,90)
     */BL2A/NEL,NC(8,110)
     */BLTEOR/ ITEOR
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */AL1A/LHX
     */AL2/XYTX(2,550)
     */AL4/XX(3),CX(3),ISTRX,PRECX,NNODX
     */AL16/MZDX(110),AMZX(110),NZX,NZ1X
     */AL21/NEPX,NEP1X,LTELX(90),LSTX(90),PPX(2,90),PTX(2,90)
     */AL2A/NCX(8,110)
     */ajint/jnkontx,jnkelx,jkelx(50),jksx(50),jkfx(50)
C
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
      PT(J,I)=PPX(J,I)
    4 CONTINUE
C
      RETURN
      END
      SUBROUTINE XABND(K,NNOD,MM,NN,NUF,IM,IPROB)
C----------------------------------------------------------------------
C    èOÑèPOÉPAMMA îOPMàPOBAHàü ÉPAHàóHõX ìCãOBàâ èO èEPEMEôEHàüM
C----------------------------------------------------------------------
      COMMON
     */AL16/MZD(110),AMZ(110),NZ,NZ1
C
      DIMENSION MM(4),NN(5),NUF(12,30)
C
      IF(K.EQ.1.OR.K.EQ.2) GO TO 10
      NZ=1
      IF(IM.EQ.0) NCENTR=(NNOD-1-2*MM(1)-1-7)/2+1
      IF(IM.EQ.1) NCENTR=(NNOD-1-2*MM(1)-1-9)/2+1
      MZD(1)=NCENTR
      IF(IPROB.EQ.3) MZD(1)=1
      GO TO 50
C
   10 IF(K.EQ.2) GO TO 20
      NZ=0
      DO 1 L=1,3
    1 NZ=NZ+MM(L)
      NZ=2*NZ+1
      DO 2 L=1,NZ
    2 MZD(L)=L
      GO TO 50
C
   20 NZ=0
      DO 3 L=1,4
    3 NZ=NZ+(2*NN(L)+1)
      N=0
      DO 4 L=1,4
      NNN=NN(L)*2+1
      NF=3*(L-1)+1
      DO 4 M=1,NNN
      N=N+1
    4 MZD(N)=NUF(NF,M)
C
   50 RETURN
      END
