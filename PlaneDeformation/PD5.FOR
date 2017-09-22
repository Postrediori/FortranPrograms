      SUBROUTINE PDDDD(II,KK)
C----------------------------------------------------------------------
C    èOCTPOEHàE MATPàñ ÜECTKOCTà TEOPàà MàKPOÑEîOPMAñàà
C     II-HOMEP ùãEMEHTA, KK- HOMEP ìáãA àHTEÉPàPOBAHàü
C----------------------------------------------------------------------
      COMMON
     */BLS/ SUM(10000),FP(10000)
     */BL2A/NEL,NC(8,110)
     */BL4/ X(3),C(3),NNOD
     */blita/m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
     */bltcom/tcom(1500000)
      DIMENSION R(3),DR(3)
C
      NN=2*NNOD
      N1=NN+(II-1)*32+(KK-1)*8
      DO I=1,3
        DR(I)=FP(N1+4+I)
        R(I)=SUM(N1+4+I)
      END DO
      DK=FP(N1+8)
      N3=(II-1)*4+KK
      DO I=1,LP
        T(I)=TCOM(LP*(N3-1)+I)
      END DO
C
      CALL PDDDP(II,KK,r,dr,dk)
      RETURN
      END SUBROUTINE PDDDD
C*****
      SUBROUTINE PDDDP(II,KK, R,DR,DK)
C----------------------------------------------------------------------
C       èOCTPOEHàE MATPàñõ KACATEãúHOâ ÜECTKOCTà
C         B ìÉãOBOâ TOóKE TPAEKTOPàà HAÉPìÜEHàü
C       Ñãü èãOCKOâ áAÑAóà TEOPàà MàKPOÑEîOPMAñàà
C----------------------------------------------------------------------
      COMMON
     */BL9/ DP(4,4),DD(3,3),INS
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1,ANU
     */BLF/F(3),A2,AM1
     */BLTEOR/ITEOR
     */blita/m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
C
      DIMENSION G(3,3),D1(3,3),G1(3,3),R(3),DR(3)
C
      CALL PDIN3(II,KK,r,dr,dk,    F,G,OM)
      IF(NTA.GT.0) GOTO 10
C----------
C     ìèPìÉOE ÑEîOPMàPOBAHàE
C
      DO I=1,4
        DO J=1,4
          DP(I,J)=0.0
        END DO
      END DO
C
      DO I=1,3
      DO J=1,3
        Z=0.0
        IF(I.EQ.J) Z=1.0
          DD(I,J)=T0/ES*Z
        END DO
      END DO
C
      CALL PDGUK
C      IF(II.EQ.1.AND.KK.EQ.1) PRINT 789,((dp(i,j),j=1,4),i=1,4)
C----------
      RETURN
C----------
C     èãACTàóECKOE ÑEîOPMàPOBAHàE Ñãü TEOPàà MàKPOÑEîOPMAñàà
C
   10 A1=AA
      A2=A1*RR0
      A3=A1*AK0
      AM1=A3/(A1+A3*OM)
C
C     (BõóàCãEHàE MATPàñõ G B èOÑBàÜHOM ÑEBàATOPHOM ÅAáàCE)
C
      DO I=1,3
        DO J=1,3
          G(I,J)=G(I,J)-AM1*F(I)*F(J)
        END DO
      END DO
C
      DO I=1,3
        DO J=1,3
          Z=0.0
          IF(I.EQ.J) Z=1.0
          G1(I,J)=Z+A2*G(I,J)/A1
        END DO
      END DO
C
      CALL MATR3(G1)
C
      DO I=1,3
        DO J=1,3
          D1(I,J)=0.0
          DO LL=1,3
            D1(I,J)=D1(I,J)+G1(I,LL)*G(LL,J)/A1
          END DO
        END DO
      END DO
C
      DO I=1,3
        DO J=1,3
          Z=0.0
          IF(I.EQ.J) Z=1.0
          DD(I,J)=T0/ES*(Z-D1(I,J))
        END DO
      END DO
C
C (îOPMàPOBAHàE MATPàñõ CBüáà MEÜÑì TEHáOPAMà HAèPüÜEHàü à ÑEîOPMAñàà
C                 Ñãü CãìóAü èãOCKOâ ÑEîOPMAñàà)
C
      OK=T0/ES*(1.0+ANU)/(1.0-2*ANU)/3.0
      DP(1,1)=(OK+DD(1,1)*2./3.)
      DP(1,2)=(OK-DD(1,1)*1./3.+DD(1,2)*1./1.732)
      DP(1,3)=(DD(1,3)*2./1.732)
      DP(1,4)=(OK-DD(1,1)*1./3.-DD(1,2)*1./1.732)
      DP(2,1)=(OK-DD(1,1)*1./3.+DD(2,1)*1./1.732)
      DP(2,2)=(OK-DD(2,1)*1./3.464+DD(2,2)*1./2.+DD(1,1)*1./6.-
     *        DD(1,2)*1./3.464)
      DP(2,3)=(DD(2,3)-DD(1,3)*1./1.732)
      DP(2,4)=(OK-DD(2,1)*1./3.464-DD(2,2)*1./2.+DD(1,1)*1./6.+
     *        DD(1,2)*1./3.464)
      DP(3,1)=(DD(3,1)*1./1.732)
      DP(3,2)=(-DD(3,1)*1./3.464+DD(3,2)*1./2.)
      DP(3,3)=DD(3,3)
      DP(3,4)=(-DD(3,1)*1./3.464-DD(3,2)*1./2.)
      DP(4,1)=(OK-DD(1,1)*1./3.-DD(2,1)*1./1.732)
      DP(4,2)=(OK+DD(1,1)*1./6.-DD(1,2)*1./3.464+DD(2,1)*1./3.464-
     *        DD(2,2)*1./2.)
      DP(4,3)=(-DD(1,3)*1./1.732-DD(2,3))
      DP(4,4)=(OK+DD(1,1)*1./6.+DD(1,2)*1./3.464+DD(2,1)*1./3.464+
     *        DD(2,2)*1./2.)
C
C  (MATPàñA CBüáà DP(I,J) Ñãü CãìóAü èãOCKOÉO HAèPüÜEHHOÉO COCTOüHàü)
C
      IF(INS.EQ.0) GO TO 20
      DO I=1,3
        DO J=1,3
          DP(I,J)=DP(I,J)-DP(I,4)*DP(4,J)/DP(4,4)
        END DO
      END DO
      DO I=1,3
        DP(4,I)=0.0
        DP(I,4)=0.0
      END DO
C----------
   20	CONTINUE
c      IF(II.EQ.1.AND.KK.EQ.1) PRINT 789,((dp(i,j),j=1,4),i=1,4)
 789  FORMAT(2x,4f12.2)
      RETURN
      END SUBROUTINE PDDDP
C*****
      SUBROUTINE PDGUK
C----------------------------------------------------------------------
C  ìèPìÉAü CBüáú MEÜÑì TEHáOPAMà HAèPüÜEHàü à ÑEîOPMAñàà
C----------------------------------------------------------------------
      COMMON
     */MDF/ T0,ES,AA,AN,RR0,AK0,ETA,AN1,ANU
     */BL9/ DP(4,4),DD(3,3),INS
C
      E=T0/ES*(1.0+ANU)
      ANU1=ANU
      IF(INS) 20,20,30
C     <= 0
  20  ANU1=ANU/(1.0-ANU)
C     > 0
  30  DP(1,1)=1.0
      DP(1,2)=ANU1
      DP(1,4)=ANU1
      DP(2,2)=1.0
      DP(2,1)=ANU1
      DP(2,4)=ANU1
      DP(3,3)=(1.0-ANU1)/2.0
      DP(4,1)=ANU1
      DP(4,2)=ANU1
      DP(4,4)=1.0
C
      DO I=1,4
        DO J=1,4
          DP(I,J)=DP(I,J)*E/((1.0+ANU)*(1.0-ANU1))
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PDGUK
C*****
      SUBROUTINE PDDES(DE,  DS,DR,de1,DK)
C----------------------------------------------------------------------
C  èPOñEÑìPA èOCTPOEHàü CKOPOCTEâ àáMEHEHàü HAèPüÜEàâ-DS(4)
C     èAPAMETPA KAèèA-DK à AKTàBHõX HAèPüÜEHàâ -DR(3)
C         èO áAÑAHHõM CKOPOCTüM ÑEîOPMAñàà-DE(4)
C----------------------------------------------------------------------
      DIMENSION DS(4),DE(4),DR(3),DE1(3),DEE(3)
C
      COMMON
     */BL9/DP(4,4),DD(3,3),INS
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1,ANU
     */BLF/ F(3),A2,AM1
     */BLTEOR/ITEOR
     */blita/m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
C
C     ìèPìÉOE ÑEîOPMàPOBAHàE
C
      IF(NTA.GT.0) GOTO 10
      do i=1,3
         de1(i)=0.0
      end do
      DO I=1,4
        DS(I)=0.0
        DO J=1,4
          DS(I)=DS(I)+DP(I,J)*DE(J)
        END DO
      END DO
      IF(INS.EQ.1) DS(4)=0.0
           CALL PDTEVE(DE,DR)
      DK=0.0
      RETURN
C
C     èãACTàóECKOE ÑEîOPMàPOBAHàE
C
   10 DK=0.0
      DO I=1,3
        DE1(I)=0.0
        CALL PDTEVE(DE,DEE)
        DO J=1,3
          DE1(I)=DE1(I)-DD(I,J)*DEE(J)
        END DO
        DR(I)=DEE(I)-A2*DE1(I)
        DK=DK+F(I)*DR(I)*AM1
      END DO
      DO I=1,4
        DS(I)=0.0
        DO J=1,4
          DS(I)=DS(I)+DP(I,J)*DE(J)
        END DO
      END DO
c      PRINT 789,((dp(i,j),j=1,4),i=1,4)
 789  FORMAT(2x,'ds',4f12.2)
      IF(INS.EQ.1) DS(4)=0.0
      RETURN
      END SUBROUTINE PDDES
C*****
      SUBROUTINE PDTEVE(TEN,VEC)
C----------------------------------------------------------------------
C  èPOñEÑìPA èEPEBOÑA KOMèOHEHT TEáOPA HAèPüÜEHàü (ÑEîOPMAñàà) B
C           KOMèOHEHTõ BEKTOPA HAèPüÜEHàü (ÑEîOPMAñàà)
C           èPà OÅOôEHHO-èãOCKOM HAèPüÜEHHOM COCTOüHàà
C
C         TEN(1)=SIGX,TEN(2)=SIGY,TEN(3)=SIGXY,TEN(4)=SIGZ
C        (TEN(1)=EPSX,TEN(2)=EPSY,TEN(3)=EPSXY,TEN(4)=EPSZ)
C----------------------------------------------------------------------
      DIMENSION TEN(4),VEC(3)
C
      TEN0=(TEN(1)+TEN(2)+TEN(4))/3.
      VEC(1)=1.732/1.414*(TEN(1)-TEN0)
      VEC(2)=1.414/2.*(TEN(1)+2*TEN(2)-3*TEN0)
      VEC(3)=1.414*TEN(3)/2.0
C
      RETURN
      END SUBROUTINE PDTEVE
