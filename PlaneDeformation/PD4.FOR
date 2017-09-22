      SUBROUTINE PDELM(II,IG)
C----------------------------------------------------------------------
C       II - HOMEP ‹EMEHTA
C       IG - Pˆ‡HAK B›PO†„EHHOCTˆ ‹EMEHTA
C----------------------------------------------------------------------
      COMMON
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL15/LEC(16),JIN(8)
     */BL7/XY(8,2),AK(16,16)
C
      DO K=1,8
        JIN(K)=NC(K,II)
      END DO
      DO I=1,8
        K=2*I
        LEC(K-1)=JIN(I)*2-1
        LEC(K)=LEC(K-1)+1
      END DO
      CALL PDDGNR(IG,JIN)
      DO K=1,8
        J=JIN(K)
        XY(K,1)=XYT(1,J)
        XY(K,2)=XYT(2,J)
      END DO
      RETURN
      END SUBROUTINE PDELM
C*****
      SUBROUTINE PDSTIF(II,IG)
C----------------------------------------------------------------------
C     ”OPMˆPOBAHˆE MATPˆ–› †ETKOCTˆ ‹EMEHTA  AK(16,16)
C----------------------------------------------------------------------
      COMMON
     */BL4/X(3),C(3),NNOD
     */BL8/B(4,16),DET,R,Z,AN(8),yak(2,2)
     */BL9/DP(4,4),DD(3,3),INS
     */BL7/XY(8,2),AK(16,16)
C
      DIMENSION CC(16,4)
C
C    ”OPMˆPOBAHˆE MATPˆ–› †ECTKOCTˆ
C
      DO I=1,16
        DO J=1,16
          AK(I,J)=0.0
        END DO
      END DO
C
      LD=4
      IF(INS.EQ.1) LD=3
      KK=1
C
      DO M=1,2
        DO N=1,2
C
          CALL PDCOB(X(M),X(N),IG)
          CALL PDDDD(II,KK)
c      IF(II.EQ.1.AND.KK.EQ.1) PRINT 789,((dp(i,j),j=1,4),i=1,4)
c 789  FORMAT(2x,'fo',4f12.2)
C
          KK=KK+1
C
          DO I=1,16
            DO J=1,LD
              CC(I,J)=0.0
              DO K=1,LD
                CC(I,J)=CC(I,J)+B(K,I)*DP(K,J)
              END DO
            END DO
          END DO
C
          DO I=1,16
            DO J=1,16
              DO K=1,LD
                AK(I,J)=AK(I,J)+CC(I,K)*B(K,J)*C(M)*C(N)*DET
              END DO
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE PDSTIF
C*****
      SUBROUTINE PDCOB(XI,ET,IG)
C----------------------------------------------------------------------
C     ”OPMˆPOBAHˆE MATPˆ–› „ˆ””EPEH–ˆPOBAHˆŸ
C----------------------------------------------------------------------
      COMMON
     */BL7/XY(8,2),AK(16,16)
     */BL8/B(4,16),DET,R,Z,AN(8),yak(2,2)
C
      DIMENSION XEN(2,8),DNE(2),DNY(2),AJ(2,2)
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
      DO I=1,7,2
        AN(I)=0.25*AN(I)
        AN(I+1)=0.5*AN(I+1)
        XEN(1,I)=0.25*XEN(1,I)
        XEN(2,I)=0.25*XEN(2,I)
      END DO
C
      IF(IG) 14,14,11
C     IG > 0
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
C     IG <= 0
   14 R=0.0
      Z=0.0
      T=0.0
      DT=0.0
C     END IF
      DO I=1,8
        R=R+AN(I)*XY(I,1)
        Z=Z+AN(I)*XY(I,2)
      END DO
C
      DO I=1,2
        DO J=1,2
          AJ(I,J)=0.0
          DO K=1,8
            AJ(I,J)=AJ(I,J)+XEN(I,K)*XY(K,J)
          END DO
        END DO
      END DO
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
C
      DO I=1,2
        DO J=1,2
          AJ(I,J)=AJ(I,J)/DET
        END DO
      END DO
C
      DO I=1,4
        DO J=1,16
          B(I,J)=0.0
        END DO
      END DO
C
      DO K=1,8
        DNE(1)=XEN(1,K)
        DNE(2)=XEN(2,K)
        DO I=1,2
          DNY(I)=0.0
          DO J=1,2
            DNY(I)=DNY(I)+AJ(I,J)*DNE(J)
          END DO
        END DO
        I1=2*K-1
        I2=2*K
C
        B(1,I1)=DNY(1)
        B(2,I2)=DNY(2)
        B(3,I1)=DNY(2)
        B(3,I2)=DNY(1)
      END DO
C
   70 RETURN
      END SUBROUTINE PDCOB
C*****
      SUBROUTINE PDDGNR(IG,JIN)
C----------------------------------------------------------------------
C PEOPA‡OBAHˆE  —ET›PEX“ƒO‹œH›X  KOHE—H›X  ‹EMEHTOB  B  TPE“ƒO‹œH›E
C----------------------------------------------------------------------
      DIMENSION JIN(8)
C
      IG=0
C
      DO 10 I=1,7,2
      IF(JIN(I)-JIN(I+1))10,20,10
C     <0 || >0
   10 CONTINUE
C
      RETURN
C     = 0
   20 IG=I+5
      IF(IG.GT.8) IG=IG-8
C
      RETURN
      END SUBROUTINE PDDGNR
C*****
      SUBROUTINE PDPRES(JS,P1,P2,G1,G2,AH)
C----------------------------------------------------------------------
C        ”OPMˆPOBAHˆE HAƒP“‡Kˆ, „ECTB“ž™E‰ HA ‹EMEHT
C     P1,P2 - HOPMA‹œHAŸ HAƒP“‡KA B “ƒ‹OB›X “‡‹AX ‹EMEHTA
C     G1,G2 - KACATE‹œHAŸ HAƒP“‡KA B “ƒ‹OB›X “‡‹AX ‹EMEHTA
C        JS - HOMEP “‡‹A B ‹OKA‹œHO‰ H“MEPA–ˆˆ Pˆ OXO„E
C             POTˆB —ACOBO‰ CTPE‹Kˆ
C----------------------------------------------------------------------
      COMMON /BL7/XY(8,2),AK(16,16)
C
      DIMENSION AH(16),XP(3,2),AN(3),XIN(3),X(3),C(3)
C
      X(1)=-0.77459667
      X(2)=-X(1)
      X(3)=0.0
      C(1)=0.55555556
      C(2)=0.55555556
      C(3)=0.88888889
C
      DO I=1,3
        J=I+JS-1
        IF(J.GT.8) J=1
        XP(I,1)=XY(J,1)
        XP(I,2)=XY(J,2)
      END DO
C
      DO M=1,3
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
        G=G1+0.5*(G2-G1)*(1.0+XI)
        R=0.0
        YS=0.0
        XS=0.0
        DO I=1,3
          R=R+AN(I)*XP(I,1)
          YS=YS+XIN(I)*XP(I,2)
          XS=XS+XIN(I)*XP(I,1)
        END DO
        DO I=1,3
          J=(I+JS-1)*2
          IF(J.GT.16) J=2
          AH(J-1)=AH(J-1)+AN(I)*C(M)*(-YS*P+XS*G)
          AH(J)=AH(J)+AN(I)*C(M)*(XS*P+YS*G)
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PDPRES
C*****
      SUBROUTINE PDSTR(K)
C----------------------------------------------------------------------
C    O„C—ET HAPŸ†EHˆ‰ B “‡‹AX ˆHTEƒPˆPOBAHˆŸ
C     K - HOMEP “‡‹A ˆHTEƒPˆPOBAHˆŸ
C----------------------------------------------------------------------
      COMMON
     */BL9/DP(4,4),DD(3,3),INS
     */BL14/EPS(4,4),SIG(4,4),UL(16)
C
      LD=4
      IF(INS.EQ.1) LD=3
      DO I=1,LD
        SIG(I,K)=0.0
        DO J=1,LD
          SIG(I,K)=SIG(I,K)+EPS(J,K)*DP(I,J)
        END DO
      END DO
C
      IF(INS.EQ.1) SIG(4,K)=0.0
C
      RETURN
      END SUBROUTINE PDSTR
C*****
      SUBROUTINE PDDEF(U,E)
C----------------------------------------------------------------------
C     ”OPMˆPOBAHˆE „E”OPMA–ˆ‰ B “‡‹AX ˆHTEƒPˆPOBAHˆŸ
C----------------------------------------------------------------------
      COMMON
     */BL8/B(4,16),DET,R,Z,AN(8),yak(2,2)
     */BL9/DP(4,4),DD(3,3),INS
     */MDF/T0,ES,AA,ANNN,RR0,AK0,ETA,AN1,ANU
C
      DIMENSION U(16),E(4)
C
      DO I=1,4
        E(I)=0.0
        DO J=1,16
          E(I)=E(I)+B(I,J)*U(J)
        END DO
      END DO
c      print *,'eps',e
C
      IF(INS.EQ.1) E(4)=-ANU*(E(1)+E(2))
C
      RETURN
      END SUBROUTINE PDDEF
C*****
      SUBROUTINE PDNODY(W,NNOD)
C----------------------------------------------------------------------
C    PO–E„“PA KCTPAO‹Ÿ–ˆˆ O‹E‰ ‡A„AHH›X B “‡‹AX ˆHTEƒPˆPOBAHˆŸ
C                   HA “‡‹› KOHE—HOƒO ‹EMEHTA
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
      DO I=1,NNOD
        SN(I)=0.0
        LK(I)=0
      END DO
C
      DO K=1,NEL
        DO I=1,8
          NPT(I)=NC(I,K)
        END DO
        DO I=1,4
          N=NPT(2*I-1)
          C=0.0
          DO J=1,4
            II=4*(K-1)+J
            C=C+T(I,J)*W(II)
          END DO
          SN(N)=SN(N)+C
          LK(N)=LK(N)+1
        END DO
      END DO
C
      DO K=1,NEL
        DO I=1,8
          NPT(I)=NC(I,K)
        END DO
        DO I=1,4
          N=NPT(2*I-1)
          J=2*I+1
          IF(J.GT.8)J=1
          N2=NPT(J)
          N1=NPT(2*I)
          SN(N1)=(SN(N)/LK(N)+SN(N2)/LK(N2))/2.
          LK(N1)=1
        END DO
      END DO
C
      DO I=1,NNOD
        W(I)=SN(I)/LK(I)
      END DO
C
      RETURN
      END SUBROUTINE PDNODY
C*****
      SUBROUTINE PDBSUB(L1,N1,I1)
C
      COMMON
     */BL1/A(400000)
     */BL1A/LH,A1(1100)
C
      IDK=(N1-1)*LH+1
      ID1=IDK-I1*LH
      M1=I1+1
      N2=L1+N1
      MM=LH-1
C
      DO I=M1,N1
        S=0.0
        IR=ID1
        DO M=1,MM
          IR=IR+LH
          IF(IR-IDK)5,5,20
C         > 0
   20     EXIT
C         <= 0
    5     II=ID1+M
          JJ=N2-I+M
          S=S+A(II)*A1(JJ)
        END DO
        II=N2-I
        A1(II)=(A1(II)-S)/A(ID1)
        ID1=ID1-LH
      END DO
C
      RETURN
      END SUBROUTINE PDBSUB
C*****
      SUBROUTINE PDSTRT(L1,N1,N,IB)
C
      COMMON/BL1/A(400000)
     *      /BL1A/LH,A1(1100)
C
      IDK=(N1-1)*LH+1
      ID=1
      IH1=LH-1
      L2=L1+N-1
C
      DO I=L1,L2
        ID1=ID
        DO M1=1,IH1
          ID1=ID1+LH
          IF(ID1-IDK) 5,5,30
C         > 0
   30     EXIT
C         <= 0
    5     IF(IB) 7,7,15
C         <=0
    7     JJ=ID+M1
          IF(A(JJ).EQ.0.0)GO TO 15
          MM=IH1-M1+1
          DO M=1,MM
            II=ID1+M-1
            KK=JJ+M-1
            A(II)=A(II)-A(JJ)/A(ID)*A(KK)
          END DO
C         > 0
   15     II=I+M1
          JJ=ID+M1
          A1(II)=A1(II)-A(JJ)/A(ID)*A1(I)
        END DO
        ID=ID+LH
      END DO
C
      RETURN
      END SUBROUTINE PDSTRT
