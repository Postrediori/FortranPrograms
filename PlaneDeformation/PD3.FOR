      SUBROUTINE PDQQ
C---------------------------------------------------------------------
C     �O��PO�PAMMA ��� B���C�EH�� HA�PAB�����X
C              KOC�H�COB  QQ(4,LP)
C---------------------------------------------------------------------
      DIMENSION   b(3),r(3)
C
      COMMON
     *     /blqq/    QQ(4,27000)
     *     /blita/   m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
      DATA B/0.555556,0.888889,0.555556/,
     *     R/-0.774597,0.0,0.774597/
C
      p=3.14159
      m1=5
      m2=5
      m3=5
      nq=3
      lp=m1*m2*m3*nq**3
C
      PM1=P/M1/2.
      PM2=P/M2/2.
      PM3=P/M3/2.
C
      NTA=0
      OM=0.0
      LL=0
C
      DO N=1,NQ
        DO J=1,NQ
          DO I=1,NQ
            BB=B(I)*B(N)*B(J)*PM1*PM2*PM3*P*2.0
            DO K=1,M1
              TE1=PM1*(R(N)+2.*K-1.)
              DO L=1,M2
                TE2=PM2*(R(J)+2.*L-1.)
                DO MI=1,M3
                  TE3=PM3*(R(I)+2.*MI-1.)
                  LL=LL+1
                  STE1=SIN(TE1)
                  STE2=SIN(TE2)
                  QQ(1,LL)=COS(TE1)
                  QQ(2,LL)=STE1*COS(TE2)
                  QQ(3,LL)=STE1*STE2*COS(TE3)
                  QQ(4,LL)=STE1**3 * STE2**2 * SIN(TE3)*BB
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE PDQQ
C*****
      SUBROUTINE PDIN3(II,KK,R,DR,DK,  F,G,OM)
C--------------------------------------------------------------------
C                �O��PO�PAMMA ��� B���C�EH�� �HTE�PA�OB
C                 �O �OBEPXHOCT� ��T�MEPHO� ���EPC�EP� .
C
C                       O��ACT� �HTE�P�POBAH��
C                             0 < X1 < PI/2
C                             0 < X2 < PI
C                             0 < X3 < PI
C             PA���BAETC� CETKO� HA M1*M2*M3  KBA�PATOB.
C         �HTE�P�POBAH�E - �O �OPM��E �A�CCA C  NQ  ���AM�.
C--------------------------------------------------------------------
      DIMENSION   DR(3),F(3),G(3,3),Q(3),R(3)
C
      COMMON    /blqq/  QQ(4,27000)
      COMMON    /blita/   m1,m2,m3,nq,lp,nta,ita(3375),tt(3375)
C
      DO i=1,3
        f(i)=0.0
        DO j=1,3
          G(I,J)=0.0
        END DO
      END DO
C
      NTA=0
      OM=0.0
      LL=0
C
      DO N=1,NQ
        DO J=1,NQ
          DO I=1,NQ
            DO K=1,M1
              DO L=1,M2
                DO MI=1,M3
                  LL=LL+1
                  Q(1)=QQ(1,LL)
                  Q(2)=QQ(2,LL)
                  Q(3)=QQ(3,LL)
                  DOMG=QQ(4,LL)
                  QR=Q(1)*DR(1)+Q(2)*DR(2)+Q(3)*DR(3)
                  QS=(Q(1)*R(1)+Q(2)*R(2)+Q(3)*R(3))*1.00001
                  IF((QS.LT.tt(LL)).OR.(QR.LT.DK)) CYCLE
C
C                 �OPM�POBAH�E HA�PAB�EH�� AKT�BHO�O M�KPO�E�OPM�POBAH��
C
                  NTA=NTA+1
                  ITA(NTA)=LL
                  OM=OM+DOMG
                  DO K1=1,3
                    FQ=Q(K1)*DOMG
                    F(K1)=F(K1)+FQ
                    DO K2=1,K1
                      G(K1,K2)=G(K1,K2)+Q(K2)*FQ
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      G(3,3)=G(3,3)
      G(1,3)=G(3,1)
      G(2,3)=G(3,2)
      G(1,2)=G(2,1)
      RETURN
      END SUBROUTINE PDIN3
C*****
      SUBROUTINE PDST3(dr,DK)
C---------------------------------------------------------------------
C     �O��PO�PAMMA ��� B���C�EH�� �HTEHC�BHOCT� PA�PE�A���X
C              HA�P��EH�� �OC�E ���OBO� TO�K�
C     O��ACT� �HTE�P�POBAH��  0 < X1 < PI
C                             0 < X2 < PI
C                             0 < X3 < PI
C        PA���BAETC� CETKO� HA (2*M1)*(2*M2) KBA�PATOB.
C----------------------------------------------------------------------
      DIMENSION DR(3),Q(3)
      COMMON
     *        /blqq/  QQ(4,27000)
     *        /blita/ m1,m2,m3,nq,lp,nta,ita(3375),tt(3375)
     *        /mdf/   t0,es,aa,an,rr,ak0,eta,an1
C
      p=3.14159
      LPA=1
      NQ2=NQ/2
      LL=0
      DO N=1,NQ
        DO J=1,NQ
          DO I=1,NQ
            DO K=1,M1
              DO L=1,M2
                DO MI=1,M3
                  LL=LL+1
                  Q(1)=QQ(1,LL)
                  Q(2)=QQ(2,LL)
                  Q(3)=QQ(3,LL)
                  QDr=0.0
                  DO JJ=1,3
                    QDr=QDr+Q(JJ)*Dr(JJ)
                  END DO
C
C  O��ACT� OCTA��H�X HA�PAB�EH��
C
                  tt(LL)=tt(LL)+DK
                  IF(ITA(LPA).NE.LL.OR.LPA.GT.NTA) CYCLE
C
C   O��ACT� AKT�BH�X HA�PAB�EH��
C
                  tt(LL)=tt(LL)+QDr-DK
                  LPA=LPA+1
                  LPR=LP-LL+1
                  tt(LPR)=tt(LPR)-ETA*QDr+ETA*DK
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE PDST3
C*****
      SUBROUTINE MATR3(a)
C----------------------------------------------------------------------
C          �OCTPOEH�E MATP��� O�PATHO� K MATP��E A(3,3)
C               PA�ME�AETC� HA MECTE MATP��� A(3,3)
C----------------------------------------------------------------------
      DIMENSION a(3,3)
      d11=a(2,2)*a(3,3)-a(3,2)**2
      d12=a(3,2)*a(3,1)-a(2,1)*a(3,3)
      d13=a(2,1)*a(3,2)-a(3,1)*a(2,2)
      d=a(1,1)*d11+a(1,2)*d12+a(1,3)*d13
      d22=a(1,1)*a(3,3)-a(3,1)**2
      d23=a(3,1)*a(1,2)-a(1,1)*a(3,2)
      d33=a(1,1)*a(2,2)-a(1,2)**2
      a(1,1)=d11/d
      a(1,2)=d12/d
      a(1,3)=d13/d
      a(2,1)=a(1,2)
      a(2,2)=d22/d
      a(2,3)=d23/d
      a(3,1)=a(1,3)
      a(3,2)=a(2,3)
      a(3,3)=d33/d
      RETURN
      END SUBROUTINE MATR3

