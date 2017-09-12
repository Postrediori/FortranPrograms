C----------------------------------------------------------------------
C             èPOÉPAMMA PEòEHàü ÉPAHàóHõX áAÑAó TEOPàà èãÄëíàóçéëíà
C                      èPà AHTàèãOCKOâ ÑEîOPMAñàà
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(500),RP2(500),RP(500)
     */MDF/T0,ES,AA,AN,RR,AK0,ETA,AN1
     */BLS/SUM(5500),FP(5500)
     */BLITA/M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLTCOM/TCOM(400000)
     */BLTEOR/ITEOR
     */BLALFA/ALFA(4,110),AL
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
     */blj/zjsr,zjupr
C
      DIMENSION A0(3,10),LA(10),DVS(3)
C
      OPEN (8,FILE='$APDFIZ')
      READ (8,10) AA,AN,RR,AK0,ETA,T0,ES,AN1
  10  FORMAT(/8F10.4)
      READ (8,20) K2
  20  FORMAT(/I5/)
      READ (8,30) ((A0(I,J),I=1,3),LA(J),J=1,K2)
  30  FORMAT(3F10.5,I5)
      CLOSE(8)
C
      DO I=1,3
         DO J=1,K2
            A0(I,J)=A0(I,J)*(sqrt(2.0)*T0)
         END DO
      END DO
C
      OPEN(6,FILE='RESULT')
      WRITE (6,40) RR,AN1,AK0,AN,AA,ETA,T0,ES
  40  FORMAT(/33X,'èAPAMETPõ MOÑEãà'
     *   /19X,'RR=',F8.4,'AN1=',F8.4,'AK0=',F8.4,'AN=',F8.4
     *   /19X,'AA=',F8.4,'ETA=',F8.4,'T0=',F8.4,'ES=',F8.4)
C
      CALL APDDAT
      READ(5,1) IBUIL,NEE,NU
    1 FORMAT(///14X,I1/14X,I3/14X,I1)
      CLOSE(5)
C
      if(ITEOR.EQ.0) WRITE (6,101)
      IF(ITEOR.EQ.1) WRITE (6,102)
      IF(ITEOR.EQ.2) WRITE (6,103)
      IF(ITEOR.EQ.3) WRITE (6,104)
  101 FORMAT(/35x,'ìèêìÉéÖ êÖòÖçàÖ')
  102 FORMAT(/24X,'PACóET èO TEOPàà MàKPOÑEîOPMAñàà')
  103 FORMAT(/24X,'PACóET èO ÑÖîéêåÄñàéççéâ íÖéêàà')
  104 FORMAT(/27X,'PACóET èO TEOPàà íÖóÖçàü')
C
      WRITE (6,100) NEL,NNOD,LH
  100 FORMAT(//29X,'OÅôàE  èAPAMETPõ MKù'/
     *19X,'KOãàóECTBO äéçÖóçõï ùãEMEHTOB',6X,'NEL =',I3/
     *19X,'KOãàóECTBO  ìáãOB  Kù ',13X,'NNOD=',I4/
     *19X,'òàPàHA ãEHTõ MATPàñõ',15X,'LH  =',I4//)
C
      IF(JNKEL.NE.0) THEN
         ZJUPR =0
         ZJSR=0
      END IF
C
      DO I=1,NEL
         DO J=1,4
            WWW(J,I)=0.0
            ALFA(J,I)=0.00001
         END DO
      END DO
C
      NSUM=NNOD+20*NEL
      DO 80 I=1,NSUM
      FP(I)=0.0
   80 SUM(I)=1.0E-10
C
      IF(ITEOR.EQ.1) THEN
         CALL APQQ
         LPZ=LP*NEL*4
         DO I=1,LPZ
            TCOM(I)=ES
         END DO
      END IF
C
      CALL APDUPR(A0,K2)
C
      IF(ITEOR.EQ.0) GO TO 200
C
      OPEN (31,FILE='DIAGRAM.DAT')
      OPEN (32,FILE='TRAJEPS1.DAT')
      OPEN (322,FILE='TRAJEPS2.DAT')
      OPEN (323,FILE='TRAJEPS3.DAT')
      OPEN (324,FILE='TRAJEPS4.DAT')
      OPEN (325,FILE='TRAJEPS5.DAT')
      OPEN (326,FILE='TRAJEPS6.DAT')
      OPEN (327,FILE='TRAJEPS7.DAT')
      OPEN (328,FILE='TRAJEPS8.DAT')
      OPEN (33,FILE='TRAJSIG.DAT')
C
      WRITE(6,50)
   50 FORMAT(//27X,'áAÑAHA TPAEKTOPàü HAÉPìÜEHàü'
     *     //25X,' P1        P2         K-BO òAÉOB')
      WRITE(6,60) ((A0(I,J),I=1,2),LA(J),J=1,K2)
   60 FORMAT(/20X,2F10.3,I13)
C
      DO I=1,NEL
         DO J=1,4
            WWW(J,I)=0.0
            ALAMBD(J,I)=0.0
            ALFA(J,I)=0.00001
         END DO
      END DO
C
      NSUM=NNOD+20*NEL
      DO 85 I=1,NSUM
      FP(I)=0.0
   85 SUM(I)=1.0E-10
C
      K1=1
      KP=LA(K2)-1
C
      DO LOAD=1,KP
           PRINT 120,LOAD
  120      FORMAT(/37X,'òAÉ',I5)
                CALL APDSE(LOAD,DVS,LA,K1,A0)
           DO I=1,NNOD
                RP(I)=RP1(I)*DVS(1)+RP2(I)*DVS(2)
           END DO
                CALL APDSOL(NSUM)
                CALL APDRES(0,0,0,0,0,0,IBUIL,NEE,NU)
      END DO
C
      CALL APDRES(1,1,1,1,1,1,IBUIL,NEE,NU)
  200 CONTINUE
      CLOSE(31)
      CLOSE(32)
      CLOSE(322)
      CLOSE(323)
      CLOSE(324)
      CLOSE(325)
      CLOSE(326)
      CLOSE(327)
      CLOSE(328)
      CLOSE(33)
      STOP
      END
C
      SUBROUTINE APQQ
C-----------------------------------------------------------------------
C          BCèOMAÉATEãúHAü èPOñEÑìPA BõóàCãEHàü MACCàBOB QQ à DOM
C        (QQ - MACCàB HAèPABãEHàâ B ÑEBàATOPHOM èPOCTPAHCTBE ALFA)
C-----------------------------------------------------------------------
      COMMON
     */BLITA/M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLQQ/QQ(2,900),DOM(900)
      DIMENSION B(3),R(3)
      DATA B/0.555556,0.888889,0.555556/,
     *     R/-0.774597,0.0,0.774597/,
     *     P/3.14159/
      M1=5
      M2=5
      M3=1
      NQ=3
      LP=M1*M2*M3*NQ**2
      PM1=P/M1/2.
      PM2=P/M2/2.
      PM3=2.0
      LL=0
      DO N=1,NQ
         DO J=1,NQ
            BB=B(N)*B(J)*PM1*PM2*PM3*P*2.0
            DO K=1,M1
               TE1=PM1*(R(N)+2.*K-1.)
               DO L=1,M2
                  TE2=PM2*(R(J)+2.*L-1.)
                  LL=LL+1
                  STE1=SIN(TE1)
                  STE2=SIN(TE2)
                  QQ(1,LL)=COS(TE1)
                  QQ(2,LL)=STE1*COS(TE2)
                  DOM(LL)=STE1**3*STE2**2*BB
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
C
      SUBROUTINE APIN2(r,dr,dk,   F,G,OM)
C-----------------------------------------------------------------------
C            èOÑèPOÉPAMMA Ñãü BõóàCãEHàü àHTEÉPAãOB
C                 èO èOBEPXHOCTà ÉàèEPCîEPõ
C-----------------------------------------------------------------------
      COMMON
     */BLITA/M1,M2,M3,NQ,LP,NTA,ita(900),T(900)
     */BLQQ/QQ(2,900),DOM(900)
      DIMENSION r(2),DR(2),F(2),G(2,2)
      DO I=1,2
         F(I)=0.0
         DO J=1,2
            G(I,J)=0.0
         END DO
      END DO
      NTA=0
      do i=1,lp
         ita(i)=0
      end do
      OM=0.0
      DO I=1,LP
         QR=QQ(1,I)*DR(1)+QQ(2,I)*DR(2)
         QS=(QQ(1,I)*R(1)+QQ(2,I)*R(2))*1.00001
         IF((QS.LT.T(I)).OR.(QR.LT.DK))   GOTO 1
C     îOPMàPOBAHàE HAèPABãEHàâ AKTàBHOÉO MàKPOÑEîOPMàPOBAHàü
         NTA=NTA+1
         ITA(NTA)=I
         OM=OM+DOM(I)
         DO K1=1,2
            FQ=QQ(K1,I)*DOM(I)
            F(K1)=F(K1)+FQ
            DO K2=1,K1
               G(K1,K2)=G(K1,K2)+QQ(K2,I)*FQ
            END DO
         END DO
    1    CONTINUE
      END DO
      G(1,2)=G(2,1)
      RETURN
      END
C
      SUBROUTINE APDIN2M(R,DR,AK,ALMAX,  AL,F,G,OM)
C----------------------------------------------------------------------
C     OèPEÑEãEHàE F(I),G(I,J),OM èPà MOHOTOHHOM HAÉPìÜEHàà
C----------------------------------------------------------------------
      COMMON /BLTEOR/ITEOR
      DIMENSION   DR(2),F(2),G(2,2),R(2),ANU(2)
C
      DATA PI2/9.869804/
C
      CALL PDFT(R,DR,AK,almax, AL,F1)
      OM=0.0
      DO 5 I=1,2
      F(I)=0.0
      DO 5 J=1,2
    5 G(I,J)=0.0
C
      IF(AL.LT.almax)  RETURN
C
      C1=1.0-COS(AL)
      C2=1.0-COS(2.0*AL)
      C3=1.0-COS(3.0*AL)
      C4=1.0-COS(4.0*AL)
      C5=1.0-COS(5.0*AL)
      OM=(3.0*C1-C3/3.0)*PI2/2.0
      Fal=(C2-C4/4.0)*PI2/4.0
      Aal=(C1+C3/6.0-C5/10.0)*PI2/4.0
      Bal=(5.0*C1-C3/1.2+C5/10.0)*PI2/16.0
      IF(ITEOR.EQ.3) Bal=0.0
      DO I=1,2
         ANU(I)=R(I)/F1
         F(I)=Fal*ANU(I)
      END DO
      DO I=1,2
         DO J=1,2
            G(I,J)=(Aal-Bal)*ANU(I)*ANU(J)
            IF(I.EQ.J) G(I,J)=G(I,J)+Bal
         END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE PDFT(R,DR,AK,almax,  AL,F1)
C----------------------------------------------------------------------
C    OèPEÑEãEHàE èAPAMETPA AL=ALFA
C  (R-AKTàBHõE HAèPüÜEHàü, AK-KAèèA)
C----------------------------------------------------------------------
      DIMENSION R(2),DR(2)
      COMMON
     */MDF/T0,ES,AA,AN,RR0,AK0,ETA,AN1
C
      AL=0.0
      F1=SQRT(R(1)*R(1)+R(2)*R(2))
      A1=ES+AK
C
      IF(A1.GE.F1) RETURN
C
      AL=ACOS(A1/F1)
      if(al.lt.almax) return
      IF((R(1)*DR(1)+R(2)*DR(2)).LT.0.0) AL=0.0
      RETURN
      END
C
      SUBROUTINE APDDDD(II,KK)
C----------------------------------------------------------------------
C    èOCTPOEHàE MATPàñ ÜECTKOCTà TEOPàà MàKPOÑEîOPMAñàà
C     II-HOMEP ùãEMEHTA, KK- HOMEP ìáãA àHTEÉPàPOBAHàü
C----------------------------------------------------------------------
      COMMON
     */BLS/SUM(5500),FP(5500)
     */BL4/X(3),C(3),NNOD
     */BL2A/NEL,NC(8,110)
     */BL9/DP(2,2),DD(2,2)
     */BLITA/M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLTCOM/TCOM(400000)
     */BLTEOR/ITEOR
      DIMENSION R(2),DR(2)
C
      N1=NNOD+(II-1)*20+(KK-1)*5
           DO 10 J=1,2
           DR(J)=FP(N1+2+J)
  10       R(J)=SUM(N1+2+J)
           DK=FP(N1+5)
           AK=SUM(N1+5)
C
      IF(ITEOR.EQ.1) THEN
           N3=(II-1)*4+KK
           DO I=1,LP
               T(I)=TCOM(LP*(N3-1)+I)
           END DO
      END IF
C
      CALL APDDP(II,KK,R,DR,AK,DK)
      RETURN
      END
C
      SUBROUTINE APMATR(A)
C----------------------------------------------------------------------
C          èOCTPOEHàE MATPàñõ OÅPATHOâ K MATPàñE A(2,2)
C               PAáMEôAETCü HA MECTE MATPàñõ A(2,2)
C----------------------------------------------------------------------
      DIMENSION A(2,2)
C
      D=A(1,1)*A(2,2)-A(1,2)*A(2,1)
C
      D1=A(1,1)
      A(1,1)=A(2,2)/D
      A(1,2)=-A(1,2)/D
      A(2,1)=A(1,2)
      A(2,2)=D1/D
C
      RETURN
      END
C
      SUBROUTINE APST2(DR,DK)
C-----------------------------------------------------------------------
C     èOÑèPOÉPAMMA Ñãü BõóàCãEHàü àHTEHCàBHOCTà PAáPEòAûôàX
C              HAèPüÜEHàâ èOCãE ìÉãOBOâ TOóKà.
C-----------------------------------------------------------------------
      COMMON
     */BLITA/ M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLQQ/ QQ(2,900),DOM(900)
     */MDF/ T0,ES,AA,AN,RR0,AK0,ETA,AN1
C
      DIMENSION DR(2)
      LPA=1
C
      DO 1 I=1,LP
            QDS=0.0
            DO 2 JJ=1,2
    2       QDS=QDS+QQ(JJ,I)*DR(JJ)
C   OÅãACTú OCTAãúHõX HAèPABãEHàâ
      T(I)=T(I)+DK
      IF(ITA(LPA).NE.I.OR.LPA.GT.NTA)   GOTO 1

C   OÅãACTú AKTàBHõX HAèPABãEHàâ
      T(I)=T(I)+QDS-DK
      LPA=LPA+1
C   OÅãACTú èPOTàBOèOãOÜHõX HAèPABãEHàâ
      LPR=LP-I+1
      T(LPR)=T(LPR)-ETA*QDS+ETA*DK
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE APDFF(N)
      COMMON
     */BL1A/LH,AA1(1100)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(500),RP2(500),RP(500)
     */BL15/JIN(8)
     */BLS/SUM(5500),FP(5500)
     */BLITA/ M1,M2,M3,NQ,LP,NTA,ITA(900),T(900)
     */BLTCOM/TCOM(400000)
     */BLTEOR/ITEOR
     */MDF/T0,ES,AA,AN,RR,AK0,ETA,AN1
     */BLALFA/ALFA(4,110),AL
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
C
      DIMENSION DUL(8),DS(2),DR(2),de1(2),UL(8),DE(2),R(2)
C
C      DO 10 I=1,N
C  10  FP(I)=0.0
      NUMITER=3
      LL=0
  20  LL=LL+1
      L2=0
      DO 30 I=1,NNOD
  30  AA1(I)=RP(I)
         CALL APDFOR(NNOD,1)
         CALL APDTRT(1,NNOD,NNOD,-1)
         CALL APDSUB(1,NNOD,0)
      DO 40 II=1,NEL
         CALL APDELM(II,IG)
      DO 50 I=1,8
      J=JIN(I)
      UL(I)=SUM(J)
  50  DUL(I)=AA1(J)
      KK=1
      DO 40 M=1,2
      DO 40 N2=1,2
         CALL APDCOB(X(M),X(N2),IG)
         CALL APDDEF(DUL,DE)
         CALL APDDDD (II,KK)
         CALL APDES(ii,kk,DE,  DS,DR,de1,DK)
      N1=NNOD+(II-1)*20+(KK-1)*5
      DO J=1,2
         FP(N1+J)=DS(J)
         FP(N1+2+J)=DR(J)
      END DO
      FP(N1+5)=DK
      IF(LL.NE.NUMITER) GO TO 45
      txz=sum(n1+1)
      tyz=sum(n1+2)
      dtxz=ds(1)
      dtyz=ds(2)
      dexz=de(1)
      deyz=de(2)
      www(kk,ii)=www(kk,ii)+((txz+dtxz/2.0)*dexz+(tyz+dtyz/2.0)*deyz)*2.
      alambd(kk,ii)=alambd(kk,ii)+sqrt(de1(1)*de1(1)+de1(2)*de1(2))
      IF(ITEOR.EQ.1) THEN
           CALL APST2(DR,DK)
           N3=(II-1)*4+KK
           DO I=1,LP
              TCOM(LP*(N3-1)+I)=T(I)
           END DO
      END IF
      IF(ITEOR.EQ.2.OR.ITEOR.EQ.3) THEN
         R(1)=SUM(N1+3)+DR(1)
         R(2)=SUM(N1+4)+DR(2)
         AK=SUM(N1+5)+DK
         F1=SQRT(R(1)*R(1)+R(2)*R(2))
         A1=ES+AK
         AL=0
         IF(A1.LT.F1) AL=ACOS(A1/F1)
         IF(AL.GT.ALFA(KK,II)) ALFA(KK,II)=AL
      END IF
   45 KK=KK+1
      L2=L2+NTA
   40 CONTINUE
      IF(LL.LT.numiter) GOTO 20
      DO 100 I=1,NNOD
  100 FP(I)=AA1(I)
      IF(ITEOR.EQ.1) PRINT 110,L2
  110 FORMAT(20X,'óàCãO AKTàBHõX HAèPABãEHàâ PABHO  ',I6/)
C
      RETURN
      END
c
      SUBROUTINE APDRES(iii,IPR,IU,ISXZ,ISYZ,ISI,IBUIL,NEE,NU)
C----------------------------------------------------------------------
C          èPOñEÑìPA BõBOÑA PEáìãúTATOB PEòEHàü áAÑAóà
C----------------------------------------------------------------------
      COMMON
     */BL1/A(60000)
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL7/XY(8,2),AKT(8,8)
     */BL8/B(2,8),DET,R,Z,yak(2,2)
     */BL14/EPS(2,4),S(2,4),UL(8)
     */BL15/JIN(8)
     */BLS/SUM(5500),FP(5500)
     */BLTEOR/ ITEOR
     */MDF/ T0,ES,AA,AN,RR0,AK0,ETA,AN1
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
     */blj/zjsr,zjupr
C
      DIMENSION E(2),U(550),SXZ(550),SYZ(550),SI(550),rrr(4),zzz(4),
     *          alamb(550),tau(2,4),xksi(4),xeta(4),yksi(4),yeta(4)
C
      OPEN(71,FILE='U.DAT')
      WRITE(71,*) '"X","Y","U"'
      DO I=1,NNOD
         U(I)=SUM(I)
         WRITE(71,*) XYT(1,I),',',XYT(2,I),',',U(I)
      END DO
      CLOSE(71)
C
      do i=1,2*jnkont
         zjint(i)=0.0
      end do
      jelem=1
c      print *,jnkont,jnkel
c      print *,(jkel(j),j=1,jnkel)
C      print *,(jks(j),j=1,jnkel)
C      print *,(jkf(j),j=1,jnkel)
C
      OPEN(71,FILE='SI.DAT')
      WRITE(71,*) '"X","Y","SI"'
      OPEN(72,FILE='SXZ.DAT')
      WRITE(72,*) '"X","Y","SXZ"'
      OPEN(73,FILE='SYZ.DAT')
      WRITE(73,*) '"X","Y","SYZ"'
      OPEN(74,FILE='LAMDA.DAT')
      WRITE(74,*) '"X","Y","LAMBDA"'
      DO II=1,NEL
         CALL APDELM(II,IG)
         DO L=1,8
            J=JIN(L)
            UL(L)=SUM(J)
         END DO
         KK=1
         DO N1=1,2
            DO N2=1,2
               CALL APDCOB(X(N1),X(N2),IG)
               rrr(2*(n1-1)+n2)=r
               zzz(2*(n1-1)+n2)=z
               xksi(kk)=yak(1,1)
               xeta(kk)=yak(2,1)
               yksi(kk)=yak(1,2)
               yeta(kk)=yak(2,2)
               CALL APDDEF(UL,E)
               DO JJ=1,2
                  EPS(JJ,KK)=E(JJ)
               END DO
               KK=KK+1
            END DO
         END DO
         DO KK=1,4
            N1=NNOD+(II-1)*20+(KK-1)*5
            KN=4*(II-1)+KK
            SXZ(KN)=SUM(N1+1)
            SYZ(KN)=SUM(N1+2)
            tau(1,kk)=sxz(kn)
            tau(2,kk)=syz(kn)
            sii=sqrt(sxz(kn)*sxz(kn)+syz(kn)*syz(kn))
            write(71,*) rrr(kk),',',zzz(kk),',',sii/(sqrt(2.0)*t0)
            write(72,*) rrr(kk),',',zzz(kk),',',sxz(kn)/(sqrt(2.0)*t0)
            write(73,*) rrr(kk),',',zzz(kk),',',syz(kn)/(sqrt(2.0)*t0)
            write(74,*) rrr(kk),',',zzz(kk),',',alambd(kk,ii)/es
         END DO
         if(jelem.le.jnkel.and.ii.eq.jkel(jelem)) then
            jn=jks(jelem)
            w1=www(1,ii)
            w2=www(3,ii)
            exz1=eps(1,1)
            exz2=eps(1,3)
            txz1=tau(1,1)
            txz2=tau(1,3)
            tyz1=tau(2,1)
            tyz2=tau(2,3)
            xksi1=xksi(1)
            xksi2=xksi(3)
            yksi1=yksi(1)
            yksi2=yksi(3)
            zjint(jn)=zjint(jn)+
     +                ((w1-txz1*exz1*2.0)*yksi1+tyz1*exz1*2.0*xksi1)+
     +                ((w2-txz2*exz2*2.0)*yksi2+tyz2*exz2*2.0*xksi2)
            w1=www(2,ii)
            w2=www(4,ii)
            exz1=eps(1,2)
            exz2=eps(1,4)
            txz1=tau(1,2)
            txz2=tau(1,4)
            tyz1=tau(2,2)
            tyz2=tau(2,4)
            xksi1=xksi(2)
            xksi2=xksi(4)
            yksi1=yksi(2)
            yksi2=yksi(4)
            zjint(jn+jnkont)=zjint(jn+jnkont)+
     +                ((w1-txz1*exz1*2.0)*yksi1+tyz1*exz1*2.0*xksi1)+
     +                ((w2-txz2*exz2*2.0)*yksi2+tyz2*exz2*2.0*xksi2)
            jelem=jelem+1
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.nee) then
            write(31,*) sqrt(eps(2,nu)*eps(2,nu)+eps(1,nu)*eps(1,nu)),
     *	              sqrt(tau(2,nu)*tau(2,nu)+tau(1,nu)*tau(1,nu))
         end if
c         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.nee) then
c            write(32,*) eps(1,nu)/es,eps(2,nu)/es
c         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.25) then
            write(32,*) eps(1,1)/es,eps(2,1)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.25) then
            write(322,*) eps(1,3)/es,eps(2,3)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.30) then
            write(323,*) eps(1,1)/es,eps(2,1)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.30) then
            write(324,*) eps(1,3)/es,eps(2,3)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.49) then
            write(325,*) eps(1,1)/es,eps(2,1)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.49) then
            write(326,*) eps(1,3)/es,eps(2,3)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.54) then
            write(327,*) eps(1,1)/es,eps(2,1)/es
         end if
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.54) then
            write(328,*) eps(1,3)/es,eps(2,3)/es
         end if
c
         if(ipr.eq.0.and.ibuil.eq.1.and.ii.eq.nee) then
            t02=sqrt(2.0)*t0
            write(33,*) tau(1,nu)/t02,tau(2,nu)/t02
         end if
      END DO
      CLOSE(71)
      CLOSE(72)
      CLOSE(73)
      CLOSE(74)
      if(jnkel.ne.0) then
         zjsr=0.0
         do i=1,jnkont-1
            zjsr=zjsr+zjint(i)
            zjsr=zjsr+zjint(i+jnkont)
         end do
         zjsr=zjsr/(2*(jnkont-1))
         IF(ITEOR.NE.0.and.iii.eq.1) THEN
            do i=1,2*jnkont
               zjint(i)=zjint(i)/zjupr
            end do
            zjsr=zjsr/zjupr
         END IF
         do i=1,jnkont
            print *,'kontur=',2*i-1,'   Jint=',zjint(i),'  ',(zjint(i)-
     *             zjsr)/zjsr*100,' %'
            print *,'kontur=',2*i,'   Jint=',zjint(i+jnkont),' ',
     *            (zjint(i+jnkont)-zjsr)/zjsr*100,' %'
         end do
         print *,'  Jstab=',zjsr
         if(ipr.eq.1) then
                  OPEN(85,FILE='JINT.DAT')
                  do i=1,jnkont
                     write(85,*) 2*i-1,zjint(i),zjsr
                     write(85,*) 2*i,zjint(i+jnkont),zjsr
                  end do
                  CLOSE(85)
         end if
      end if
      do i=1,4
         do j=1,nel
           alamb(i+4*(j-1))=alambd(i,j)/es
         end do
      end do
      call apdnod(alamb,nnod)
      CALL APDNOD(SXZ,NNOD)
      CALL APDNOD(SYZ,NNOD)
      SIMAX=0.0
      DO I=1,NNOD
         SI(I)=SQRT(SXZ(I)*SXZ(I)+SYZ(I)*SYZ(I))
         IF(SI(I).GT.SIMAX) SIMAX=SI(I)
      END DO
      PRINT *,'  SIMAX=',SIMAX/(sqrt(2.0)*t0)
      DO I=1,NNOD
         SI(I)=SI(I)/(sqrt(2.0)*T0)
      END DO
      if((iii.eq.0.and.iteor.ne.0).or.(iii.eq.1.and.iteor.eq.0))
     * go to 500
      IF(IPR.EQ.0) GOTO 401
      WRITE (6,200)
  200 FORMAT(//31X,'PEáìãúTATõ  PACóETA')
      WRITE (6,220)
  220 FORMAT(/1X,80('-')/1X,'I HOM.I KOOPÑàH. I KOOPÑàH. I',
     *' HAèPüÜEHàEI HAèPüÜEHàEI  èEPEMEô.  I àHTEH. HAèP.I'/
     *1X,'I ìá. I  ìáãA  X I  ìáãA  Y I    SXZ    I    SYZ    I',
     *'  ìáãA  U   I     SI      I'/
     *1X,80('-'))
      DO 250 I=1,NNOD
  250 WRITE (6,230)I,XYT(1,I),XYT(2,I),SXZ(I),SYZ(I),U(I),SI(I)
  230 FORMAT(1X,'I ',I3,' I  ',F7.2,' I  ',F7.2,' I ',E9.3,
     *' I ',E9.3,' I ',E10.3,' I',F12.4,' I')
      WRITE (6,231)
  231 FORMAT(1X,80('-'))
      do i=1,jnkont
         write(6,*) 'kontur=',i,'  ','Jint=',zjint(i)
      end do
  401 OPEN(7,FILE='XY.SAM')
      DO 402 JI=1,NNOD
  402 WRITE (7,403) XYT(1,JI),XYT(2,JI)
  403 FORMAT(5X,2F10.3)
      CLOSE(7)
      OPEN(7,FILE='NCM.SAM')
      JI1=1
      DO 404 JI=1,NEL
  404 WRITE (7,405) (NC(I,JI),I=1,8),JI1
  405 FORMAT(5X,9I6)
      CLOSE(7)
      OPEN(7,FILE='INFO.SAM')
      WRITE(7,406)
  406 FORMAT(' kel  nnod  mz  np  nf  lh  istr ips nsym space'
     *       '   xsym      ysym')
      WRITE(7,407) NEL,NNOD
  407 FORMAT(2I5,'   0   0   0  318   2   0   0   .1    .0000   .0000'
     */'  mzt net not'/'   0   0   0'/'    0 = iseb')
      CLOSE(7)
C
      ila=1
      IF(IU.NE.1) GOTO 410
      OPEN(7,FILE='U.SAM')
      WRITE(7,501) (U(I),I=1,NNOD)
      CLOSE(7)
  410 IF(ISXZ.NE.1) GOTO 420
      OPEN(7,FILE='SXZ.SAM')
      WRITE(7,501) (SXZ(I),I=1,NNOD)
      CLOSE(7)
  420 IF(ISYZ.NE.1) GOTO 430
      OPEN(7,FILE='SYZ.SAM')
      WRITE(7,501) (SYZ(I),I=1,NNOD)
      CLOSE(7)
  430 IF(ISI.NE.1) GOTO 440
      OPEN(7,FILE='SI.SAM')
      WRITE(7,501) (SI(I),I=1,NNOD)
      CLOSE(7)
  440 IF(ILA.NE.1) GOTO 460
      OPEN(7,FILE='LAM.SAM')
      WRITE(7,501) (alamb(I),I=1,NNOD)
      CLOSE(7)
  460 CONTINUE
  501 FORMAT(5X,F10.4)
C
  500 RETURN
      END
      SUBROUTINE APDUPR(A0,K2)
C----------------------------------------------------------------------
C     èOCTPOEHàE ìèPìÉOÉO PEòEHàü áAÑAóà à MACòTAÅàPOBAHàE
C         PEáìãúTATA C ìóETOM ìCãOBàü èãACTàóHOCTà
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(500),RP2(500),RP(500)
     */BL14/ESP(2,4),SIG(2,4),UL(8)
     */BL15/JIN(8)
     */BLS/SUM(5500),FP(5500)
     */MDF/T0,ES,ADA,AN,RR0,AK0,ETA,AN1
     */BLTEOR/ ITEOR
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */blj/zjsr,zjupr
C
      DIMENSION E(2),A0(3,10),s(2),r(2)
C
      DO 10 I=1,NNOD
   10 A1(I)=0.0
C
C     èPABAü óACTú
C
      DO 20 I=1,NNOD
   20 A1(I)=RP1(I)*A0(1,K2)+RP2(I)*A0(2,K2)
C
      CALL APDFOR(NNOD,1)
      CALL APDTRT(1,NNOD,NNOD,-1)
      CALL APDSUB(1,NNOD,0)
C
C     èOCTPOEHàE èOãEâ èEPEMEôEHàâ, HAèPüÜEHàâ à ÑEîOPMAñàâ
C
      DO 30 J=1,NNOD
   30 SUM(J)=A1(J)
C
      DO 40 II=1,NEL
C
      CALL APDELM(II,IG)
C
      DO 50 I=1,8
      J=NC(I,II)
   50 UL(I)=A1(J)
C
      KK=1
      DO 40 N1=1,2
      DO 40 N2=1,2
      CALL APDCOB(X(N1),X(N2),IG)
      CALL APDDEF(UL,E)
C
      DO 60 JJ=1,2
   60 ESP(JJ,KK)=E(JJ)
C
      CALL APDDDD(II,KK)
      CALL APDSTR(KK)
C
      DO 70 J=1,2
      S(J)=SIG(J,KK)
   70 R(J)=S(J)*ES/T0
C
      www(kk,ii)=www(kk,ii)+(s(1)*r(1)+s(2)*r(2))
C
      M1=NNOD+(II-1)*20+(KK-1)*5
           SUM(M1+1)=S(1)
           SUM(M1+2)=S(2)
           SUM(M1+3)=R(1)
           SUM(M1+4)=R(2)
           SUM(M1+5)=0.0
   82         CONTINUE
      KK=KK+1
   40 CONTINUE
C
      CALL APDRES(0,1,1,1,1,1,0,0,0)
          ZJUPR=ZJSR
C
      RETURN
      END
