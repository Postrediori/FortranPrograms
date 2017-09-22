      PROGRAM PD
C----------------------------------------------------------------------
C             èPOÉPAMMA PEòEHàü ÉPAHàóHõX áAÑAó TEOPàà
C    MàKPOÑEîOPMAñàà èPà èãOCKOâ ÑEîOPMAñàà àãà èãOCKOM HAèPüÜEHHOM
C                COCTOüHàà à èPà CãOÜHOM HAÉPìÜEHàà
C----------------------------------------------------------------------
      COMMON
     */BL1A/LH,A1(1100)
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL6/RP1(1100),RP2(1100),RP(1100)
     */BL9/DP(4,4),DD(3,3),INS
     */MDF/T0,ES,AA,AN,RR,AK0,ETA,AN1,ANU
     */BLS/SUM(10000),FP(10000)
     */BLTEOR/ITEOR
     */blita/   m1,m2,m3,nq,lp,nta,ita(3375),t(3375)
     */bltcom/ tcom(1500000)
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
C
      DIMENSION A0(3,10),LA(10),DVS(3)
C
C    TPAEKTOPàü BHEòHEÉO HAÉPìÜEHàü áAÑAETCü OTHECEHHOâ
C     K àHTEHCàBHOCTà HAÉPìáKà, ÑOCTàÉAEMOâ HA èPEÑEãE
C            TEKìóECTà B OÑHOâ àá TOóEK OÅãACTà
C
      open (31,file='trajec')
      OPEN(8,FILE='$PDFIZ')
      read (8,*)
      READ (8,*) AA,AN,RR,AK0,ETA,T0,ES,AN1,ANU
      READ (8,20) K2,LTWO
  20  FORMAT(/2I5/)
      READ (8,30) ((A0(I,J),I=1,3),LA(J),J=1,K2)
  30  FORMAT(3F10.5,I5)
      CLOSE(8)
      do i=1,3
         do j=1,k2
            a0(i,j)=a0(i,j)*t0
         end do
      end do
C
      OPEN(6,FILE='RESULT')
      WRITE (6,40) RR,AN1,AK0,AN,AA,ETA,T0,ES,ANU
  40  FORMAT(//33X,'èAPAMETPõ MOÑEãà'
     *   /15X,'RR=',F7.4,'   AN1=',F7.3,'   AK0=',F7.3,'   AN=',F7.3
     *   /15X,'AA=',F7.3,'   ETA=',F7.3,'   T0 =',F7.1,'   ES=',F7.6,
     *   /35X,'ANU=',F5.2)
C
      CALL PDDATA(LTWO)
      CALL PDRESU(ins,ip,isx,isy,isxy,isi,iux,iuy)
C
      IF(ITEOR.EQ.0) WRITE (6,105)
      IF(ITEOR.EQ.1) WRITE (6,101)
      IF(ITEOR.EQ.2) WRITE (6,102)
      IF(ITEOR.EQ.3) WRITE (6,103)
  101 FORMAT(/25X,'PACóET èO TEOPàà MàKPOÑEîOPMAñàà')
  102 FORMAT(/25X,'PACóET èO ÑEîOPMAñàOHHOâ TEOPàà')
  103 FORMAT(/29X,'PACóET èO TEOPàà TEóEHàü')
  105 FORMAT(/34X,'ìèêìÉéÖ êÖòÖçàÖ')
C
      LH=LH*2
      NN=NNOD*2
C
      WRITE (6,100) NEL,NNOD,LH
  100 FORMAT(//31X,'OÅôàE  èAPAMETPõ MKù'/
     *19X,'KOãàóECTBO èãOCKàX ùãEMEHTOB',7X,'NEL =',I3/
     *19X,'KOãàóECTBO  ìáãOB  Kù ',13X,'NNOD=',I4/
     *19X,'òàPàHA ãEHTõ MATPàñõ',15X,'LH  =',I4//)
C
      do i=1,nel
         do j=1,4
            www(j,i)=0.0
         end do
      end do
C
  401 OPEN(7,FILE='XY.SAM')
      DO JI=1,NNOD
         WRITE (7,403) XYT(1,JI),XYT(2,JI)
      END DO
  403 FORMAT(5X,2F10.3)
      CLOSE(7)
      OPEN(7,FILE='NCM.SAM')
      JI1=1
      DO JI=1,NEL
         WRITE (7,405) (NC(I,JI),I=1,8),JI1
      END DO
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
      NSUM=NN+32*NEL
      DO I=1,NSUM
        FP(I)=0.0
        SUM(I)=0.0
      END DO
C
      CALL PDQQ
C
      LPZ=LP*NEL*4
      DO I=1,LPZ
         TCOM(I)=ES
      END DO
C
      CALL PDUPR(ip,isx,isy,isxy,isi,iux,iuy,A0,K2)
C
      IF (ITEOR.EQ.0) GO TO 145
C
      WRITE(6,50)
   50 FORMAT(/25X,'áAÑAHA TPAEKTOPàü HAÉPìÜEHàü'
     *     //22X,'   P1        P2      K-BO òAÉOB')
      WRITE (6,60) ((A0(I,J),I=1,2),LA(J),J=1,K2)
   60 FORMAT(/18X,2F10.3,I10)
C
      do i=1,nel
         do j=1,4
            www(j,i)=0.0
            alambd(j,i)=0.0
         end do
      end do
C
      NSUM=NN+32*NEL
      DO I=1,NSUM
        FP(I)=0.0
        SUM(I)=0.0
      END DO
C
      K1=1
      KP=LA(K2)-1
C
      DO LOAD=1,KP
           PRINT 120,LOAD
  120      FORMAT(30X,'òAÉ',I3)
                CALL PDSE(LOAD,DVS,LA,K1,A0)
           DO I=1,NN
             RP(I)=RP1(I)*DVS(1)+RP2(I)*DVS(2)
           END DO
                CALL PDSOL(NSUM)
                CALL PDRES(0,0,0,0,0,0,0,0)
      END DO
C
      CALL PDRES(1,ip,isx,isy,isxy,isi,iux,iuy)
C
  145 CONTINUE
      close(31)
      CLOSE(9,STATUS='DELETE')
      STOP
      END PROGRAM PD
C*****
      SUBROUTINE PDUPR(ip,isx,isy,isxy,isi,iux,iuy,A0,K2)
C----------------------------------------------------------------------
C     èOCTPOEHàE ìèPìÉOÉO PEòEHàü áAÑAóà Ñãü HAÉPìáKà B KOHEóHOâ
C       TOóKE TPAEKTOPàà BHEòHEÉO HAÉPìÜEHàü à MACòTAÅàPOBAHàE
C             PEáìãúTATA C ìóETOM ìCãOBàü èãACTàóHOCTà
C----------------------------------------------------------------------
      COMMON
     */BL4/X(3),C(3),NNOD
     */BL1A/LH,A1(1100)
     */BL2A/NEL,NC(8,110)
     */BL6/RP1(1100),RP2(1100),RP(1100)
     */BL7/XY(8,2),AK(16,16)
     */BL14/EPS(4,4),SIG(4,4),UL(16)
     */BL15/LEC(16),JIN(8)
     */BLS/SUM(10000),FP(10000)
     */MDF/T0,ES,AA,AN,RR,AK0,EPSS,AN1,ANU
     */BLTEOR/ ITEOR
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
C
      DIMENSION A0(3,10),E(4),S(4),R(3)
C
      NN=NNOD*2
C
      DO I=1,NN
        A1(I)=0.0
      END DO
C
C     èPABAü óACTú
C
      DO I=1,NN
        A1(I)=RP1(I)*A0(1,K2)+RP2(I)*A0(2,K2)
      END DO
C
      CALL PDFOR(NN,0)
      CALL PDSTRT(1,NN,NN,-1)
      CALL PDBSUB(1,NN,0)
C
C     èOCTPOEHàE èOãEâ èEPEMEôEHàâ,ÑEîOPMAñàâ,HAèPüÜEHàâ
C
      DO I=1,NN
        SUM(I)=A1(I)
      END DO
C
      DO II=1,NEL
c
        CALL PDELM(II,IG)
C
        DO I=1,16
          J=LEC(I)
          UL(I)=A1(J)
        END DO
C
        KK=1
        DO M=1,2
          DO N=1,2
C
            CALL PDCOB(X(M),X(N),IG)
            CALL PDDEF(UL,E)
C
            DO I=1,4
              EPS(I,KK)=E(I)
            END DO
C
            CALL PDDDD(II,KK)
            CALL PDSTR(KK)
C
            DO J=1,4
              S(J)=SIG(J,KK)
            END DO
            CALL PDTEVE(E,R)
C
            nwww=2*(m-1)+n
            www(nwww,ii)=www(nwww,ii)+
     +                   (s(1)*e(1)+s(2)*e(2)+2.0*s(3)*e(3))/2.0
C
            M1=NN+(II-1)*32+(KK-1)*8
            DO J=1,4
              SUM(M1+J)=S(J)
            END DO
            DO J=1,3
              SUM(M1+4+J)=R(J)
            END DO
            SUM(M1+8)=0.0
            KK=KK+1
          END DO
        END DO
      END DO
      CALL PDRES(0,ip,isx,isy,isxy,isi,iux,iuy)
      RETURN
      END SUBROUTINE PDUPR
C*****
      SUBROUTINE PDRES(iii,IP,ISX,ISY,ISXY,ISI,IUX,IUY)
C----------------------------------------------------------------------
C        èPOñEÑìPA BõBOÑA PEáìãúTATOB PACóETA HA èEóATú
C----------------------------------------------------------------------
      COMMON
     */BL1/A(400000)
     */BL2/XYT(2,550)
     */BL2A/NEL,NC(8,110)
     */BL4/X(3),C(3),NNOD
     */BL8/B(4,16),DET,R,Z,AN(8),yak(2,2)
     */BL14/EPS(4,4),SIG(4,4),UL(16)
     */BL15/LEC(16),JIN(8)
     */BLS/SUM(10000),FP(10000)
     */BLTEOR/ ITEOR
     */MDF/T0,ES,AA,ANN,RR,AK0,EPSS,AN1,ANU
     */bljint/jnkont,jnkel,jkel(50),jks(50),jkf(50),www(4,110),zjint(10)
     */bllamb/alambd(4,110)
C
      DIMENSION UX(550),UY(550),SX(550),SY(550),SZ(550),SXY(550),
     *          E(4),SG1(550),SG2(550),SI(550),alamb(550),rrr(4),zzz(4),
     *          sigma(3,4),xksi(4),xeta(4),yksi(4),yeta(4),dvx(4),
     *          dvy(4),zjinty(10),zji(10)
C
      OPEN(76,FILE='SURFER.DAT')
      WRITE(76,*) '"X","Y","SI"'
      OPEN(77,FILE='LAMBDA.DAT')
      WRITE(77,*) '"X","Y","LAMBDA"'
      OPEN(71,FILE='SX.DAT')
      WRITE(71,*) '"X","Y","SX"'
      OPEN(72,FILE='SY.DAT')
      WRITE(72,*) '"X","Y","SY"'
      OPEN(73,FILE='SXY.DAT')
      WRITE(73,*) '"X","Y","SXY"'
      OPEN(74,FILE='SI.DAT')
      WRITE(74,*) '"X","Y","SI"'
      NN=NNOD*2
      DO I=1,NNOD
         UX(I)=SUM(2*I-1)
         UY(I)=SUM(2*I)
      END DO
      do i=1,2*jnkont
         zjint(i)=0.0
         zjinty(i)=0.0
      end do
      jelem=1
      DO II=1,NEL
           CALL PDELM(II,IG)
           DO I=1,16
              J=LEC(I)
              UL(I)=SUM(J)
           END DO
           KK=1
           DO M=1,2
              DO N=1,2
                 CALL PDCOB(X(M),X(N),IG)
                 rrr(2*(m-1)+n)=r
                 zzz(2*(m-1)+n)=z
                 xksi(kk)=yak(1,1)
                 xeta(kk)=yak(2,1)
                 yksi(kk)=yak(1,2)
                 yeta(kk)=yak(2,2)
                 CALL PDDEF(UL,E)
                 DO I=1,4
                    EPS(I,KK)=E(I)
                 END DO
                 dvx(kk)=0.0
                 do i=2,16,2
                    dvx(kk)=dvx(kk)+b(1,i-1)*ul(i)
                 end do
                 dvy(kk)=0.0
                 do i=1,16,2
                    dvy(kk)=dvy(kk)+b(2,i+1)*ul(i)
                 end do
                 KK=KK+1
              END DO
           END DO
           DO KK=1,4
              N1=NN+(II-1)*32+(KK-1)*8
              KN=4*(II-1)+KK
              SX(KN)=SUM(N1+1)
              SY(KN)=SUM(N1+2)
              SXY(KN)=SUM(N1+3)
              SZ(KN)=SUM(N1+4)
              sigma(1,kk)=sx(kn)
              sigma(2,kk)=sy(kn)
              sigma(3,kk)=sxy(kn)
              SG=SQRT(SXY(kn)*SXY(kn)+(SX(kn)-SY(kn))**2/4.0)
              S1=(SX(kn)+SY(kn))/2.0+SG
              S2=(SX(kn)+SY(kn))/2.0-SG
              sii=SQRT((S1-S2)**2+S1**2+S2**2)/1.414
              write(76,*) rrr(kk),',',zzz(kk),',',sii/(sqrt(2.0)*t0)
              write(77,*) rrr(kk),',',zzz(kk),',',alambd(kk,ii)/es
              write(71,*) rrr(kk),',',zzz(kk),',',sx(kn)
              write(72,*) rrr(kk),',',zzz(kk),',',sy(kn)
              write(73,*) rrr(kk),',',zzz(kk),',',sxy(kn)
              write(74,*) rrr(kk),',',zzz(kk),',',sii
           END DO
           if(jelem.le.jnkel.and.ii.eq.jkel(jelem)) then
              jn=jks(jelem)
              w1=www(1,ii)
              w2=www(3,ii)
              dux1=eps(1,1)
              dux2=eps(1,3)
              duy1=eps(2,1)
              duy2=eps(2,3)
              dvx1=dvx(1)
              dvx2=dvx(3)
              dvy1=dvy(1)
              dvy2=dvy(3)
              sxx1=sigma(1,1)
              sxx2=sigma(1,3)
              syy1=sigma(2,1)
              syy2=sigma(2,3)
              sxy1=sigma(3,1)
              sxy2=sigma(3,3)
              xksi1=xksi(1)
              xksi2=xksi(3)
              yksi1=yksi(1)
              yksi2=yksi(3)
              zjint(jn)=zjint(jn)+
     +                  ( (w1-sxx1*dux1-sxy1*dvx1)*yksi1+
     +                    (sxy1*dux1+syy1*dvx1)*xksi1 )+
     +                    ( (w2-sxx2*dux2-sxy2*dvx2)*yksi2+
     +                    (sxy2*dux2+syy2*dvx2)*xksi2 )
              zjinty(jn)=zjinty(jn)+
     +                  ( (syy1*duy1+sxy1*dvy1-w1)*xksi1-
     -                    (sxy1*duy1+sxx1*dvy1)*yksi1 )+
     +                    ( (syy2*duy2+sxy2*dvy2-w2)*xksi2-
     -                    (sxy2*duy2+sxx2*dvy2)*yksi2 )
              zji(jn)=sqrt(zjint(jn)*zjint(jn)+
     +                  zjinty(jn)*zjinty(jn))
              w1=www(2,ii)
              w2=www(4,ii)
              dux1=eps(1,2)
              dux2=eps(1,4)
              duy1=eps(2,2)
              duy2=eps(2,4)
              dvx1=dvx(2)
              dvx2=dvx(4)
              dvy1=dvy(2)
              dvy2=dvy(4)
              sxx1=sigma(1,2)
              sxx2=sigma(1,4)
              syy1=sigma(2,2)
              syy2=sigma(2,4)
              sxy1=sigma(3,2)
              sxy2=sigma(3,4)
              xksi1=xksi(2)
              xksi2=xksi(4)
              yksi1=yksi(2)
              yksi2=yksi(4)
              zjint(jn+jnkont)=zjint(jn+jnkont)+
     +                 ( (w1-sxx1*dux1-sxy1*dvx1)*yksi1+
     +                   (sxy1*dux1+syy1*dvx1)*xksi1 )+
     +                   ( (w2-sxx2*dux2-sxy2*dvx2)*yksi2+
     +                   (sxy2*dux2+syy2*dvx2)*xksi2 )
              zjinty(jn+jnkont)=zjinty(jn+jnkont)+
     +                  ( (syy1*duy1+sxy1*dvy1-w1)*xksi1-
     -                    (sxy1*duy1+sxx1*dvy1)*yksi1 )+
     +                    ( (syy2*duy2+sxy2*dvy2-w2)*xksi2-
     -                    (sxy2*duy2+sxx2*dvy2)*yksi2 )
              zji(jn+jnkont)=sqrt(zjint(jn+jnkont)*zjint(jn+jnkont)+
     +			             zjinty(jn+jnkont)*zjinty(jn+jnkont))
      		jelem=jelem+1
         end if
         if(ip.eq.0.and.ii.eq.10) then
             write(31,*) eps(2,1),sigma(2,1)
         end if
      END DO
      close(76)
      close(77)
      close(71)
      close(72)
      close(73)
      close(74)
      if(jnkel.ne.0) then
         zjsr=0.0
         do i=1,jnkont-1
            zjsr=zjsr+zji(i)
            zjsr=zjsr+zji(i+jnkont)
         end do
         zjsr=zjsr/(2*(jnkont-1))
         do i=1,jnkont
            print *,'kontur=',2*i-1,'   Jint=',zji(i),'  ',(zji(i)-
     *             zjsr)/zjsr*100,' %'  	     
            print *,'kontur=',2*i,'   Jint=',zji	(i+jnkont),' ',
     *            (zji(i+jnkont)-zjsr)/zjsr*100,' %'  	     
         end do
         print *,''  
         print *,'                     Jstab=',zjsr
      end if
c      do i=1,jnkont
c         print *,'kontur=',2*i-1,'   Jint=',zjint(i)
c         print *,'kontur=',2*i,'   Jint=',zjint(i+jnkont)
c      end do
      do i=1,4
         do j=1,nel
           alamb(i+4*(j-1))=alambd(i,j)/es
         end do
      end do
      call pdnody(alamb,nnod)
      CALL PDNODY(SX,NNOD)
      CALL PDNODY(SY,NNOD)
      CALL PDNODY(SXY,NNOD)
      CALL PDNODY(SZ,NNOD)
      SIMAX=0.0
      DO I=1,NNOD
         SG=SQRT(SXY(I)*SXY(I)+(SX(I)-SY(I))**2/4.0)
         SG1(I)=(SX(I)+SY(I))/2.0+SG
         SG2(I)=(SX(I)+SY(I))/2.0-SG
         SI(I)=SQRT((SG1(I)-SG2(I))**2+SG1(I)**2+SG2(I)**2)/1.414
         IF(SIMAX.LT.SI(I)) SIMAX=SI(I)
      END DO
      PRINT *,'                     SImax=',SIMAX/(sqrt(2.0)*t0)
      print *,''
      DO I=1,NNOD
         SI(I)=SI(I)/(sqrt(2.0)*T0)
      END DO
C
      if((iii.eq.0.and.iteor.ne.0).or.(iii.eq.1.and.iteor.eq.0))
     * go to 500
      IF(IP.EQ.0) GOTO 400
      WRITE(6,200)
  200 FORMAT(//30X,'PEáìãúTATõ  PACóETA')
      WRITE(6,220)
  220 FORMAT(/1X,104('-')/1X,'I HOM.I KOOPÑàH. I KOOPÑàH. I',
     *' HAèPüÜEHàEI HAèPüÜEHàEI HAèPüÜEHàEI àçíÖçëàÇH.I',
     *'  èEPEMEô.  I  èEPEMEô.  I'/
     *1X,'I ìá. I  ìáãA  X I',
     *'  ìáãA  Y I    SX     I    SY     I    SXY    I',
     *'    SI     I  ìáãA  UX  I  ìáãA  UY  I'/
     *1X,104('-'))
      DO 250 I=1,NNOD
  250 WRITE(6,230) I,XYT(1,I),XYT(2,I),SX(I),SY(I),SXY(I),SI(I),
     *UX(I),UY(I)
  230 FORMAT(1X,'I ',I3,' I  ',F7.2,' I  ',F7.2,' I ',E9.3,
     *' I ',E9.3,' I ',E9.3,' I ',E9.3,' I ',E10.3,
     *' I ',E10.3,' I')
      WRITE(6,231)
  231 FORMAT(1X,104('-'))
C
  400 IF(ISX.NE.1) GOTO 410
      OPEN(7,FILE='SX.SAM')
      WRITE(7,501) (SX(I),I=1,NNOD)
      CLOSE(7)
  410 IF(ISY.NE.1) GOTO 420
      OPEN(8,FILE='SY.SAM')
      WRITE(8,501) (SY(I),I=1,NNOD)
      CLOSE(8)
  420 IF(ISXY.NE.1) GOTO 430
      OPEN(9,FILE='SXY.SAM')
      WRITE(9,501) (SXY(I),I=1,NNOD)
      CLOSE(9)
  430 IF(ISI.NE.1) GOTO 440
      OPEN(10,FILE='SI.SAM')
      WRITE(10,501) (SI(I),I=1,NNOD)
      CLOSE(10)
  440 IF(IUX.NE.1) GOTO 450
      OPEN(11,FILE='UX.SAM')
      WRITE(11,502) (UX(I),I=1,NNOD)
      CLOSE(11)
  450 IF(IUY.NE.1) GOTO 460
      OPEN(12,FILE='UY.SAM')
      WRITE(12,502) (UY(I),I=1,NNOD)
      CLOSE(12)
  460 CONTINUE
  501 FORMAT(5X,F10.4)
  502 FORMAT(5X,E10.3)
C
  500 RETURN
      END SUBROUTINE PDRES
C*****
      SUBROUTINE XABND(K,NNOD,MM,NN,NUF,IM,IPROB)
C----------------------------------------------------------------------
C    èOÑèPOÉPAMMA îOPMàPOBAHàü ÉPAHàóHõX ìCãOBàâ èO èEPEMEôEHàüM
C----------------------------------------------------------------------
      COMMON
     */AL16/MZD(110),AMZ(110),NZ,NZ1
      DIMENSION MM(4),NN(5),NUF(12,30)
C
      IF(K.EQ.1.OR.K.EQ.2) GO TO 10
c      NZ=3
      IF(IM.EQ.0) NCENTR=(NNOD-1-2*MM(1)-1-7)/2+1
      IF(IM.EQ.1) NCENTR=(NNOD-1-2*MM(1)-1-9)/2+1
c      MZD(1)=2*NCENTR-1
c	MZD(2)=2*NCENTR
c      IF(IM.EQ.0) kol=4
c      IF(IM.EQ.1) kol=5
c	mzd(3)=2*(ncentr+2*kol+2*mm(3))
      IF(IM.EQ.0) kol=9
      if(im.eq.1) kol=11
      nz=kol+2*mm(3)+1
      MZD(1)=2*NCENTR-1
      do i=1,nz-1
         mzd(1+i)=2*(ncentr+i-1)
      end do
      IF(IPROB.EQ.3) THEN
         MZD(1)=1
         MZD(2)=2
      END IF
      GO TO 50
C
   10 IF(K.EQ.2) GO TO 20
      NZ=0
      DO L=1,3
        NZ=NZ+MM(L)
      END DO
      nz=2*nz+1
      DO L=1,NZ
        MZD(L)=2*L
      END DO
      NZ=NZ+1
      mzd(nz)=1
c	nz=23
c	mzd(1)=288
c	mzd(2)=290
c	mzd(3)=292
c	mzd(4)=294
c	mzd(5)=296
c	mzd(6)=298
c	mzd(7)=300
cc	mzd(8)=302
c	mzd(9)=304
c	mzd(10)=306
c	mzd(11)=308
c	mzd(12)=310
c	mzd(13)=312
c	mzd(14)=314
c	mzd(15)=316
c	mzd(16)=318
c	mzd(17)=320
c	mzd(18)=322
c	mzd(19)=324
c	mzd(20)=326
c	mzd(21)=328
cc	mzd(22)=330
c	mzd(23)=287
      GO TO 50
C
   20 NZ=0
      DO L=1,4
        NZ=NZ+2*NN(L)
      END DO
      NZ=NZ+2
      N=0
      DO L=1,4
        NNN=NN(L)*2+1
        NF=3*(L-1)+1
        DO M=2,NNN
          N=N+1
          m1=m
          if(iprob.ne.3.and.l.eq.3.and.m.eq.2) m1=1
          MZD(N)=2*NUF(NF,M1)-1-2
        END DO
      END DO
      mzd(n+1)=nnod*2-1
      mzd(n+2)=nnod*2
   50 RETURN
      END SUBROUTINE XABND
C*****
      subroutine pdresu(ins,ip,isx,isy,isxy,isi,iux,iuy)
      read(5,1) ins,ip,isx,isy,isxy,isi,iux,iuy
   1  format(///14x,i1////14x,i1/14x,i1/14x,i1/14x,i1/14x,i1
     */14x,i1/14x,i1)
      close(5)
      return
      end subroutine pdresu
