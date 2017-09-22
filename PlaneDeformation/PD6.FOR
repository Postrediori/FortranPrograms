      SUBROUTINE XADATA
C----------------------------------------------------------------------
C    OCHOBHAü èOÑèPOÉPAMMA PAáÅàBKà OÅãACTà HA Kù
C----------------------------------------------------------------------
      COMMON
     */AL1A/LH
     */AL2/XYT(2,550)
     */AL4/XX(3),C(3),ISTR,PREC,NNOD
     */AL7/XY(8,2),AKT(8,8)
     */AL16/MZD(110),AMZ(110),NZ,NZ1
     */AL21/NEP,NEP1,LTEL(90),LST(90), PP(2,90),PT(2,90)
     */AL2A/NC(8,110)
     */ajint/jnkont,jnkel,jkel(50),jks(50),jkf(50)
     */atop/itope(8),itopu
C
      DIMENSION MM(4),NN(5),X(4,12),Y(4,12),NSF(12),ISF(5,6),
     *          NLD(12,10),NUFD(12,30),NL(12,10),NUF(12,30),
     *          NLKL(12,10),NUFKL(12,30)
C
      OPEN(5,FILE='$VVOD')
      READ(5,20) ITEOR,IPROB,NVER1,NVER2,NGOR1,NGOR2,
     *           XZ1,XZ2,YZ1,YZ2,ZONE,ALFA,X0,Y0,R,IGRA,IM
   20 FORMAT(//////14X,I1////14X,I1//////14X,I2/14X,I2/14X,I2/14X,I2///
     *       14X,F6.2/14X,F6.2/14X,F6.2/14X,F6.2///14X,F6.2///14X,F6.2
     *       ///14X,F6.2/14X,F6.2/14X,F6.2//////14X,I1////14X,I1)
C
      IF (IPROB.EQ.4) CALL X4EL(iteor)
      IF (IPROB.EQ.4) GO TO 1006
      LN=4
      LM=3
      NF=LM*LN
C
      MM(1)=NVER1
      MM(2)=2
      IF(IPROB.EQ.3) MM(2)=4
      MM(3)=NVER2
C
      NN(1)=NGOR1
      NN(2)=1
      IF(IPROB.EQ.3) NN(2)=2
      NN(3)=1
      IF(IPROB.EQ.3) NN(3)=2
      NN(4)=NGOR2
C
      DO J=1,LN
        X(1,3*J-2)=0.0
        X(2,3*J-2)=XZ1
        X(3,3*J-2)=XZ1
        X(4,3*J-2)=0.0
C
        X(1,3*J-1)=XZ1
        X(2,3*J-1)=XZ1+ZONE
        X(3,3*J-1)=XZ1+ZONE
        X(4,3*J-1)=XZ1
C
        X(1,3*J)=XZ1+ZONE
        X(2,3*J)=XZ1+ZONE+XZ2
        X(3,3*J)=XZ1+ZONE+XZ2
        X(4,3*J)=XZ1+ZONE
      END DO
C
      DO I=1,LM
        Y(1,I)=0.0
        Y(2,I)=0.0
        Y(3,I)=YZ1
        Y(4,I)=YZ1
C
        Y(1,3+I)=YZ1
        Y(2,3+I)=YZ1
        Y(3,3+I)=YZ1+ZONE/2.
        Y(4,3+I)=YZ1+ZONE/2.
C
        Y(1,6+I)=YZ1+ZONE/2.
        Y(2,6+I)=YZ1+ZONE/2.
        Y(3,6+I)=YZ1+ZONE
        Y(4,6+I)=YZ1+ZONE
C
        Y(1,9+I)=YZ1+ZONE
        Y(2,9+I)=YZ1+ZONE
        Y(3,9+I)=YZ1+ZONE+YZ2
        Y(4,9+I)=YZ1+ZONE+YZ2
      END DO
C
      KLIN=0
      IF(IPROB.NE.2) ALFA=0.0
      IF(IPROB.NE.2) GO TO 21
      KLIN=1
      ALFA=ALFA*3.141592/180.0
           Y(3,4)=Y(3,4)-ZONE/2.*TAN(ALFA/2.)
           Y(4,4)=Y(4,4)-(ZONE/2.+XZ1)*TAN(ALFA/2.)
           Y(2,7)=Y(2,7)+ZONE/2.*TAN(ALFA/2.)
           Y(1,7)=Y(1,7)+(ZONE/2.+XZ1)*TAN(ALFA/2.)
   21 CONTINUE
C
      DO I=1,NF
        NSF(I)=1
      END DO
      NSF(5)=2
      IF(IPROB.EQ.3) NSF(5)=4
      NSF(8)=3
      IF(IPROB.EQ.3) NSF(8)=5
C
      LN2=LN+2
      LM2=LM+2
      DO J=1,LN2
        DO I=1,LM2
          ISF(I,J)=0
        END DO
      END DO
C
      LN1=LN+1
      LM1=LM+1
      DO J=2,LN1
        DO I=2,LM1
           ISF(J,I)=(I-1)+(J-2)*LM
        END DO
      END DO
C
      DO I=1,LM
        DO J=1,LN
          KOLVO=NN(J)
          NF1=(J-1)*LM+I
          DO K=1,KOLVO
            NLD(NF1,K)=MM(I)
          END DO
        END DO
      END DO
C
           NLD(5,1)=16
           NLD(8,1)=16
           IF(IM.EQ.1) NLD(5,1)=20
           IF(IM.EQ.1) NLD(8,1)=20
               IF(IPROB.NE.3) GO TO 22
               NLD(5,1)=28
               NLD(5,2)=8
               NLD(8,1)=8
               NLD(8,2)=28
               IF(IM.NE.0) GO TO 22
               NLD(5,1)=22
               NLD(5,2)=6
               NLD(8,1)=6
               NLD(8,2)=22
   22          CONTINUE
C
      DO I=1,LM
        DO J=1,LN
          KOLVO=2*NN(J)+1
          NF1=(J-1)*LM+I
          DO K=1,KOLVO,2
            NUFD(NF1,K)=2*MM(I)
          END DO
          DO L=2,KOLVO,2
            NUFD(NF1,L)=MM(I)
          END DO
        END DO
      END DO
C
           NUFD(5,1)=4
           NUFD(5,2)=34
           NUFD(5,3)=16
           NUFD(8,1)=16
           NUFD(8,2)=34
           NUFD(8,3)=4
           IF(IM.EQ.1) NUFD(5,2)=44
           IF(IM.EQ.1) NUFD(5,3)=20
           IF(IM.EQ.1) NUFD(8,1)=20
           IF(IM.EQ.1) NUFD(8,2)=44
                IF(IPROB.NE.3) GO TO 23
                NUFD(5,1)=8
                NUFD(5,2)=33
                NUFD(5,3)=31
                NUFD(5,4)=35
                NUFD(5,5)=17
                NUFD(8,1)=17
                NUFD(8,2)=35
                NUFD(8,3)=31
                NUFD(8,4)=33
                NUFD(8,5)=8
                IF(IM.NE.0) GO TO 23
                NUFD(5,2)=25
                NUFD(5,3)=25
                NUFD(5,4)=27
                NUFD(5,5)=13
                NUFD(8,1)=13
                NUFD(8,2)=27
                NUFD(8,3)=25
                NUFD(8,4)=25
   23           CONTINUE
C
      CALL XREG( MM,NN,LM,LN,NF,X,Y,ISF,NSF,NUFD,NLD,ZONE,
     *          KLIN,ALFA,X0,Y0,R,IM,   NUF,NL,NEL,NNOD)
C
      if(iprob.ne.3) then
         if(im.eq.0) then
            itope(1)=nl(5,1)+3
            itope(2)=nl(5,1)+7
            itope(3)=nl(5,1)+11
            itope(4)=nl(5,1)+15
            itope(5)=nl(8,1)+3
            itope(6)=nl(8,1)+7
            itope(7)=nl(8,1)+11
            itope(8)=nl(8,1)+15
            itopu=nuf(8,1)+8
                     else
            itope(1)=nl(5,1)+4
            itope(2)=nl(5,1)+9
            itope(3)=nl(5,1)+14
            itope(4)=nl(5,1)+19
            itope(5)=nl(8,1)+4
            itope(6)=nl(8,1)+9
            itope(7)=nl(8,1)+14
            itope(8)=nl(8,1)+19
            itopu=nuf(8,1)+10
         end if
      end if
c	print *,(itope(i),i=1,8),itopu
C
      jnkel=0 
      if(iprob.eq.1) then
         jnkel=16
         if(im.eq.1) jnkel=20
         jside=4
         do i=1,jnkel
            jkel(i)=nl(5,1)+i-1
            jkel(i+jnkel)=nl(8,1)+i-1
         end do
         jnkont=jnkel/4
         do i=1,jside
            do j=1,jnkont
                jkf((i-1)*jnkont+j)=1
                jkf((i-1)*jnkont+j+jnkel)=1
                jks((i-1)*jnkont+j)=j
                jks((i-1)*jnkont+j+jnkel)=j
            end do
         end do
         jnkel=2*jnkel
c         print *,jnkont,jnkel
c         print *,(jkel(j),j=1,jnkel)
c         print *,(jks(j),j=1,jnkel)
c         print *,(jkf(j),j=1,jnkel)
      end if
C
      NNN=NNOD
C
      IF(IPROB.eq.1.or.iprob.eq.2) then
        CALL XCRACK(NL,NUF,NEL,NNOD,MM,NN,IM)
C
        IF(IPROB.eq.2) then
          LNKL=LN/2
          LN2=LNKL+2
          LM2=LM+2
          DO J=1,LN2
            DO I=1,LM2
              ISF(I,J)=0
            END DO
          END DO
          LN1=LNKL+1
          LM1=LM+1
          DO J=2,LN1
            DO I=2,LM1
              ISF(J,I)=(I-1)+(J-2)*LM
            END DO
          END DO
          NFKL=LM*LNKL
C
          CALL XREG( MM,NN,LM,LNKL,NFKL,X,Y,ISF,NSF,NUFD,NLD,ZONE,
     *               KLIN,ALFA,X0,Y0,R,IM,   NUFKL,NLKL,NELKL,NNODKL)
        end if
      end if
C
      LH=0
      DO I=1,NEL
        IF(NC(8,I).EQ.0) CYCLE
        MAX=0
        MIN=10000
        DO J=1,8
          IF (NC(J,I)) 204,205,203
  203     K=1
          GO TO 206
  204     K=-1
          NC(J,I)=-NC(J,I)
          GO TO 206
  205     NC(J,I)=NC(J-1,I)+K
  206     JJ=NC(J,I)
          IF(JJ.GT.MAX) MAX=JJ
          IF(JJ.LT.MIN) MIN=JJ
          XY(J,1)=XYT(1,JJ)
          XY(J,2)=XYT(2,JJ)
        END DO
        IF(MAX-MIN.GT.LH) LH=MAX-MIN
        DO J=2,8,2
          II=J+1
          IF(II.GT.8) II=1
          JJ=NC(J,I)
          IF((XY(J,1).NE.0.0).OR.(XY(J,2).NE.0.0)) CYCLE
          DO K=1,2
            XYT(K,JJ)=(XY(J-1,K)+XY(II,K))/2.0
          END DO
        END DO
      END DO
C
      IF(IPROB.EQ.1) CALL XTOP(IM,NL)
C
      LH=LH+1
C
      CALL XNAGR(NL,MM,NN,NEL)
C
      CALL XABND(IGRA,NNOD,MM,NN,NUF,IM,IPROB)
C
      CALL XBLOCK(NEL,ITEOR)
C
C     PRINT 1000,(J,(NC(I,J),I=1,8),J=1,NEL)
C     PRINT 1001,(J,(XYT(I,J),I=1,2),J=1,NNOD)
c      PRINT 1005,(J,LTEL(J),LST(J),PP(1,J),PT(1,J),J=1,NEP)
c      PRINT 1005,(J,LTEL(NEP+J),LST(NEP+J),PP(1,NEP+J),PT(1,NEP+J),
c     *	J=1,NEP1)
c      PRINT 1002,NZ,(MZD(I),I=1,NZ)
  401 OPEN(7,FILE='XY.SAM')
      DO 402 JI=1,NNOD
  402 WRITE (7,403) XYT(1,JI),XYT(2,JI)
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
 1000 FORMAT(9I5)
 1001 FORMAT(I5,2F10.2)
 1002 FORMAT(10I5)
 1005 format(3i5,2f8.2)
 1006 RETURN
      END SUBROUTINE XADATA
C*****
      BLOCK DATA
      COMMON/AL2/XYT(2,550)
      COMMON /UZLY/ K,T(25),AR(25)
      DATA K/5/,T/-.9061798,-.5384693,0.,.5384693,.9061798,20*0./,
     *AR/.2369269,.4786287,.5688889,.4786287,.2369269,20*0./
      DATA XYT/1100*0.0/
      END BLOCK DATA
C*****
      SUBROUTINE X4EL(iteor)
      COMMON
     */AL1A/LH
     */AL2/XYT(2,550)
     */AL4/XX(3),C(3),ISTR,PREC,NNOD
     */AL7/XY(8,2),AKT(8,8)
     */AL16/MZD(110),AMZ(110),NZ,NZ1
     */AL21/NEP,NEP1,LTEL(90),LST(90), PP(2,90),PT(2,90)
     */AL2A/NC(8,110)
C
      DIMENSION X(4),Y(4),NL(10),NU(30)
C
      P0=1.0
C
      M=2
      N=2
C
      X(1)=0.0
      X(2)=40.0
      X(3)=40.0
      X(4)=0.0
C
      Y(1)=0.0
      Y(2)=0.0
      Y(3)=40.0
      Y(4)=40.0
C
      NL(1)=1
      NL(2)=3
C
      NU(1)=1
      NU(2)=6
      NU(3)=9
      NU(4)=14
      NU(5)=17
C
      CALL XELR1(M,N,X,Y,NU,NL)
C
      NEL=4
      NNOD=21
      NNN=NNOD
C
      NEP=4
C
      LTEL(1)=1
      LTEL(2)=2
      LTEL(3)=3
      LTEL(4)=4
C
      LST(1)=1
      LST(2)=1
      LST(3)=5
      LST(4)=5
C
      PP(1,1)=P0
      PP(2,1)=P0
      PP(1,2)=P0
      PP(2,2)=P0
      PP(1,3)=-P0
      PP(2,3)=-P0
      PP(1,4)=-P0
      PP(2,4)=-P0
C
      NEP1=4
C
      LTEL(5)=1
      LTEL(6)=2
      LTEL(7)=3
      LTEL(8)=4
C
      LST(5)=7
      LST(6)=3
      LST(7)=7
      LST(8)=3
C
      PP(1,5)=P0
      PP(2,5)=P0
      PP(1,6)=-P0
      PP(2,6)=-P0
      PP(1,7)=P0
      PP(2,7)=P0
      PP(1,8)=-P0
      PP(2,8)=-P0
C
      NZ=1
      MZD(1)=11
C
      LH=0
C
      DO I=1,NEL
        IF(NC(8,I).EQ.0) CYCLE
        MAX=0
        MIN=10000
        DO J=1,8
          IF (NC(J,I)) 204,205,203
  203     K=1
          GO TO 206
  204     K=-1
          NC(J,I)=-NC(J,I)
          GO TO 206
  205     NC(J,I)=NC(J-1,I)+K
  206     JJ=NC(J,I)
          IF(JJ.GT.MAX) MAX=JJ
          IF(JJ.LT.MIN) MIN=JJ
          XY(J,1)=XYT(1,JJ)
          XY(J,2)=XYT(2,JJ)
        END DO
        IF(MAX-MIN.GT.LH) LH=MAX-MIN
        DO J=2,8,2
          II=J+1
          IF(II.GT.8) II=1
          JJ=NC(J,I)
          IF((XY(J,1).NE.0.0).OR.(XY(J,2).NE.0.0)) CYCLE
          DO K=1,2
            XYT(K,JJ)=(XY(J-1,K)+XY(II,K))/2.0
          END DO
        END DO
      END DO
C
      LH=LH+1
C
      CALL XBLOCK(NEL,ITEOR)
C
C     PRINT 1000,(J,(NC(I,J),I=1,8),J=1,NEL)
C     PRINT 1001,(J,(XYT(I,J),I=1,2),J=1,NNOD)
C      PRINT 1002,(J,LTEL(J),LST(J),PT(1,J),J=1,NEP)
C      PRINT 1002,(J,LTEL(NEP+J),LST(NEP+J),PT(1,NEP+J),J=1,NEP1)
 1000 FORMAT(9I5)
 1001 FORMAT(I5,2F10.2)
 1002 FORMAT(3I5,F5.2)
C
      RETURN
      END SUBROUTINE X4EL
C*****
      SUBROUTINE XCRACK(NL,NUF,NEL,NNOD,MM,NN,IM)
C----------------------------------------------------------------------
C                 èPOñEÑìPA îOPMàPOBAHàü TPEôàHõ
C----------------------------------------------------------------------
      COMMON
     */AL2/ XYT(2,550)
     */AL2A/ NC(8,110)
C
      REAL XYT1(2,550)/1100*0.0/
      INTEGER NC1(8,110)/880*0/
C
      DIMENSION NL(12,10),NUF(12,30),MM(4),NN(5)
C
      DO I=1,NEL
        DO J=1,8
          NC1(J,I)=NC(J,I)
        END DO
      END DO
C
      DO I=1,NNOD
        DO J=1,2
          XYT1(J,I)=XYT(J,I)
        END DO
      END DO
C
      NDEL=2*MM(1)+8
      IF(IM.EQ.1) NDEL=2*MM(1)+10
      NE=NL(7,1)
      NE1=NL(8,1)+3
      NE2=NL(8,1)+7
      NE3=NL(8,1)+11
      NE4=NL(8,1)+15
      NE5=NL(9,1)+MM(3)-1
      NE6=NE5+1
      NALL=2*MM(1)+2*MM(3)+17
      IF(IM.EQ.0) GO TO 21
      NE1=NL(8,1)+4
      NE2=NL(8,1)+9
      NE3=NL(8,1)+14
      NE4=NL(8,1)+19
      NALL=2*MM(1)+2*MM(3)+21
   21 CONTINUE
C
      DO I=NE6,NEL
        DO J=1,8
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      NE11=NL(8,1)-1
      DO I=NE,NE11
        DO J=1,3
          NC1(J,I)=NC(J,I)+NALL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      NE11=NL(8,1)-1
      DO I=NE,NE11
        DO J=4,8
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      DO I=NL(8,1),NE1
        DO J=3,5
          IF(I.EQ.NE1.AND.J.EQ.5) CYCLE
          NC1(J,I)=NC(J,I)+NALL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      DO I=NL(8,1),NE1
        DO J=1,8
          IF(J.EQ.3.OR.J.EQ.4.OR.J.EQ.5) CYCLE
          IF(I.EQ.NE1.AND.(J.EQ.6.OR.J.EQ.7)) CYCLE
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      NE11=NE1+1
      DO I=NE11,NE3
        DO J=1,8
          IF(I.EQ.NE2.AND.((J.EQ.5).OR.(J.EQ.6).OR.(J.EQ.7))) CYCLE
          IF(I.EQ.NE3.AND.((J.EQ.5).OR.(J.EQ.6).OR.(J.EQ.7))) CYCLE
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      NE11=NE3+1
      DO I=NE11,NE4
        DO J=1,8
          IF(J.EQ.1.OR.J.EQ.7.OR.J.EQ.8) CYCLE
          IF(I.EQ.NE4.AND.(J.EQ.5.OR.J.EQ.6)) CYCLE
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      DO I=NL(9,1),NE5
        DO J=1,8
          IF(J.EQ.1.OR.J.EQ.2.OR.J.EQ.3) CYCLE
          NC1(J,I)=NC(J,I)+NDEL
          NC1JI=NC1(J,I)
          NCJI=NC(J,I)
          DO K=1,2
            XYT1(K,NC1JI)=XYT(K,NCJI)
          END DO
        END DO
      END DO
C
      NNOD=NNOD+NDEL
C
      DO I=1,NEL
        DO J=1,8
          NC(J,I)=NC1(J,I)
        END DO
      END DO
C
      DO I=1,NNOD
        DO J=1,2
          XYT(J,I)=XYT1(J,I)
        END DO
      END DO
C
      NUF(7,1)=NUF(7,1)+NALL
      NUF(7,2)=NUF(7,2)+NDEL
      NUF(7,3)=NUF(7,3)+NDEL
C
      NN4=NN(4)*2+1
      DO I=1,NN4
        NUF(10,I)=NUF(10,I)+NDEL
      END DO
C
      RETURN
      END SUBROUTINE XCRACK
C*****
      SUBROUTINE XELR1(M,N,X,Y,NU,NL)
C----------------------------------------------------------------------
C èPOñEÑìPA PAáÅàBKà óETõPEXìÉOãúHOâ OÅãACTà CETKOâ KOHEóHõX ùãEMEHTOB
C----------------------------------------------------------------------
      COMMON
     */AL2/XYT(2,550)
     */AL2A/NC(8,110)
C
      DIMENSION X(4),Y(4),NU(30),NL(10)
C
      NN=(2*N+1)*(M+1)+(N+1)*M
C
      DO J=1,N
        NL1=NL(J)
        NJ1=NU(2*J-1)
        NJ2=NU(2*J)
        NJ3=NU(2*J+1)
        DO I=1,M
          NEL=NL1+I-1
          NJ=NJ1+2*(I-1)
          NC(1,NEL)=NJ
          NC(2,NEL)=NJ+1
          NC(3,NEL)=NJ+2
          NJ=NJ2+I-1
          NC(4,NEL)=NJ+1
          NC(8,NEL)=NJ
          NJ=NJ3+2*(I-1)
          NC(7,NEL)=NJ
          NC(6,NEL)=NJ+1
          NC(5,NEL)=NJ+2
        END DO
      END DO
C
      M2=M+1
      N2=N+1
C
      DO J=1,N2
        NJ1=NU(2*J-1)
        DO I=1,M2
          NJ=NJ1+2*(I-1)
          X1=X(1)+(X(4)-X(1))/N*(J-1)
          Y1=Y(1)+(Y(4)-Y(1))/N*(J-1)
          Y2=Y(2)+(Y(3)-Y(2))/N*(J-1)
          X2=X(2)+(X(3)-X(2))/N*(J-1)
          XYT(1,NJ)=X1+(X2-X1)/M*(I-1)
          XYT(2,NJ)=Y1+(Y2-Y1)/M*(I-1)
        END DO
      END DO
C
      RETURN
      END SUBROUTINE XELR1
C*****
      SUBROUTINE XELR2(NF,XZONE,YZONE,X,Y,NL,NUF,KLIN,ALFA)
C----------------------------------------------------------------------
C   PAáÅàEHàE HA KOHEóHõE ùãEMEHTõ HàÜHEÉO îPAÉMEHTA èOÑ TPEôàHOâ
C----------------------------------------------------------------------
      COMMON
     */AL2/XYT(2,550)
     */AL2A/NC(8,110)
C
      DIMENSION Y(4,12),X(4,12),NUF(12,30),NL(12,10)
C
      NEK=NL(NF,1)
      N1=NUF(NF,1)
      N2=NUF(NF,2)
      N3=NUF(NF,3)
C
      NC(1,NEK)=N3
      NC(2,NEK)=N2
      NC(3,NEK)=N1
      NC(4,NEK)=N2+4
      NC(5,NEK)=N2+5
      NC(6,NEK)=N2+1
      NC(7,NEK)=N3+2
      NC(8,NEK)=N3+1
C
      NC(1,NEK+1)=N3+2
      NC(2,NEK+1)=N2+1
      NC(3,NEK+1)=N2+5
      NC(4,NEK+1)=N2+6
      NC(5,NEK+1)=N2+7
      NC(6,NEK+1)=N2+2
      NC(7,NEK+1)=N3+4
      NC(8,NEK+1)=N3+3
C
      NC(1,NEK+2)=N3+4
      NC(2,NEK+2)=N2+2
      NC(3,NEK+2)=N2+7
      NC(4,NEK+2)=N2+8
      NC(5,NEK+2)=N2+9
      NC(6,NEK+2)=N2+3
      NC(7,NEK+2)=N3+6
      NC(8,NEK+2)=N3+5
C
      NC(1,NEK+3)=N3+6
      NC(2,NEK+3)=N2+3
      NC(3,NEK+3)=N2+9
      NC(4,NEK+3)=N2+10
      NC(5,NEK+3)=N3+8
      NC(6,NEK+3)=N3+8
      NC(7,NEK+3)=N3+8
      NC(8,NEK+3)=N3+7
C
      NC(1,NEK+4)=N1
      NC(2,NEK+4)=N1+1
      NC(3,NEK+4)=N1+2
      NC(4,NEK+4)=N2+14
      NC(5,NEK+4)=N2+15
      NC(6,NEK+4)=N2+11
      NC(7,NEK+4)=N2+5
      NC(8,NEK+4)=N2+4
C
      NC(1,NEK+5)=N2+5
      NC(2,NEK+5)=N2+11
      NC(3,NEK+5)=N2+15
      NC(4,NEK+5)=N2+16
      NC(5,NEK+5)=N2+17
      NC(6,NEK+5)=N2+12
      NC(7,NEK+5)=N2+7
      NC(8,NEK+5)=N2+6
C
      NC(1,NEK+6)=N2+7
      NC(2,NEK+6)=N2+12
      NC(3,NEK+6)=N2+17
      NC(4,NEK+6)=N2+18
      NC(5,NEK+6)=N2+19
      NC(6,NEK+6)=N2+13
      NC(7,NEK+6)=N2+9
      NC(8,NEK+6)=N2+8
C
      NC(1,NEK+7)=N2+9
      NC(2,NEK+7)=N2+13
      NC(3,NEK+7)=N2+19
      NC(4,NEK+7)=N2+20
      NC(5,NEK+7)=N3+8
      NC(6,NEK+7)=N3+8
      NC(7,NEK+7)=N3+8
      NC(8,NEK+7)=N2+10
C
      NC(1,NEK+8)=N1+2
      NC(2,NEK+8)=N1+3
      NC(3,NEK+8)=N1+4
      NC(4,NEK+8)=N2+24
      NC(5,NEK+8)=N2+25
      NC(6,NEK+8)=N2+21
      NC(7,NEK+8)=N2+15
      NC(8,NEK+8)=N2+14
C
      NC(1,NEK+9)=N2+15
      NC(2,NEK+9)=N2+21
      NC(3,NEK+9)=N2+25
      NC(4,NEK+9)=N2+26
      NC(5,NEK+9)=N2+27
      NC(6,NEK+9)=N2+22
      NC(7,NEK+9)=N2+17
      NC(8,NEK+9)=N2+16
C
      NC(1,NEK+10)=N2+17
      NC(2,NEK+10)=N2+22
      NC(3,NEK+10)=N2+27
      NC(4,NEK+10)=N2+28
      NC(5,NEK+10)=N2+29
      NC(6,NEK+10)=N2+23
      NC(7,NEK+10)=N2+19
      NC(8,NEK+10)=N2+18
C
      NC(1,NEK+11)=N2+19
      NC(2,NEK+11)=N2+23
      NC(3,NEK+11)=N2+29
      NC(4,NEK+11)=N2+30
      NC(5,NEK+11)=N3+8
      NC(6,NEK+11)=N3+8
      NC(7,NEK+11)=N3+8
      NC(8,NEK+11)=N2+20
C
      NC(1,NEK+12)=N1+4
      NC(2,NEK+12)=N2+34
      NC(3,NEK+12)=N3+16
      NC(4,NEK+12)=N3+15
      NC(5,NEK+12)=N3+14
      NC(6,NEK+12)=N2+33
      NC(7,NEK+12)=N2+25
      NC(8,NEK+12)=N2+24
C
      NC(1,NEK+13)=N2+25
      NC(2,NEK+13)=N2+33
      NC(3,NEK+13)=N3+14
      NC(4,NEK+13)=N3+13
      NC(5,NEK+13)=N3+12
      NC(6,NEK+13)=N2+32
      NC(7,NEK+13)=N2+27
      NC(8,NEK+13)=N2+26
C
      NC(1,NEK+14)=N2+27
      NC(2,NEK+14)=N2+32
      NC(3,NEK+14)=N3+12
      NC(4,NEK+14)=N3+11
      NC(5,NEK+14)=N3+10
      NC(6,NEK+14)=N2+31
      NC(7,NEK+14)=N2+29
      NC(8,NEK+14)=N2+28
C
      NC(1,NEK+15)=N2+29
      NC(2,NEK+15)=N2+31
      NC(3,NEK+15)=N3+10
      NC(4,NEK+15)=N3+9
      NC(5,NEK+15)=N3+8
      NC(6,NEK+15)=N3+8
      NC(7,NEK+15)=N3+8
      NC(8,NEK+15)=N2+30
C
      XYT(1,N1)=X(1,NF)
      XYT(1,N1+2)=X(1,NF)+XZONE/2.
      XYT(1,N1+4)=X(2,NF)
      XYT(1,N2+5)=X(1,NF)+XZONE/4.
      XYT(1,N2+15)=X(1,NF)+XZONE/2.
      XYT(1,N2+25)=X(2,NF)-XZONE/4.
      XYT(1,N2+7)=X(1,NF)+3.*XZONE/8.
      XYT(1,N2+17)=X(1,NF)+XZONE/2.
      XYT(1,N2+27)=X(1,NF)+5.*XZONE/8.
      XYT(1,N2+9)=X(1,NF)+7*XZONE/16.
      XYT(1,N2+19)=X(1,NF)+XZONE/2.
      XYT(1,N2+29)=X(1,NF)+9*XZONE/16.
      XYT(1,N3)=X(4,NF)
      XYT(1,N3+2)=X(4,NF)+XZONE/4.
      XYT(1,N3+4)=X(4,NF)+3*XZONE/8.
      XYT(1,N3+6)=X(4,NF)+7*XZONE/16.
      XYT(1,N3+8)=X(4,NF)+XZONE/2.
      XYT(1,N3+10)=X(4,NF)+9*XZONE/16.
      XYT(1,N3+12)=X(4,NF)+5*XZONE/8.
      XYT(1,N3+14)=X(3,NF)-XZONE/4.
      XYT(1,N3+16)=X(3,NF)
C
      XYT(2,N1)=Y(1,NF)
      XYT(2,N1+2)=Y(1,NF)
      XYT(2,N1+4)=Y(1,NF)
      XYT(2,N2+5)=Y(1,NF)+YZONE/4.
      XYT(2,N2+15)=Y(1,NF)+YZONE/4.
      XYT(2,N2+25)=Y(1,NF)+YZONE/4.
      XYT(2,N2+7)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+17)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+27)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+9)=Y(1,NF)+7*YZONE/16.
      XYT(2,N2+19)=Y(1,NF)+7*YZONE/16.
      XYT(2,N2+29)=Y(1,NF)+7*YZONE/16.
      XYT(2,N3)=Y(3,NF)         -KLIN*XZONE/2.*TAN(ALFA/2.)*8./8.
      XYT(2,N3+2)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*4./8.
      XYT(2,N3+4)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*2./8.
      XYT(2,N3+6)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*1./8.
      XYT(2,N3+8)=Y(3,NF)
      XYT(2,N3+10)=Y(3,NF)
      XYT(2,N3+12)=Y(3,NF)
      XYT(2,N3+14)=Y(3,NF)
      XYT(2,N3+16)=Y(3,NF)
C
      RETURN
      END SUBROUTINE XELR2
C*****
      SUBROUTINE XELR21(NF,XZONE,YZONE,X,Y,NL,NUF,KLIN,ALFA)
C----------------------------------------------------------------------
C   PAáÅàEHàE HA KOHEóHõE ùãEMEHTõ îPAÉMEHTA èOÑ TPEôàHOâ
C----------------------------------------------------------------------
      COMMON
     */AL2/XYT(2,550)
     */AL2A/NC(8,110)
C
      DIMENSION Y(4,12),X(4,12),NUF(12,30),NL(12,10)
C
      NEK=NL(NF,1)
      N1=NUF(NF,1)
      N2=NUF(NF,2)
      N3=NUF(NF,3)
C
      NC(1,NEK)=N3
      NC(2,NEK)=N2
      NC(3,NEK)=N1
      NC(4,NEK)=N2+5
      NC(5,NEK)=N2+6
      NC(6,NEK)=N2+1
      NC(7,NEK)=N3+2
      NC(8,NEK)=N3+1
C
      NC(1,NEK+1)=N3+2
      NC(2,NEK+1)=N2+1
      NC(3,NEK+1)=N2+6
      NC(4,NEK+1)=N2+7
      NC(5,NEK+1)=N2+8
      NC(6,NEK+1)=N2+2
      NC(7,NEK+1)=N3+4
      NC(8,NEK+1)=N3+3
C
      NC(1,NEK+2)=N3+4
      NC(2,NEK+2)=N2+2
      NC(3,NEK+2)=N2+8
      NC(4,NEK+2)=N2+9
      NC(5,NEK+2)=N2+10
      NC(6,NEK+2)=N2+3
      NC(7,NEK+2)=N3+6
      NC(8,NEK+2)=N3+5
C
      NC(1,NEK+3)=N3+6
      NC(2,NEK+3)=N2+3
      NC(3,NEK+3)=N2+10
      NC(4,NEK+3)=N2+11
      NC(5,NEK+3)=N2+12
      NC(6,NEK+3)=N2+4
      NC(7,NEK+3)=N3+8
      NC(8,NEK+3)=N3+7
C
      NC(1,NEK+4)=N3+8
      NC(2,NEK+4)=N2+4
      NC(3,NEK+4)=N2+12
      NC(4,NEK+4)=N2+13
      NC(5,NEK+4)=N3+10
      NC(6,NEK+4)=N3+10
      NC(7,NEK+4)=N3+10
      NC(8,NEK+4)=N3+9
C
      NC(1,NEK+5)=N1
      NC(2,NEK+5)=N1+1
      NC(3,NEK+5)=N1+2
      NC(4,NEK+5)=N2+18
      NC(5,NEK+5)=N2+19
      NC(6,NEK+5)=N2+14
      NC(7,NEK+5)=N2+6
      NC(8,NEK+5)=N2+5
C
      NC(1,NEK+6)=N2+6
      NC(2,NEK+6)=N2+14
      NC(3,NEK+6)=N2+19
      NC(4,NEK+6)=N2+20
      NC(5,NEK+6)=N2+21
      NC(6,NEK+6)=N2+15
      NC(7,NEK+6)=N2+8
      NC(8,NEK+6)=N2+7
C
      NC(1,NEK+7)=N2+8
      NC(2,NEK+7)=N2+15
      NC(3,NEK+7)=N2+21
      NC(4,NEK+7)=N2+22
      NC(5,NEK+7)=N2+23
      NC(6,NEK+7)=N2+16
      NC(7,NEK+7)=N2+10
      NC(8,NEK+7)=N2+9
C
      NC(1,NEK+8)=N2+10
      NC(2,NEK+8)=N2+16
      NC(3,NEK+8)=N2+23
      NC(4,NEK+8)=N2+24
      NC(5,NEK+8)=N2+25
      NC(6,NEK+8)=N2+17
      NC(7,NEK+8)=N2+12
      NC(8,NEK+8)=N2+11
C
      NC(1,NEK+9)=N2+12
      NC(2,NEK+9)=N2+17
      NC(3,NEK+9)=N2+25
      NC(4,NEK+9)=N2+26
      NC(5,NEK+9)=N3+10
      NC(6,NEK+9)=N3+10
      NC(7,NEK+9)=N3+10
      NC(8,NEK+9)=N2+13
C
      NC(1,NEK+10)=N1+2
      NC(2,NEK+10)=N1+3
      NC(3,NEK+10)=N1+4
      NC(4,NEK+10)=N2+31
      NC(5,NEK+10)=N2+32
      NC(6,NEK+10)=N2+27
      NC(7,NEK+10)=N2+19
      NC(8,NEK+10)=N2+18
C
      NC(1,NEK+11)=N2+19
      NC(2,NEK+11)=N2+27
      NC(3,NEK+11)=N2+32
      NC(4,NEK+11)=N2+33
      NC(5,NEK+11)=N2+34
      NC(6,NEK+11)=N2+28
      NC(7,NEK+11)=N2+21
      NC(8,NEK+11)=N2+20
C
      NC(1,NEK+12)=N2+21
      NC(2,NEK+12)=N2+28
      NC(3,NEK+12)=N2+34
      NC(4,NEK+12)=N2+35
      NC(5,NEK+12)=N2+36
      NC(6,NEK+12)=N2+29
      NC(7,NEK+12)=N2+23
      NC(8,NEK+12)=N2+22
C
      NC(1,NEK+13)=N2+23
      NC(2,NEK+13)=N2+29
      NC(3,NEK+13)=N2+36
      NC(4,NEK+13)=N2+37
      NC(5,NEK+13)=N2+38
      NC(6,NEK+13)=N2+30
      NC(7,NEK+13)=N2+25
      NC(8,NEK+13)=N2+24
C
      NC(1,NEK+14)=N2+25
      NC(2,NEK+14)=N2+30
      NC(3,NEK+14)=N2+38
      NC(4,NEK+14)=N2+39
      NC(5,NEK+14)=N3+10
      NC(6,NEK+14)=N3+10
      NC(7,NEK+14)=N3+10
      NC(8,NEK+14)=N2+26
C
      NC(1,NEK+15)=N1+4
      NC(2,NEK+15)=N2+44
      NC(3,NEK+15)=N3+20
      NC(4,NEK+15)=N3+19
      NC(5,NEK+15)=N3+18
      NC(6,NEK+15)=N2+43
      NC(7,NEK+15)=N2+32
      NC(8,NEK+15)=N2+31
C
      NC(1,NEK+16)=N2+32
      NC(2,NEK+16)=N2+43
      NC(3,NEK+16)=N3+18
      NC(4,NEK+16)=N3+17
      NC(5,NEK+16)=N3+16
      NC(6,NEK+16)=N2+42
      NC(7,NEK+16)=N2+34
      NC(8,NEK+16)=N2+33
C
      NC(1,NEK+17)=N2+34
      NC(2,NEK+17)=N2+42
      NC(3,NEK+17)=N3+16
      NC(4,NEK+17)=N3+15
      NC(5,NEK+17)=N3+14
      NC(6,NEK+17)=N2+41
      NC(7,NEK+17)=N2+36
      NC(8,NEK+17)=N2+35
C
      NC(1,NEK+18)=N2+36
      NC(2,NEK+18)=N2+41
      NC(3,NEK+18)=N3+14
      NC(4,NEK+18)=N3+13
      NC(5,NEK+18)=N3+12
      NC(6,NEK+18)=N2+40
      NC(7,NEK+18)=N2+38
      NC(8,NEK+18)=N2+37
C
      NC(1,NEK+19)=N2+38
      NC(2,NEK+19)=N2+40
      NC(3,NEK+19)=N3+12
      NC(4,NEK+19)=N3+11
      NC(5,NEK+19)=N3+10
      NC(6,NEK+19)=N3+10
      NC(7,NEK+19)=N3+10
      NC(8,NEK+19)=N2+39
C
      XYT(1,N1)=X(1,NF)
      XYT(1,N1+2)=X(1,NF)+XZONE/2.
      XYT(1,N1+4)=X(2,NF)
      XYT(1,N2+6)=X(1,NF)+XZONE/4.
      XYT(1,N2+19)=X(1,NF)+XZONE/2.
      XYT(1,N2+32)=X(2,NF)-XZONE/4.
      XYT(1,N2+8)=X(1,NF)+3.*XZONE/8.
      XYT(1,N2+21)=X(1,NF)+XZONE/2.
      XYT(1,N2+34)=X(1,NF)+5.*XZONE/8.
      XYT(1,N2+10)=X(1,NF)+7*XZONE/16.
      XYT(1,N2+23)=X(1,NF)+XZONE/2.
      XYT(1,N2+36)=X(1,NF)+9*XZONE/16.
      XYT(1,N2+12)=X(1,NF)+15*XZONE/32.
      XYT(1,N2+25)=X(1,NF)+XZONE/2.
      XYT(1,N2+38)=X(1,NF)+17*XZONE/32.
      XYT(1,N3)=X(4,NF)
      XYT(1,N3+2)=X(4,NF)+XZONE/4.
      XYT(1,N3+4)=X(4,NF)+3*XZONE/8.
      XYT(1,N3+6)=X(4,NF)+7*XZONE/16.
      XYT(1,N3+8)=X(4,NF)+15*XZONE/32.
      XYT(1,N3+10)=X(4,NF)+XZONE/2.
      XYT(1,N3+12)=X(4,NF)+17*XZONE/32.
      XYT(1,N3+14)=X(4,NF)+9*XZONE/16.
      XYT(1,N3+16)=X(4,NF)+5*XZONE/8.
      XYT(1,N3+18)=X(3,NF)-XZONE/4.
      XYT(1,N3+20)=X(3,NF)
C
      XYT(2,N1)=Y(1,NF)
      XYT(2,N1+2)=Y(1,NF)
      XYT(2,N1+4)=Y(1,NF)
      XYT(2,N2+6)=Y(1,NF)+YZONE/4.
      XYT(2,N2+19)=Y(1,NF)+YZONE/4.
      XYT(2,N2+32)=Y(1,NF)+YZONE/4.
      XYT(2,N2+8)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+21)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+34)=Y(1,NF)+3.*YZONE/8.
      XYT(2,N2+10)=Y(1,NF)+7*YZONE/16.
      XYT(2,N2+23)=Y(1,NF)+7*YZONE/16.
      XYT(2,N2+36)=Y(1,NF)+7*YZONE/16.
      XYT(2,N2+12)=Y(1,NF)+15*YZONE/32.
      XYT(2,N2+25)=Y(1,NF)+15*YZONE/32.
      XYT(2,N2+38)=Y(1,NF)+15*YZONE/32.
      XYT(2,N3)=Y(3,NF)         -KLIN*XZONE/2.*TAN(ALFA/2.)*8./8.
      XYT(2,N3+2)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*4./8.
      XYT(2,N3+4)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*2./8.
      XYT(2,N3+6)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*1./8.
      XYT(2,N3+8)=Y(3,NF)       -KLIN*XZONE/2.*TAN(ALFA/2.)*1./16.
      XYT(2,N3+10)=Y(3,NF)
      XYT(2,N3+12)=Y(3,NF)
      XYT(2,N3+14)=Y(3,NF)
      XYT(2,N3+16)=Y(3,NF)
      XYT(2,N3+18)=Y(3,NF)
      XYT(2,N3+20)=Y(3,NF)
C
      RETURN
      END SUBROUTINE XELR21

