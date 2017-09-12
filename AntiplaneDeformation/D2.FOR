      SUBROUTINE XELR3(NF,XZONE,YZONE,X,Y,NL,NUF,KLIN,ALFA)
C----------------------------------------------------------------------
C   PAáÅàEHàE HA KOHEóHõE ùãEMEHTõ BEPXHEÉO îPAÉMEHTA HAÑ TPEôàHOâ
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
      NC(4,NEK)=N1+1
      NC(5,NEK)=N1+2
      NC(6,NEK)=N2+1
      NC(7,NEK)=N2+5
      NC(8,NEK)=N2+4
C
      NC(1,NEK+1)=N2+5
      NC(2,NEK+1)=N2+1
      NC(3,NEK+1)=N1+2
      NC(4,NEK+1)=N1+3
      NC(5,NEK+1)=N1+4
      NC(6,NEK+1)=N2+2
      NC(7,NEK+1)=N2+7
      NC(8,NEK+1)=N2+6
C
      NC(1,NEK+2)=N2+7
      NC(2,NEK+2)=N2+2
      NC(3,NEK+2)=N1+4
      NC(4,NEK+2)=N1+5
      NC(5,NEK+2)=N1+6
      NC(6,NEK+2)=N2+3
      NC(7,NEK+2)=N2+9
      NC(8,NEK+2)=N2+8
C
      NC(1,NEK+3)=N2+9
      NC(2,NEK+3)=N2+3
      NC(3,NEK+3)=N1+6
      NC(4,NEK+3)=N1+7
      NC(5,NEK+3)=N1+8
      NC(6,NEK+3)=N1+8
      NC(7,NEK+3)=N1+8
      NC(8,NEK+3)=N2+10
C
      NC(1,NEK+4)=N3+2
      NC(2,NEK+4)=N3+1
      NC(3,NEK+4)=N3
      NC(4,NEK+4)=N2+4
      NC(5,NEK+4)=N2+5
      NC(6,NEK+4)=N2+11
      NC(7,NEK+4)=N2+15
      NC(8,NEK+4)=N2+14
C
      NC(1,NEK+5)=N2+15
      NC(2,NEK+5)=N2+11
      NC(3,NEK+5)=N2+5
      NC(4,NEK+5)=N2+6
      NC(5,NEK+5)=N2+7
      NC(6,NEK+5)=N2+12
      NC(7,NEK+5)=N2+17
      NC(8,NEK+5)=N2+16
C
      NC(1,NEK+6)=N2+17
      NC(2,NEK+6)=N2+12
      NC(3,NEK+6)=N2+7
      NC(4,NEK+6)=N2+8
      NC(5,NEK+6)=N2+9
      NC(6,NEK+6)=N2+13
      NC(7,NEK+6)=N2+19
      NC(8,NEK+6)=N2+18
C
      NC(1,NEK+7)=N2+19
      NC(2,NEK+7)=N2+13
      NC(3,NEK+7)=N2+9
      NC(4,NEK+7)=N2+10
      NC(5,NEK+7)=N1+8
      NC(6,NEK+7)=N1+8
      NC(7,NEK+7)=N1+8
      NC(8,NEK+7)=N2+20
C
      NC(1,NEK+8)=N3+4
      NC(2,NEK+8)=N3+3
      NC(3,NEK+8)=N3+2
      NC(4,NEK+8)=N2+14
      NC(5,NEK+8)=N2+15
      NC(6,NEK+8)=N2+21
      NC(7,NEK+8)=N2+25
      NC(8,NEK+8)=N2+24
C
      NC(1,NEK+9)=N2+25
      NC(2,NEK+9)=N2+21
      NC(3,NEK+9)=N2+15
      NC(4,NEK+9)=N2+16
      NC(5,NEK+9)=N2+17
      NC(6,NEK+9)=N2+22
      NC(7,NEK+9)=N2+27
      NC(8,NEK+9)=N2+26
C
      NC(1,NEK+10)=N2+27
      NC(2,NEK+10)=N2+22
      NC(3,NEK+10)=N2+17
      NC(4,NEK+10)=N2+18
      NC(5,NEK+10)=N2+19
      NC(6,NEK+10)=N2+23
      NC(7,NEK+10)=N2+29
      NC(8,NEK+10)=N2+28
C
      NC(1,NEK+11)=N2+29
      NC(2,NEK+11)=N2+23
      NC(3,NEK+11)=N2+19
      NC(4,NEK+11)=N2+20
      NC(5,NEK+11)=N1+8
      NC(6,NEK+11)=N1+8
      NC(7,NEK+11)=N1+8
      NC(8,NEK+11)=N2+30
C
      NC(1,NEK+12)=N1+16
      NC(2,NEK+12)=N2+34
      NC(3,NEK+12)=N3+4
      NC(4,NEK+12)=N2+24
      NC(5,NEK+12)=N2+25
	NC(6,NEK+12)=N2+33
      NC(7,NEK+12)=N1+14
      NC(8,NEK+12)=N1+15
C
      NC(1,NEK+13)=N1+14
      NC(2,NEK+13)=N2+33
      NC(3,NEK+13)=N2+25
      NC(4,NEK+13)=N2+26
      NC(5,NEK+13)=N2+27
      NC(6,NEK+13)=N2+32
      NC(7,NEK+13)=N1+12
      NC(8,NEK+13)=N1+13
C
      NC(1,NEK+14)=N1+12
      NC(2,NEK+14)=N2+32
      NC(3,NEK+14)=N2+27
      NC(4,NEK+14)=N2+28
      NC(5,NEK+14)=N2+29
      NC(6,NEK+14)=N2+31
      NC(7,NEK+14)=N1+10
      NC(8,NEK+14)=N1+11
C
      NC(1,NEK+15)=N1+10
      NC(2,NEK+15)=N2+31
      NC(3,NEK+15)=N2+29
      NC(4,NEK+15)=N2+30
      NC(5,NEK+15)=N1+8
      NC(6,NEK+15)=N1+8
      NC(7,NEK+15)=N1+8
      NC(8,NEK+15)=N1+9
C
      XYT(1,N3)=X(1,NF)
      XYT(1,N3+2)=X(1,NF)+XZONE/2.
      XYT(1,N3+4)=X(2,NF)
      XYT(1,N2+5)=X(1,NF)+XZONE/4.
      XYT(1,N2+15)=X(1,NF)+XZONE/2.
      XYT(1,N2+25)=X(2,NF)-XZONE/4.
      XYT(1,N2+7)=X(1,NF)+3.*XZONE/8.
      XYT(1,N2+17)=X(1,NF)+XZONE/2.
      XYT(1,N2+27)=X(1,NF)+5.*XZONE/8.
      XYT(1,N2+9)=X(1,NF)+7*XZONE/16.
      XYT(1,N2+19)=X(1,NF)+XZONE/2.
      XYT(1,N2+29)=X(1,NF)+9*XZONE/16.
      XYT(1,N1)=X(4,NF)
      XYT(1,N1+2)=X(4,NF)+XZONE/4.
      XYT(1,N1+4)=X(4,NF)+3*XZONE/8.
      XYT(1,N1+6)=X(4,NF)+7*XZONE/16.
      XYT(1,N1+8)=X(4,NF)+XZONE/2.
      XYT(1,N1+10)=X(4,NF)+9*XZONE/16.
      XYT(1,N1+12)=X(4,NF)+5*XZONE/8.
      XYT(1,N1+14)=X(3,NF)-XZONE/4.
      XYT(1,N1+16)=X(3,NF)
C
      XYT(2,N3)=Y(4,NF)
      XYT(2,N3+2)=Y(4,NF)
      XYT(2,N3+4)=Y(4,NF)
      XYT(2,N2+5)=Y(4,NF)-YZONE/4.
      XYT(2,N2+15)=Y(4,NF)-YZONE/4.
      XYT(2,N2+25)=Y(4,NF)-YZONE/4.
      XYT(2,N2+7)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+17)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+27)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+9)=Y(4,NF)-7*YZONE/16.
      XYT(2,N2+19)=Y(4,NF)-7*YZONE/16.
      XYT(2,N2+29)=Y(4,NF)-7*YZONE/16.
      XYT(2,N1)=Y(2,NF)           +KLIN*XZONE/2.*TAN(ALFA/2.)*8./8.
      XYT(2,N1+2)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*4./8.
      XYT(2,N1+4)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*2./8.
      XYT(2,N1+6)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*1./8.
      XYT(2,N1+8)=Y(2,NF)
      XYT(2,N1+10)=Y(2,NF)
      XYT(2,N1+12)=Y(2,NF)
      XYT(2,N1+14)=Y(2,NF)
      XYT(2,N1+16)=Y(2,NF)
C
      RETURN
      END
      SUBROUTINE XELR31(NF,XZONE,YZONE,X,Y,NL,NUF,KLIN,ALFA)
C----------------------------------------------------------------------
C   PAáÅàEHàE HA KOHEóHõE ùãEMEHTõ BEPXHEÉO îPAÉMEHTA HAÑ TPEôàHOâ
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
      NC(4,NEK)=N1+1
      NC(5,NEK)=N1+2
      NC(6,NEK)=N2+1
      NC(7,NEK)=N2+6
      NC(8,NEK)=N2+5
C
      NC(1,NEK+1)=N2+6
      NC(2,NEK+1)=N2+1
      NC(3,NEK+1)=N1+2
      NC(4,NEK+1)=N1+3
      NC(5,NEK+1)=N1+4
      NC(6,NEK+1)=N2+2
      NC(7,NEK+1)=N2+8
      NC(8,NEK+1)=N2+7
C
      NC(1,NEK+2)=N2+8
      NC(2,NEK+2)=N2+2
      NC(3,NEK+2)=N1+4
      NC(4,NEK+2)=N1+5
      NC(5,NEK+2)=N1+6
      NC(6,NEK+2)=N2+3
      NC(7,NEK+2)=N2+10
      NC(8,NEK+2)=N2+9
C
      NC(1,NEK+3)=N2+10
      NC(2,NEK+3)=N2+3
      NC(3,NEK+3)=N1+6
      NC(4,NEK+3)=N1+7
      NC(5,NEK+3)=N1+8
      NC(6,NEK+3)=N2+4
      NC(7,NEK+3)=N2+12
      NC(8,NEK+3)=N2+11
C
      NC(1,NEK+4)=N2+12
      NC(2,NEK+4)=N2+4
      NC(3,NEK+4)=N1+8
      NC(4,NEK+4)=N1+9
      NC(5,NEK+4)=N1+10
      NC(6,NEK+4)=N1+10
      NC(7,NEK+4)=N1+10
      NC(8,NEK+4)=N2+13
C
      NC(1,NEK+5)=N3+2
      NC(2,NEK+5)=N3+1
      NC(3,NEK+5)=N3
      NC(4,NEK+5)=N2+5
      NC(5,NEK+5)=N2+6
      NC(6,NEK+5)=N2+14
      NC(7,NEK+5)=N2+19
      NC(8,NEK+5)=N2+18
C
      NC(1,NEK+6)=N2+19
      NC(2,NEK+6)=N2+14
      NC(3,NEK+6)=N2+6
      NC(4,NEK+6)=N2+7
      NC(5,NEK+6)=N2+8
      NC(6,NEK+6)=N2+15
      NC(7,NEK+6)=N2+21
	NC(8,NEK+6)=N2+20
C
      NC(1,NEK+7)=N2+21
      NC(2,NEK+7)=N2+15
      NC(3,NEK+7)=N2+8
      NC(4,NEK+7)=N2+9
      NC(5,NEK+7)=N2+10
      NC(6,NEK+7)=N2+16
      NC(7,NEK+7)=N2+23
      NC(8,NEK+7)=N2+22
C
      NC(1,NEK+8)=N2+23
      NC(2,NEK+8)=N2+16
      NC(3,NEK+8)=N2+10
      NC(4,NEK+8)=N2+11
      NC(5,NEK+8)=N2+12
      NC(6,NEK+8)=N2+17
      NC(7,NEK+8)=N2+25
      NC(8,NEK+8)=N2+24
C
      NC(1,NEK+9)=N2+25
      NC(2,NEK+9)=N2+17
      NC(3,NEK+9)=N2+12
      NC(4,NEK+9)=N2+13
      NC(5,NEK+9)=N1+10
      NC(6,NEK+9)=N1+10
      NC(7,NEK+9)=N1+10
      NC(8,NEK+9)=N2+26
C
      NC(1,NEK+10)=N3+4
      NC(2,NEK+10)=N3+3
      NC(3,NEK+10)=N3+2
      NC(4,NEK+10)=N2+18
      NC(5,NEK+10)=N2+19
	NC(6,NEK+10)=N2+27
      NC(7,NEK+10)=N2+32
      NC(8,NEK+10)=N2+31
C
      NC(1,NEK+11)=N2+32
      NC(2,NEK+11)=N2+27
      NC(3,NEK+11)=N2+19
      NC(4,NEK+11)=N2+20
      NC(5,NEK+11)=N2+21
      NC(6,NEK+11)=N2+28
      NC(7,NEK+11)=N2+34
      NC(8,NEK+11)=N2+33
C
      NC(1,NEK+12)=N2+34
      NC(2,NEK+12)=N2+28
      NC(3,NEK+12)=N2+21
      NC(4,NEK+12)=N2+22
      NC(5,NEK+12)=N2+23
      NC(6,NEK+12)=N2+29
      NC(7,NEK+12)=N2+36
      NC(8,NEK+12)=N2+35
C
      NC(1,NEK+13)=N2+36
      NC(2,NEK+13)=N2+29
      NC(3,NEK+13)=N2+23
      NC(4,NEK+13)=N2+24
      NC(5,NEK+13)=N2+25
      NC(6,NEK+13)=N2+30
      NC(7,NEK+13)=N2+38
      NC(8,NEK+13)=N2+37
C
      NC(1,NEK+14)=N2+38
      NC(2,NEK+14)=N2+30
      NC(3,NEK+14)=N2+25
      NC(4,NEK+14)=N2+26
      NC(5,NEK+14)=N1+10
      NC(6,NEK+14)=N1+10
      NC(7,NEK+14)=N1+10
      NC(8,NEK+14)=N2+39
C
      NC(1,NEK+15)=N1+20
      NC(2,NEK+15)=N2+44
      NC(3,NEK+15)=N3+4
      NC(4,NEK+15)=N2+31
      NC(5,NEK+15)=N2+32
      NC(6,NEK+15)=N2+43
      NC(7,NEK+15)=N1+18
      NC(8,NEK+15)=N1+19
C
      NC(1,NEK+16)=N1+18
      NC(2,NEK+16)=N2+43
      NC(3,NEK+16)=N2+32
      NC(4,NEK+16)=N2+33
      NC(5,NEK+16)=N2+34
      NC(6,NEK+16)=N2+42
      NC(7,NEK+16)=N1+16
      NC(8,NEK+16)=N1+17
C
      NC(1,NEK+17)=N1+16
      NC(2,NEK+17)=N2+42
      NC(3,NEK+17)=N2+34
      NC(4,NEK+17)=N2+35
      NC(5,NEK+17)=N2+36
      NC(6,NEK+17)=N2+41
      NC(7,NEK+17)=N1+14
      NC(8,NEK+17)=N1+15
C
      NC(1,NEK+18)=N1+14
      NC(2,NEK+18)=N2+41
      NC(3,NEK+18)=N2+36
      NC(4,NEK+18)=N2+37
      NC(5,NEK+18)=N2+38
      NC(6,NEK+18)=N2+40
      NC(7,NEK+18)=N1+12
      NC(8,NEK+18)=N1+13
C
      NC(1,NEK+19)=N1+12
      NC(2,NEK+19)=N2+40
      NC(3,NEK+19)=N2+38
      NC(4,NEK+19)=N2+39
      NC(5,NEK+19)=N1+10
      NC(6,NEK+19)=N1+10
      NC(7,NEK+19)=N1+10
      NC(8,NEK+19)=N1+11
C
      XYT(1,N3)=X(1,NF)
      XYT(1,N3+2)=X(1,NF)+XZONE/2.
      XYT(1,N3+4)=X(2,NF)
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
      XYT(1,N1)=X(4,NF)
      XYT(1,N1+2)=X(4,NF)+XZONE/4.
      XYT(1,N1+4)=X(4,NF)+3*XZONE/8.
      XYT(1,N1+6)=X(4,NF)+7*XZONE/16.
      XYT(1,N1+8)=X(4,NF)+15*XZONE/32.
      XYT(1,N1+10)=X(4,NF)+XZONE/2.
      XYT(1,N1+12)=X(4,NF)+17*XZONE/32.
      XYT(1,N1+14)=X(4,NF)+9*XZONE/16.
      XYT(1,N1+16)=X(4,NF)+5*XZONE/8.
      XYT(1,N1+18)=X(3,NF)-XZONE/4.
      XYT(1,N1+20)=X(3,NF)
C
      XYT(2,N3)=Y(4,NF)
      XYT(2,N3+2)=Y(4,NF)
      XYT(2,N3+4)=Y(4,NF)
      XYT(2,N2+6)=Y(4,NF)-YZONE/4.
      XYT(2,N2+19)=Y(4,NF)-YZONE/4.
      XYT(2,N2+32)=Y(4,NF)-YZONE/4.
      XYT(2,N2+8)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+21)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+34)=Y(4,NF)-3.*YZONE/8.
      XYT(2,N2+10)=Y(4,NF)-7*YZONE/16.
      XYT(2,N2+23)=Y(4,NF)-7*YZONE/16.
      XYT(2,N2+36)=Y(4,NF)-7*YZONE/16.
      XYT(2,N2+12)=Y(4,NF)-15*YZONE/32.
      XYT(2,N2+25)=Y(4,NF)-15*YZONE/32.
      XYT(2,N2+38)=Y(4,NF)-15*YZONE/32.
      XYT(2,N1)=Y(2,NF)           +KLIN*XZONE/2.*TAN(ALFA/2.)*8./8.
      XYT(2,N1+2)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*4./8.
      XYT(2,N1+4)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*2./8.
      XYT(2,N1+6)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*1./8.
      XYT(2,N1+8)=Y(2,NF)         +KLIN*XZONE/2.*TAN(ALFA/2.)*1./16.
      XYT(2,N1+10)=Y(2,NF)
      XYT(2,N1+12)=Y(2,NF)
      XYT(2,N1+14)=Y(2,NF)
      XYT(2,N1+16)=Y(2,NF)
      XYT(2,N1+18)=Y(2,NF)
      XYT(2,N1+20)=Y(2,NF)
C
      RETURN
      END
      SUBROUTINE XNAGR(NL,MM,NN,NEL)
C----------------------------------------------------------------------
C        îOPMàPOBAHàE HAÉPìáKà HA ÉPAHàñE
C        PP,PT Ñãü OÅOÅôEHHO-èãOCKOÉO HÑC
C        PP  Ñãü AHTàèãOCKOâ ÑEîOPMAñàà
C----------------------------------------------------------------------
      COMMON
     */AL21/NEP,NEP1,LTEL(90),LST(90), PP(2,90),PT(2,90)
C
      DIMENSION NL(12,10),NN(5),MM(4)
C
      MMM=0
      DO 1 I=1,3
         MMM=MMM+MM(I)
    1 CONTINUE
C
      NEP=2*MMM
C
      DO 2 I=1,MMM
C  ----------------------->  îOPMàPOBAHàE HAÉPìáKà CHàáì
         LTEL(I)=I
         LST(I)=1
         PP(1,I)=1.0
         PP(2,I)=1.0
         PT(1,I)=0.0
         PT(2,I)=0.0
C  ----------------------->  îOPMàPOBAHàE HAÉPìáKà CBEPXì
            LTEL(I+MMM)=NEL-MMM+I
            LST(I+MMM)=5
            PP(1,I+MMM)=-1.0
            PP(2,I+MMM)=-1.0
            PT(1,I+MMM)=0.0
            PT(2,I+MMM)=0.0
    2 CONTINUE
C
      NNN=0
      DO 3 I=1,4
         NNN=NNN+NN(I)
    3 CONTINUE
C
      NEP1=2*NNN
C
      K=0
      DO 4 I=1,4
         NNI=NN(I)
      DO 4 J=1,NNI
C  ----------------------->  îOPMàPOBAHàE HAÉPìáKà CãEBA
         K=K+1
         L=2*K-1
         M=(I-1)*3+1
         LTEL(NEP+L)=NL(M,J)
         LST(NEP+L)=7
         PP(1,NEP+L)=1.0
         PP(2,NEP+L)=1.0
         PT(1,NEP+L)=0.0
         PT(2,NEP+L)=0.0
C  ----------------------->  îOPMàPOBAHàE HAÉPìáKà CèPABA
            L=2*K
            M=3*I
            LTEL(NEP+L)=NL(M,J)+MM(3)-1
            LST(NEP+L)=3
            PP(1,NEP+L)=-1.0
            PP(2,NEP+L)=-1.0
            PT(1,NEP+L)=0.0
            PT(2,NEP+L)=0.0
    4 CONTINUE
C
      RETURN
      END
      SUBROUTINE XREG( MM,NN,LM,LN,NF,X,Y,ISF,NSF,NUFD,NLD,ZONE,
     *                KLIN,ALFA,X0,Y0,R,IM,AKOEF,   NUF,NL,NEL,NJ)
C----------------------------------------------------------------------
C         èPOñEÑìPA CòàBAHàü KOHEóHOùãEMEHTHõX îPAÉMEHTOB
C----------------------------------------------------------------------
      COMMON
     */AL2/XYT(2,550)
     */AL2A/NC(8,110)
C
      DIMENSION NN(5),MM(4),X1(4),Y1(4),NU(30),NL1(10),
     *          X(4,12),Y(4,12),ISF(5,6),NSF(12),
     *          NUF(12,30),NUFD(12,30),NL(12,10),NLD(12,10)
C
      NJ=1
      NEL=1
      LM1=LM+2
      LN1=LN+2
      NN(LN1-1)=1
      MM(LM1-1)=1
C
      DO 10 I=2,LN1
         N1=NN(I-1)
         KK=1
         DO 20 K=1,N1
            DO 30 L=1,2
               DO 31 J=2,LM1
                  NOF =ISF(I,J)
                  NOF1=ISF(I,J-1)
                  NF1=0
                  IF(K.NE.1.OR.L.NE.1) GO TO 110
                  NOF2=ISF(I-1,J-1)
                  NOF3=ISF(I-1,J)
                  IF(NOF.NE.0.OR.NOF3.NE.0) GOTO 120
                  IF(NOF1.NE.0.OR.NOF2.NE.0) NF1=1
                  GOTO 31
C
 120              IF(NOF.NE.0)  NF1=NUFD(NOF,KK)
                  MIM=1
                  IF(NOF3.NE.0.AND.NSF(NOF3).LT.0) MIM=0
                  IF(I.NE.2) N2=2*NN(I-2)+MIM
                  IF(NOF.EQ.0.AND.NOF3.NE.0) NF1=NUFD(NOF3,N2)
                  IF(NOF.NE.0) NUF(NOF,1)=NJ
                  IF(NOF3.NE.0) NUF(NOF3,N2)=NJ
                  GOTO 31
C
 110              IF(NOF.EQ.0) GOTO 32
                  IF(L.EQ.1) GOTO 40
                  NL(NOF,K)=NEL
                  NEL=NEL+NLD(NOF,K)
                  IF(NSF(NOF).LT.0) GO TO 31
                  NF1=NUFD(NOF,KK)
                  NUF(NOF,2*K)=NJ
                  GOTO 31
C
  40              NUF(NOF,2*K-1)=NJ
                  NF1=NUFD(NOF,KK)
                  GOTO 31
C
  32              IF(NOF1.NE.0.AND.NSF(NOF1).GT.0 ) NF1=1
  31           NJ=NJ+NF1
               KK=KK+1
  30        CONTINUE
  20     CONTINUE
  10  CONTINUE
C
      NEL=NEL-1
      NJ=NJ-1
C
      DO 60 I=2,LN1
         N1=NN(I-1)
         N2=N1*2+1
      DO 60 J=2,LM1
         M1=MM(J-1)
         NOF=ISF(I,J)
         IF(NOF.EQ.0) GOTO 60
         IF(NSF(NOF).LT.0) N2=2
         DO 80 L=1,N2
  80        NU(L)=NUF(NOF,L)
            DO 90 M=1,N1
  90        NL1(M)=NL(NOF,M)
         DO 70 K=1,4
            X1(K)=X(K,NOF)
  70        Y1(K)=Y(K,NOF)
C
         KNS=IABS(NSF(NOF))
         GOTO (210,220,230,240,250),KNS
C
 210     CALL XELR1(M1,N1,X1,Y1,NU,NL1)
         GOTO 60
 220     IF(IM.EQ.0) CALL XELR2(NOF,ZONE,ZONE,X,Y,NL,NUF,KLIN,ALFA)
         IF(IM.EQ.1) CALL XELR21(NOF,ZONE,ZONE,X,Y,NL,NUF,KLIN,ALFA)
         GOTO 60
 230     IF(IM.EQ.0) CALL XELR3(NOF,ZONE,ZONE,X,Y,NL,NUF,KLIN,ALFA)
         IF(IM.EQ.1) CALL XELR31(NOF,ZONE,ZONE,X,Y,NL,NUF,KLIN,ALFA)
         GOTO 60
 240     IF(IM.EQ.1) CALL XELR4(NOF,X0,Y0,ZONE,R,AKOEF,X,Y,NL,NUF)
         IF(IM.EQ.0) CALL XELR41(NOF,X0,Y0,ZONE,R,AKOEF,X,Y,NL,NUF)
         GOTO 60
 250     IF(IM.EQ.1) CALL XELR5(NOF,X0,Y0,ZONE,R,AKOEF,X,Y,NL,NUF)
         IF(IM.EQ.0) CALL XELR51(NOF,X0,Y0,ZONE,R,AKOEF,X,Y,NL,NUF)
  60     CONTINUE
C
      RETURN
      END
      SUBROUTINE XRING(AKOEF,X0,Y0,R1,R2,R3,R4,NOM,ANGL,IQUADR,I,IPR)
C
      COMMON /AL2/XYT(2,550)
C
      IF(IQUADR.NE.1) GO TO 1
      ANGL1=ANGL/180.*3.14159
      CO=COS(ANGL1)*AKOEF
      SI=SIN(ANGL1)
      R11=R4
      R22=R3
      R33=R2
      R44=R1
      IF(IPR.EQ.0) R11=R3
      IF(IPR.EQ.0) R22=R2
      IF(IPR.EQ.0) R33=R1
      XYT(1,NOM)=X0+R11*CO
      XYT(1,NOM+1*I)=X0+R22*CO
      XYT(1,NOM+2*I)=X0+R33*CO
      IF(IPR.NE.0) XYT(1,NOM+3*I)=X0+R44*CO
      XYT(2,NOM)=Y0+R11*SI
      XYT(2,NOM+1*I)=Y0+R22*SI
      XYT(2,NOM+2*I)=Y0+R33*SI
      IF(IPR.NE.0) XYT(2,NOM+3*I)=Y0+R44*SI
      GO TO 4
C
    1 IF(IQUADR.NE.2) GO TO 2
      ANGL1=(180.-ANGL)/180.*3.14159
      CO=COS(ANGL1)*AKOEF
      SI=SIN(ANGL1)
      XYT(1,NOM)=X0-R1*CO
      XYT(1,NOM+1*I)=X0-R2*CO
      XYT(1,NOM+2*I)=X0-R3*CO
      IF(IPR.NE.0) XYT(1,NOM+3*I)=X0-R4*CO
      XYT(2,NOM)=Y0+R1*SI
      XYT(2,NOM+1*I)=Y0+R2*SI
      XYT(2,NOM+2*I)=Y0+R3*SI
      IF(IPR.NE.0) XYT(2,NOM+3*I)=Y0+R4*SI
      GO TO 4
C
    2 IF(IQUADR.NE.3) GO TO 3
	ANGL1=(ANGL-180.0)/180.0*3.14159
      CO=COS(ANGL1)*AKOEF
      SI=SIN(ANGL1)
      XYT(1,NOM)=X0-R1*CO
      XYT(1,NOM+1*I)=X0-R2*CO
      XYT(1,NOM+2*I)=X0-R3*CO
      IF(IPR.NE.0) XYT(1,NOM+3*I)=X0-R4*CO
      XYT(2,NOM)=Y0-R1*SI
      XYT(2,NOM+1*I)=Y0-R2*SI
      XYT(2,NOM+2*I)=Y0-R3*SI
      IF(IPR.NE.0) XYT(2,NOM+3*I)=Y0-R4*SI
      GO TO 4
C
    3 IF(IQUADR.NE.4) GO TO 4
      ANGL1=(360.-ANGL)/180.*3.14159
      CO=COS(ANGL1)*AKOEF
      SI=SIN(ANGL1)
      R11=R4
      R22=R3
      R33=R2
      R44=R1
      IF(IPR.EQ.0) R11=R3
      IF(IPR.EQ.0) R22=R2
      IF(IPR.EQ.0) R33=R1
      XYT(1,NOM)=X0+R11*CO
      XYT(1,NOM+1*I)=X0+R22*CO
      XYT(1,NOM+2*I)=X0+R33*CO
      IF(IPR.NE.0) XYT(1,NOM+3*I)=X0+R44*CO
      XYT(2,NOM)=Y0-R11*SI
      XYT(2,NOM+1*I)=Y0-R22*SI
      XYT(2,NOM+2*I)=Y0-R33*SI
      IF(IPR.NE.0) XYT(2,NOM+3*I)=Y0-R44*SI
C
    4 RETURN
      END
      SUBROUTINE XTOP(IM,NL)
C----------------------------------------------------------------------
C     CMEôEHàE ìáãOB Kù HA 1/4 CTOPOHõ B OKPECTHOCTà BEPòàHõ TPEôàHõ
C----------------------------------------------------------------------
      COMMON
     */AL2/XYT(2,550)
     */AL2A/NC(8,110)
C
      DIMENSION NL(12,10)
C
      N1=NL(5,1)
      N2=NL(8,1)
C
      IF(IM.EQ.1) GO TO 10
      M1=NC(5,N1+3)
      M2=NC(8,N1+3)
      M3=NC(4,N1+3)
      M4=NC(4,N1+7)
      M5=NC(4,N1+11)
      M6=NC(4,N1+15)
      M7=NC(4,N2+15)
      M8=NC(4,N2+11)
      M9=NC(4,N2+7)
      M10=NC(4,N2+3)
C
   10 IF(IM.EQ.0) GO TO 20
      M1=NC(5,N1+4)
      M2=NC(8,N1+4)
      M3=NC(4,N1+4)
      M4=NC(4,N1+9)
      M5=NC(4,N1+14)
      M6=NC(4,N1+19)
      M7=NC(4,N2+19)
      M8=NC(4,N2+14)
      M9=NC(4,N2+9)
      M10=NC(4,N2+4)
   20 CONTINUE
C
      DO 30 I=1,2
      XYT(I,M2)=(XYT(I,M2)+XYT(I,M1))/2.
      XYT(I,M3)=(XYT(I,M3)+XYT(I,M1))/2.
      XYT(I,M4)=(XYT(I,M4)+XYT(I,M1))/2.
      XYT(I,M5)=(XYT(I,M5)+XYT(I,M1))/2.
      XYT(I,M6)=(XYT(I,M6)+XYT(I,M1))/2.
      XYT(I,M7)=(XYT(I,M7)+XYT(I,M1))/2.
      XYT(I,M8)=(XYT(I,M8)+XYT(I,M1))/2.
      XYT(I,M9)=(XYT(I,M9)+XYT(I,M1))/2.
      XYT(I,M10)=(XYT(I,M10)+XYT(I,M1))/2.
   30 CONTINUE
C
      RETURN
      END
