USE Material
USE Plit
USE Modl

IMPLICIT NONE


REAL(8), EXTERNAL :: deter, Time
REAL(8) s(3), e(3), ddd(3,3), Pi, SQ2, SQ3, h_a, x
REAL(8), ALLOCATABLE :: tim(:), a0(:,:)
INTEGER ncalk, k 
INTEGER, ALLOCATABLE :: la(:)

OPEN(1,FILE='trek.dat',STATUS='OLD')
OPEN(2,FILE='TrekDeviator.dat')
READ(1,*)aa,an,rr,ak0,eta,sig_s,eps_s
WRITE(2,'(7(g15.10,1x))')aa,an,rr,ak0,eta,sig_s,eps_s

READ(1,*)K2
WRITE(2,*)k2
ALLOCATE(a0(3,k2),la(k2),tim(k2))
DO j=1,k2
	READ(1,*)A0(:,J),LA(J), Tim(j)
	s0=(A0(1,J)+A0(2,J)+A0(3,J))/3.0D0
	WRITE(2,'(3(g15.5,1x),i4,g15.5)') (A0(1,J)-A0(2,J))/SQ2, (A0(3,J)-s0)*SQ3/SQ2, 0.0D0, LA(J), Tim(j)
END DO
CLOSE(1)
CLOSE(2)

!RUNQQ (MDF_Project TrekDeviator.dat))

OPEN(1,FILE='plit.res')
OPEN(2,FILE='plit.dat',STATUS='OLD')
OPEN(3,FILE='trek.dat',STATUS='OLD')
OPEN(4,FILE='text.res')
OPEN(5,FILE='Matrix.res',STATUS='OLD')
pi=ATAN(1.0D0)*4.0D0
sq2=SQRT(2.0d0)
sq3=SQRT(3.0d0)
!
!     Plate geometrical parameters
!
READ(2,*)a_b
READ(2,*)mbig,mend
READ(2,*)nbig,nend
WRITE(4,'(2(1x,a,i2,a,i2/))')' Wave parameters m from ',mbig,' to ',mend,'                 n from ',nbig,' to ',nend
CLOSE(2)
!
!     Model constants
!


READ(3,*)aa,an,rr,ak0,eta,sig_s,eps_s,ncalk

t0=sig_s*sq2/sq3
es=eps_s*sq2/sq3*(1.0D0+0.3D0)

gupr=t0/es/2.0d0
lam=2.0d0*gupr*0.3d0/(1.0d0-2.0d0*0.3d0)
kvol=(lam+2.d0/3.d0*gupr)

k=0
DO WHILE (NOT(EOF(5)))
	k=k+1
	READ(5,*)s,e,ddd
	IF (k.GE.ncalk)THEN
		sig1=((s(1)-SQ3*s(2))/SQ2)/gupr
		sig2=-((SQ3*S(2)+S(1))/SQ2)/gupr
		g1212=(ddd(3,3))/2.0
		g1133=(kvol+ddd(2,1)/sq3-ddd(2,2)/3.d0)
		g2233=(kvol-ddd(2,1)/sq3-ddd(2,2)/3.d0)
		g3333=(kvol+ddd(2,2)/3.d0*2.d0)
		g1111=(3.d0*kvol-g1133+ddd(1,1)-ddd(1,2)/sq3)/2.d0
		g1122=(3.d0*kvol-g2233-ddd(1,1)-ddd(1,2)/sq3)/2.d0
		g2222=(g1122+ddd(1,1)+ddd(1,2)/sq3)

		g1212=g1212/gupr
		g1111=g1111/gupr
		g1122=g1122/gupr
		g1133=g1133/gupr
		g2222=g2222/gupr
		g2233=g2233/gupr
		g3333=g3333/gupr
		h_a=-1.0

		DO mpar=mbig,mend
			DO npar=nbig,nend
				CALL polovin(deter,x,0.00001d0,0.5d0)
				IF(x.GT.h_a)THEN
					h_a=x
					mres=mpar
					nres=npar
				END IF
			END DO
		END DO
		WRITE(1,'(i4,1x,2(f10.5,1x),f15.10,2(i2,1x))')k,sig1*gupr/sig_s,sig2*gupr/sig_s,h_a, nres, mres
	END IF
END DO
END

REAL(8) FUNCTION TIME()
IMPLICIT REAL*8(A-H,O-Z)
INTEGER*2 h,m,s,s100
CALL GETTIM(h,m,s,s100)
time=h*3600.+m*60.+s+s100/100.
RETURN
END
