SUBROUTINE cordano(b,c,d,x1,x2,x3,kod)
IMPLICIT NONE
COMPLEX calf(3)
REAL(8) b,c,d,x1,x2,x3,eps,p,q,det,zn,alf,be,al  
INTEGER kod
eps=1.0D-16
IF(ABS(b*b-3.0*c) .LT. eps .AND.ABS(b*b*b-27.0*d) .LT. eps)THEN
	kod=3
	x1=(-b)/3.0
	x2=x1
	x3=x1
	RETURN
END IF
p=c-b*b/3.0
q=2.0*b*b*b/27.0-b*c/3.D0+d
det=-(27.0*q*q+4.0*p*p*p)/108.D0
IF(ABS(det).LT.eps)THEN
	kod=2
	zn=q/ABS(q)
	alf=-zn*(ABS(q)/2.D0)**(1.D0/3.D0)
	x1=2.0*alf-b/3.D0
	x2=-alf-b/3.D0
	x3=x2
	RETURN
END IF
IF(det.GT.0.D0)THEN
	kod=1
	CALL rad3(CMPLX(-q/2.D0,SQRT(det)),calf)
	x1=REAL(-p/3.D0/calf(1)+calf(1))-b/3.D0
	x2=REAL(-p/3.D0/calf(2)+calf(2))-b/3.D0
	x3=REAL(-p/3.D0/calf(3)+calf(3))-b/3.D0
	RETURN
END IF
IF(det.LT.0.D0)THEN
	kod=4
	be=q/2.D0+SQRT(-det)
	zn=be/(ABS(be))
	be=-zn*(ABS(be))**(1.D0/3.D0)
	al=-p/3.D0/be
	x1=al+be-b/3.D0
	x2=-(al+be)/2.D0-b/3.D0
	x3=SQRT(3.0D0)*(al-be)/2.0
	RETURN
END IF
END


SUBROUTINE rad3(arg,res)
COMPLEX arg,res(3)
REAL(8) r,fi,pi,a,b
pi=ATAN(1.D0)*4.0D0
a=REAL(arg)
b=AIMAG(arg)
r=SQRT(a*a+b*b)
fi=ATAN(b/a)
DO i=1,3
	res(i)=CMPLX(r**(1.0/3.0)*COS((fi+2.0*pi*(i-1))/3.0),r**(1.0/3.0)*SIN((fi+2.0*pi*(i-1))/3.0))
END DO
RETURN
END


