REAL(8) FUNCTION deter(H_A)
USE PLIT
USE MODL

IMPLICIT NONE

REAL(8) h_a

COMPLEX const(3), root

REAL(8) bound(6,6), d(3,4), Pi, a, b, ak, bk, ck, dk, x1, x2, x3, det
INTEGER le(6), me(6), kod

pi=DATAN(1.0D0)*4.0D0
a=pi*mpar
b=pi*npar*a_b
!cÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÂÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÂÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
!c³     d(1,1)*lam*lam-d(1,2) ³d(1,3)                ³d(1,4)*lam            ³
!c³     d(2,1)                ³d(2,2)*lam*lam-d(2,3) ³d(2,4)*lam            ³
!c³     d(3,1)*lam            ³d(3,2)*lam            ³d(3,3)*lam*lam-d(3,4) ³
!cÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

d(1,1)=g1212-0.5*sig1
d(2,2)=g1212-0.5*sig2
d(3,3)=g3333
d(1,2)=(g1111-sig1)*a*a+(g1212+0.5*(sig2-sig1))*b*b
d(2,3)=(g2222-sig2)*b*b+(g1212+0.5*(sig1-sig2))*a*a
d(1,3)=-(g1122+g1212-0.5*(sig1+sig2))*a*b
d(1,4)=(g1133+g1212-0.5*sig1)*a
d(2,1)=d(1,3)
d(2,4)=(g2233+g1212-0.5*sig2)*b
d(3,1)=-(g1133+g1212-0.5*sig1)*a
d(3,2)=-(g2233+g1212-0.5*sig2)*b
d(3,4)=(g1212+0.5*sig1)*a*a+(g1212+0.5*sig2)*b*b
ak=d(1,1)*d(2,2)*d(3,3)
bk=-d(1,1)*d(2,2)*d(3,4)-d(1,1)*d(2,3)*d(3,3)-d(1,2)*d(2,2)*d(3,3)-d(2,2)*d(3,1)*d(1,4)-d(1,1)*d(3,2)*d(2,4)
bk=bk/ak
ck=d(1,2)*d(2,3)*d(3,3)+d(1,2)*d(2,2)*d(3,4)+d(1,1)*d(2,3)*d(3,4)+d(2,1)*d(3,2)*d(1,4)+d(1,3)*d(2,4)*d(3,1)+d(1,4)*d(3,1)*d(2,3)-d(2,1)*d(1,3)*d(3,3)+d(3,2)*d(2,4)*d(1,2)
ck=ck/ak
dk=-d(1,2)*d(2,3)*d(3,4)+d(2,1)*d(1,3)*d(3,4)
dk=dk/ak
CALL cordano(bk,ck,dk,x1,x2,x3,kod)
IF(kod.EQ.4)THEN
root=CMPLX(SQRT(x1))
CALL s_vect(root,const,d)
bound(1,1)=REAL(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,1)=REAL(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,1)=REAL(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,1)=REAL(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,1)=REAL(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,1)=REAL(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))

root=-root
CALL s_vect(root,const,d)
bound(1,2)=REAL(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,2)=REAL(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,2)=REAL(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,2)=REAL(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,2)=REAL(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,2)=REAL(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))

root=CSQRT(CMPLX(x2,x3))
CALL s_vect(root,const,d)
bound(1,3)=REAL(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,3)=REAL(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,3)=REAL(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,3)=REAL(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,3)=REAL(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,3)=REAL(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))

bound(1,4)=AIMAG(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,4)=AIMAG(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,4)=AIMAG(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,4)=AIMAG(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,4)=AIMAG(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,4)=AIMAG(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))

root=-root
CALL s_vect(root,const,d)
bound(1,5)=REAL(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,5)=REAL(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,5)=REAL(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,5)=REAL(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,5)=REAL(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,5)=REAL(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(1,6)=AIMAG(ZEXP( root*h_a)*(a*const(3)+root*const(1)))
bound(2,6)=AIMAG(ZEXP( root*h_a)*(b*const(3)+root*const(2)))
bound(3,6)=AIMAG(ZEXP( root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
bound(4,6)=AIMAG(ZEXP(-root*h_a)*(a*const(3)+root*const(1)))
bound(5,6)=AIMAG(ZEXP(-root*h_a)*(b*const(3)+root*const(2)))
bound(6,6)=AIMAG(ZEXP(-root*h_a)*(-a*const(1)*g1133-b*const(2)*g2233+root*const(3)*g3333))
CALL minv(bound,6,det,LE,ME)
deter=det
RETURN
END IF
RETURN
END

SUBROUTINE s_vect(lam,a,d)
!cÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
!c³      å®¤¨â á®¡áâ¢¥­­ë© ¢¥ªâ®à const á®®â¢¥âáâ¢ãîé¨© ª®à­î root    ³
!cÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

IMPLICIT NONE
COMPLEX lam, m(3,3), dm(3,3), a(3)
REAL(8) d(3,4), dl
INTEGER i

m(1,1)=d(1,1)*lam*lam-d(1,2)
m(1,2)=d(1,3)
m(1,3)=d(1,4)*lam
m(2,1)=d(2,1)
m(2,2)=d(2,2)*lam*lam-d(2,3)
m(2,3)=d(2,4)*lam
m(3,1)=d(3,1)*lam
m(3,2)=d(3,2)*lam
m(3,3)=d(3,3)*lam*lam-d(3,4)

dm(1,1) =  ( m(2,2) * m(3,3) - m(3,2) * m(2,3) )
dm(1,2) = -( m(2,1) * m(3,3) - m(3,1) * m(2,3) )
dm(1,3) =  ( m(2,1) * m(3,2) - m(3,1) * m(2,2) )
dm(2,1) = -( m(1,2) * m(3,3) - m(3,2) * m(1,3) )
dm(2,2) =  ( m(1,1) * m(3,3) - m(3,1) * m(1,3) )
dm(2,3) = -( m(1,1) * m(3,2) - m(3,1) * m(1,2) )
dm(3,1) =  ( m(1,2) * m(2,3) - m(2,2) * m(1,3) )
dm(3,2) = -( m(1,1) * m(2,3) - m(2,1) * m(1,3) )
dm(3,3) =  ( m(1,1) * m(2,2) - m(2,1) * m(1,2) )

a(1)=0.
a(2)=0.
a(3)=0.

DO i=1,3
	a(1)=a(1)+dm(i,1)
	a(2)=a(2)+dm(i,2)
	a(3)=a(3)+dm(i,3)
END DO
dl=MAX(ABS(a(1)),ABS(a(2)),ABS(a(3)))
DO i=1,3
	a(i)=a(i)/dl
END DO
RETURN
END
