      SUBROUTINE polovin(f,x,a,b)
      IMPLICIT REAL*8 (a-h,o-z)
      xa=a
      xb=b
      fa=f(xa)
      fb=f(xb)
      IF (fa*fb.GT.0.0d0)THEN
        WRITE(*,*)' ���� ����� �� ����'
        RETURN
      END IF
      DO i=1,1000
        xc=0.5d0*(xa+xb)
        fc=f(xc)
        IF (ABS(xa-xc).LT.1.0d-8)THEN
          x=xc
          RETURN
        END IF
        IF (fc*fa.GT.0.0d0)THEN
          xa=xc
          fa=fc
        ELSE
          xb=xc
          fb=fc
        END IF
      END DO
      x=xc
      WRITE(*,*)' �������� �筮��� �� ���⨣���'
      RETURN
      END
