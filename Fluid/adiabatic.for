      program adiabatic
        implicit none
        integer,parameter::MaxIter=5
        real,parameter::Epsilon=1e-5
        real,parameter::k=1.4
        real A,B,C
        real P,hiL,ksi,x,D
        real l1,l2!,lcr
        
        A=(k-1.0)/(k+1.0)
        B=2.0/(k+1.0)
        C=B**(-1.0/(k-1.0))
        
        P=2.0
        ksi=0.01
        x=200.0
        D=1.0
        
        hiL=2.0*k*ksi*x/(D*(k+1.0))
!        lcr=lcritical()
        
        write(*,'(''P  ='',f12.4)')P
        write(*,'(''hiL='',f12.4)')hiL
!        write(*,'(''lcr='',f12.4)')lcr
        
        call nonlinsimq(l1,l2)
        write(*,'(''l1 ='',f12.4)')l1
        write(*,'(''l2 ='',f12.4)')l2
        
      contains
        
        real function tau(l)
          implicit none
          real l
          tau=1.0-A*l**2
        end function
        
        real function pie(l)
          implicit none
          real l
          pie=tau(l)**(k/(k-1.0))
        end function
        
        real function eps(l)
          implicit none
          real l
          eps=tau(l)**(1.0/(k-1.0))
        end function
        
        real function q(l)
          implicit none
          real l
          q=C*eps(l)*l
        end function
        
        real function phi(l)
          implicit none
          real l
          phi=1.0/l**2+2.0*log(l)
        end function
        
        real function dphi(l)
          implicit none
          real l
          dphi=-2.0/l**3+2.0/l
        end function
        
        real function y(l)
          implicit none
          real l
          y=q(l)/pie(l)
        end function
        
        real function lcritical()
          implicit none
          integer i
          real xn,oldxn
          
          xn=1.0
          i=1
          do while (.true.)
            oldxn=xn
            xn=oldxn-(phi(oldxn)-(1.0+hiL))/dphi(oldxn)
            if ((abs(xn-oldxn)<Epsilon).or.(i>MaxIter)) exit
            i=i+1
          end do
          
          lcritical=xn
        end function
      
        real function F1(l1,l2)
          implicit none
          real l1,l2
          F1=phi(l1)-phi(l2)-hiL
        end function
        
        real function F2(l1,l2)
          implicit none
          real l1,l2
          F2=y(l2)-P*q(l1)
        end function
        
        real function dF1l1(l1,l2)
          implicit none
          real l1,l2
          l2=l2
          dF1l1=dphi(l1)
        end function
        
        real function dF1l2(l1,l2)
          implicit none
          real l1,l2
          l1=l1
          dF1l2=-dphi(l2)
        end function
        
        real function dF2l1(l1,l2)
          implicit none
          real l1,l2
          l2=l2
          dF2l1=-P*B*eps(l1)*(1.0-2.0*B*l1**2/tau(l1))
        end function
        
        real function dF2l2(l1,l2)
          implicit none
          real l1,l2
          l1=l1
          dF2l2=C*(1.0+2.0*A*l2**2)/tau(l2)
        end function
        
        subroutine nonlinsimq(x,y)
          implicit none
          real x,y
          real,dimension(2,2)::J
          real,dimension(2)::dx,F
          integer i
          
          x=0.1
          y=0.9
          dx(1)=0.1
          dx(2)=0.1
          
          i=0
          do while(.true.)
            J(1,1)=dF1l1(x,y)
            J(1,2)=dF1l2(x,y)
            J(2,1)=dF2l1(x,y)
            J(2,2)=dF2l2(x,y)
            F(1)=-F1(x,y)
            F(2)=-F2(x,y)
            call simq2(J,F,dx)
            x=x+dx(1)
            y=y+dx(2)
            !write(*,*)x,y
            if (sqrt(dx(1)**2+dx(2)**2)<Epsilon.or.i>MaxIter) exit
            i=i+1
          end do
        end subroutine
      
        subroutine simq2(a,b,x)
          implicit none
          real,dimension(2,2)::a
          real,dimension(2)::b,x
          real delta
          !k=a(1,1)/a(2,1)
          !x(2)=(b(1)-b(2)*k)/(a(1,2)-k*a(2,2))
          !x(1)=(b(1)-a(1,2)*x(2))/a(1,1)
          delta=a(1,1)*a(2,2)-a(2,1)*a(1,2)
          x(1)=(b(1)*a(2,2)-b(2)*a(1,2))/delta
          x(2)=(b(2)*a(1,1)-b(1)*a(2,1))/delta
        end subroutine
              
      end program
      
