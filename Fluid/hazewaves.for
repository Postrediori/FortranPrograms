      program hazewaves
        implicit none
        
        integer,parameter::N=10
        integer,parameter::M=100
        real,parameter::L=10.0
        real,parameter::sigma=0.1
        real,parameter::a=1.0
        real,parameter::rho=1.0
        real,parameter::rho_a=rho*a
        
        integer i,j
        real h,h2
        real t,tau
        real,dimension(0:N)::x
        real,dimension(N)::u1,u2,p1,p2
        real,dimension(0:N)::U,P
        
        real,dimension(2,2)::S
        real,dimension(2)::Y
        real,dimension(2)::W
        
        real bl,cl
        real br,cr
        
        h=L/N
        h2=h/2.0
        tau=sigma*h/a
        
        bl=0.0
        cl=1.0
        br=0.0
        cr=1.0
        
        do i=0,N
          x(i)=h*i
        end do
        
        do i=1,N
          u1(i)=un(x(i-1)+h2)
          p1(i)=pn(x(i-1)+h2)
        end do
        
        do j=1,M
          t=tau*j
        
          S(1,1)=1.0
          S(1,2)=-1.0/rho_a
          S(2,1)=cl
          S(2,2)=bl
          Y(1)  =u1(1)-p1(1)/rho_a
          Y(2)  =0.0
          call simq2(S,Y,W)
          U(0)=W(1)!0.0
          P(0)=W(2)!p1(1)
        
          S(1,1)=1.0
          S(1,2)=1.0/rho_a
          S(2,1)=cr
          S(2,2)=br
          Y(1)  =u1(N)+p1(N)/rho_a
          Y(2)  =0.0
          call simq2(S,Y,W)
          U(N)=W(1)!0.0
          P(N)=W(2)!p1(N)
        
          do i=1,N-1
            U(i)=((u1(i+1)+u1(i))-(p1(i+1)-p1(i))/rho_a)/2.0
            P(i)=((p1(i+1)+p1(i))-(u1(i+1)-u1(i))*rho_a)/2.0
          end do
                    
          do i=1,N
            u2(i)=u1(i)-(P(i)-P(i-1))*tau/(h*rho)
            p2(i)=p1(i)-(U(i)-U(i-1))*rho_a*a*tau/h
          end do
          
          do i=1,N
            u1(i)=u2(i)
            p1(i)=p2(i)
          end do
          
          !do i=1,N
            !write(*,'(11(f7.1))')U(0:N)
            !write(*,'(11(f7.1))')P(0:N)
            write(*,'(''   '',10(f7.1))')u1(1:N)
            write(*,'(''   '',10(f7.1))')p1(1:N)
          !end do
          write(*,'(''----'')')
          
        end do
        
      contains
      
        real function un(x)
          implicit none
          real x
          x=x
          un=0.0
        end function
        
        real function pn(x)
          implicit none
          real x
          x=x
          if (x<L/2.0) then
            pn=2.0
          else
            pn=1.0
          end if
        end function
        
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
      
