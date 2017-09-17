      program flowturn
        implicit none
        real,parameter::Eps=1e-3
        real,parameter::MaxIter=50000
        real,parameter::Pi=3.1415926
        real,parameter::k=1.4
        real,parameter::M1=3.5
        real,parameter::theta1=20.0
        real t,M2
        
        t=theta1*Pi/180
        M2=M2func(M1,t)
        
        print *, 'theta = ', theta1
        print *, 'M1 = ', M1
        print *
        print *, 'M2 = ', M2
        
      contains
      
        real function func(M)
          implicit none
          real M,a,b
          a=sqrt((k+1)/(k-1))
          b=sqrt(M**2-1)
          func=a*atan(b/a)-atan(b)
        end function
        
        real function dfunc(M)
          implicit none
          real M,a,b
          a=sqrt((k+1)/(k-1))
          b=sqrt(M**2-1)
          dfunc=b*(a**2-1)/(M*(b**2+a**2))
        end function
        
        real function M2func(M,theta)
          implicit none
          real M,theta
          real tstar,tt
          real Mi1,Mi
          integer i
          
          tstar=func(M)
          tt=theta+tstar
        
          Mi1=M
          i=0
          do while (.true.)
            Mi=Mi1
            Mi1=Mi-(func(Mi)-tt)/dfunc(Mi)
            if ((abs(Mi1-Mi)<Eps).or.(i>MaxIter))then
              exit
            end if
            i=i+1
          end do
          M2func = Mi1
        end function
      
      end program
      
