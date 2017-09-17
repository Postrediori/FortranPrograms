      program termination
        implicit none
        real,parameter::eps=1e-5
        integer,parameter::MaxIter=50000
        real,parameter::pi=3.14159265358979323846
        real,parameter::M1=2.5
        real,parameter::theta1=11.20
        real,parameter::k=1.4
        
        real beta_
        real a,b
        real sin2bmax,bmax,tantmax,tmax
        real t
        
        t = theta1*Pi/180
        
        a=(k+1)*M1**2/4-1
        b=(k+1)*(1+(k-1)*M1**2/2+(k+1)*M1**4/16)
        sin2bmax=(a+sqrt(b))/(k*M1**2)
        bmax=asin(sqrt(sin2bmax))
        
        a=((M1*sin2bmax)**2-1)
        b=1+M1**2*((k+1)/2-sin2bmax)
        tantmax=a/(tan(bmax)*b)
        tmax=atan(tantmax)
        
        print *, 'M = ', M1
        print *, 'theta1 = ', theta1
        print *
        print *, 'tmax = ', tmax*180/Pi
        
        if (t>tmax) then
          print *, 'Angle is greater than maximal'
        else
          beta_=funcbeta(M1,t)
          print *,'beta = ',beta_*180/pi
        end if
        
      contains
                   
        real function funcbeta(M,theta)
          implicit none
          real M,theta
          real mi,a,b,c
          real func,dfunc
          real xn,oldxn
          integer i
          
          mi=1/(1+(k-1)*M**2/2)
          a=(1-M**2)*mi/tan(theta)
          b=(1+(k+1)*M**2/2)*mi
          c=mi/tan(theta)
          
          xn=atan(1.0) / sqrt(M**2+1)
          
          oldxn=0
          i=0
          do while (.true.)
            oldxn=xn
            func=xn**3+a*xn**2+b*xn+c
            dfunc=3*xn**2+2*a*xn+b
            xn=oldxn-func/dfunc
            if ((abs(xn-oldxn)<eps).or.(i>MaxIter)) exit
            i=i+1
          end do
          
          funcbeta=atan(xn)
        end function
        
      end program
