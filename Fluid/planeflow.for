      program planeflow
        implicit none
        real,parameter::Eps=1e-5
        real,parameter::EpsLin=(2.+.1*10)/100.
        real,parameter::Pi=3.141592653598
        real,parameter::k=1.4
        integer,parameter::fout=100
        integer,parameter::N=100
        integer,parameter::MaxIter=10000
        integer i
        real alpha,v,Pdown,Pup
        real,dimension(N)::M,Cylin,Cx,Cy,attack
        real dM,dalpha
        real Cylin0,Cx0,Cy0
        
        v=sqrt((k+1.)/(k-1.))
        dM=.1
        dalpha=0.01*Pi/180.
        
        do i=1,N
          M(i)=1.+dM*i
          
          Cylin0=0.
          Cy0=0.
          Cx0=0.
          alpha=0.
          
          do while (.true.)
            Cylin(i)=Cylin0
            Cy(i)=Cy0
            Cx(i)=Cx0
            attack(i)=alpha*180./Pi
            alpha=alpha+dalpha
            
            Pdown=2.*(termination(M(i),alpha)-1.)/(k*M(i)**2)
            Pup=2.*(wave(M(i),alpha)-1.)/(k*M(i)**2)
            
            Cylin0=4.*alpha/sqrt(M(i)**2-1.)
            Cy0=(Pdown-Pup)*cos(alpha)
            Cx0=(Pdown-Pup)*sin(alpha)
            
            if (abs((Cy0-Cylin0)/Cy0)>=EpsLin) exit
          end do
        end do
        
        ! Results output
        open(unit=fout,file='planeflow.csv')
        write(fout,'(12x,''M;'',7x,''attack;'',
     *    8x,''Cylin;'',11x,''Cy;'',11x,''Cx;'')')
        do i=1,N
          write(fout,'(5(f13.6,'';''))')M(i),attack(i),
     *    Cylin(i),Cy(i),Cx(i)
        end do
        close(fout)
        
      contains
      
        real function ftau(M)
          implicit none
          real M
          ftau=1./(1.+.5*(k-1.)*M**2)
        end function
                
        real function fpi(M)
          implicit none
          real M
          fpi=ftau(M)**(k/(k-1.))
        end function
        
        real function func(M)
          implicit none
          real M,b
          b=sqrt(M**2-1.)
          func=v*atan(b/v)-atan(b)
        end function
        
        real function dfunc(M)
          implicit none
          real M,b
          b=sqrt(M**2-1.)
          dfunc=(v**2*M/(M**2+v**2-1.)-1./M)/b
        end function
        
        real function termination(M1,theta)
          implicit none
          real M1,theta
          real m,a,b,c,x,oldx
          real f,df
          integer i
          
          m=1./(1.+.5*(k-1.)*M1**2)
          a=(1.-M1**2)*m/tan(theta)
          b=(1.+.5*(k+1.)*M1**2)*m
          c=m/tan(theta)
          
          x=atan(1.)/sqrt(M1**2+1.)
          oldx=10.+x
          i=0
          
          do while (abs(x-oldx)>=Eps.and.i<MaxIter)
            oldx=x
            f=oldx**3+a*oldx**2+b*oldx+c
            df=3.*oldx**2+2.*a*oldx+b
            x=oldx-f/df
            i=i+1
          end do
          x=atan(x)
          termination=(2.*k*((M1*sin(x))**2)-(k-1.))/(k+1.)
        end function
        
        real function wave(M1,theta)
          implicit none
          real M1,theta
          real tstar,tt
          real Mi1,Mi
          integer i
          
          tstar=func(M1)
          tt=theta+tstar
          
          Mi1=M1
          i=0
          do while (abs(Mi-Mi1)>=Eps.and.i<MaxIter)
            Mi=Mi1
            Mi1=Mi-(func(Mi)-tt)/dfunc(Mi)
            i=i+1
          end do
          wave=fpi(Mi1)/fpi(M1)
        end function
        
      end program
