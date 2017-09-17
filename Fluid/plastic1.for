      program plastic1
        implicit none
        integer,parameter::fout=100
        real,parameter::Delta=1e-3
        real,parameter::Pi=3.1415926
        integer,parameter::n=5
        integer,parameter::nD=1000
        real,parameter::omega=1.0
        real,parameter::Ex=1.0
        real,parameter::Hx=0.1
        real,parameter::Kx=1.0
        real,parameter::ht=2*n*Pi/(omega*nD)
        real,dimension(0:nD)::t,sigma,eps,alpha,epsp,od
        integer i
        
        t(0)=0.0
        sigma(0)=0.0
        eps(0)=0.0
        alpha(0)=0.0
        epsp(0)=0.0
        od(0)=0.0
        
        do i=1,nD
          t(i)=ht*i
          sigma(i)=sig(t(i))
          eps(i)=eps(i-1)+(sigma(i)-sigma(i-1))/Et(od(i-1))
          epsp(i)=eps(i)-sigma(i)/Ex
          alpha(i)=alpha(i-1)+Kx*(epsp(i)-epsp(i-1))
          od(i)=od(i)+abs(epsp(i)-epsp(i-1))
          if(F(sigma(i),alpha(i),od(i))<0.0.or.
     *      (sigma(i)-alpha(i))*(epsp(i)-epsp(i-1))<=0.0)then
            eps(i)=eps(i-1)+(sigma(i)-sigma(i-1))/Ex
            epsp(i)=epsp(i-1)
            alpha(i)=alpha(i-1)
            od(i)=od(i-1)
          end if
        end do
        
        open(unit=fout,file='plastic1.csv')
        do i=0,nD
          write(fout,'(6(f12.4,'';''))'),t(i),sigma(i),eps(i),epsp(i),
     *      alpha(i),od(i)
        end do
        close(fout)
        
      contains
      
        real function sig(tx)
          implicit none
          real tx
          sig=2.0*sin(omega*tx)
        end function
        
        real function sigmaY(odx)
          implicit none
          real odx
          sigmaY=Hx*odx+1.0
        end function
        
        real function F(sigmax,alphax,odx)
          implicit none
          real sigmax,alphax,odx
          F=abs(sigmax-alphax)-sigmaY(odx)
        end function
        
        real function Et(odx)
          implicit none
          real odx
          Et=Ex-Ex**2/(Ex+Kx+
     *      (sigmaY(odx+Delta)-sigmaY(odx))/Delta)
        end function
      
      end program
      
