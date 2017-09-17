      program fluidcritical
        implicit none
        integer,parameter::fout=100
        real k,dk,lambda,tau_,pi_,eps_
        
        lambda=1
        k=1.1
        dk=(1.67-1.1)/25
        
        open(unit=fout,file='fluidcritical.csv')
        
        do while (k<1.67)
          tau_=tau(lambda)
          pi_=pi(lambda)
          eps_=eps(lambda)
          write(fout,'(4(f12.4,'';''))'),k,tau_,pi_,eps_
          k=k+dk
        end do
        
        close(fout)
        
      contains
        real function tau(lambda)
          implicit none
          real lambda
          tau=1-lambda**2*(k-1)/(k+1)
        end function
        
        real function pi(lambda)
          implicit none
          real lambda
          pi=tau(lambda)**(k/(k-1))
        end function
        
        real function eps(lambda)
          implicit none
          real lambda
          eps=tau(lambda)**(1/(k-1))
        end function
      end program
      
