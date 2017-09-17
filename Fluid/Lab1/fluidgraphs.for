      program fluidgraphcs
        implicit none
        real,parameter::k=1.4
        integer,parameter::fout=100
        real lambda,tau_,pi_,eps_,q_
        
        open(unit=fout,file='fluidgraphs.csv')
        
        lambda=0.1
        do while (lambda<2.4)
          tau_=tau(lambda)
          pi_=pi(lambda)
          eps_=eps(lambda)
          q_=q(lambda)
          write(fout,'(5(f12.4,'';''))'),lambda,tau_,pi_,eps_,q_
          lambda=lambda+0.1
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
        
        real function q(lambda)
          implicit none
          real lambda,q1
          q1=(2/(k+1))**(-1/(k-1))
          q=lambda*eps(lambda)*q1
        end function
      end program
