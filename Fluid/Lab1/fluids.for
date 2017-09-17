      program fluids
        implicit none
        real,parameter::k=1.4
        real,parameter::eps=1e-3
        integer,parameter::maxiter=50000
        real qc,mnlambda,mxlambda

        mnlambda=0
        mxlambda=sqrt((k+1)/(k-1))
        print*,'__________q|_________l1|_________l2'

        ! q(lambda)-qc=0
        qc=0.1
        do while (qc<1.0)
          write(*,'(3(f12.4))'),qc,lk(qc,1),lk(qc,2)
          qc=qc+0.1
        end do

      contains

        real function q(lambda)
          implicit none
          real lambda,q1,q2
          q1=(1-lambda**2*(k-1)/(k+1))**(1/(k-1))
          q2=(2/(k+1))**(-1/(k-1))
          q=lambda*q1*q2
        end function

        real function dq(lambda)
          implicit none
          real lambda,q1,q2
          q1=(2/(k+1))**(-1/(k-1))
          q2=(1-(k-1)/(k+1)*lambda**2)**(1/(k-1)-1)
          dq=(1-lambda**2)*q1*q2
        end function

        real function lk(qc,k)
          implicit none
          real qc,xn,oldxn
          integer k,i

          i=0
          if (k==1) then
            xn=0.5
          else
            xn=1.5
          endif

          do while (.true.)
            oldxn=xn
            xn=oldxn-(q(oldxn)-qc)/dq(oldxn)
            if ((abs(xn-oldxn)<eps).or.(i>maxiter)) exit
            i=i+1
          end do
          lk=xn
        end function

      end program fluids
