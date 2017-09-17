      program heat
        implicit none
        integer,parameter::n=10
        integer,parameter::m=10
        real,parameter::h=1.0
        real,parameter::dt=0.1
        integer kk, i, j
        real,dimension(n,m)::t, tt
        real q, a
        
        call InitCond
        
        kk = 0
        do while (kk<25)
          kk = kk + 1
          
          do i=2,n-1
            do j=2,m-1
              t(n,j)=t(n-1,j)
              call Calc
            end do
          end do
          
          do i=2,n-1
            do j=2,m-1
              t(i,j)=tt(i,j)
            end do
          end do
        end do
        
        write(*, '(f12.4)') (t(i, 1), i=1,n)
        
      contains
      
        subroutine Sources
          implicit none
          if ((i>=2).and.(i<=5).and.(j>=2).and.(j<=3)) then
            q = 20
          else
            q = 0
          end if
          
          if ((i>=4).and.(i<=6).and.(j>=5).and.(j<=6)) then
            q = -10
          end if
        end subroutine
        
        subroutine InitCond
          implicit none
          do i=1,n
            do j=1,m
              t(i,j) = 1
            end do
          end do
        end subroutine
        
        subroutine Calc
          implicit none
          real t1, t2
          call Sources
          
          if ((abs(i-4)<=1).and.(abs(j-4)<=1)) then
            a = 0
          else
            a = 1
          end if
          
          t1 = a * (t(i,j+1) - 2 * t(i,j) + t(i,j-1)) * dt / (h**2)
          t2 = a * (t(i+1,j) - 2 * t(i,j) + t(i-1,j)) * dt / (h**2)
          tt(i,j) = t(i,j) + t1 + t2 + q
        end subroutine
        
      end program
