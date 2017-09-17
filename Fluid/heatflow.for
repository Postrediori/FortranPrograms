      program heatflow
        implicit none
        real,parameter::l1=3.0
        real,parameter::l2=1.0
        real,parameter::h=0.01
        integer,parameter::n1=300 !15
        integer,parameter::n2=100 !5
        real,parameter::k1=2.0
        real,parameter::k2=1.0
        real,parameter::q=1
        real,parameter::eps=1e-5
        real,dimension(0:n1,0:n2)::u1
        real,dimension(0:n1,0:n2)::u2
        real un
        real m1,m2
        logical flag
        integer i,j
        integer k
        
        u1(0:n1,0:n2)=0.0
        u2=u1
        
        k=0
        do while (.true.)
          do i=1,n1-1
            do j=1,n2-1
              un=k1*(u1(i+1,j)+u1(i-1,j))+k2*(u1(i,j+1)+u1(i,j-1))
              u2(i,j)=(un-h*h*q)/(2.0*(k1+k2))
            end do
          end do
          
          !do j=0,n2
          !  print '(6(f10.3,'';''))',u2(0:n1,j)
          !end do
          
          flag=.false.
          do i=0,n1
            do j=0,n2
              if (flag) then
                if (abs(u1(i,j))>m1) m1=abs(u1(i,j))
                if (abs(u2(i,j))>m2) m2=abs(u2(i,j))
              else
                m1=abs(u1(i,j))
                m2=abs(u2(i,j))
                flag=.true.
              end if
            end do
          end do
        
          u1=u2
          
          if (abs(m2-m1)<eps) exit
          k=k+1
        end do
        
      ! flag=.false.
      ! do i=0,n1
      !   do j=0,n2
      !     if (flag) then
      !       if (u1(i,j)>m1) m1=u1(i,j)
      !     else
      !       m1=u1(i,j)
      !       flag=.true.
      !     end if
      !   end do
      ! end do    
      ! print*,m1

      ! do i=0,n1
      !   print '(6(f10.3,'';''))',u1(i,0:n2)
      ! end do    
        print '(''Iterations: '',i8)',k
        print '(''Minimal u: '',f10.3)',m1
        
      end program
