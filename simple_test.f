      
      SUBROUTINE simple_test

              use affine_transform

              real*8, DIMENSION(5) ::upar_odd, x_odd
              real*8, DIMENSION(6) ::upar_even,x_even
              real*8 tabU(0:2)
              real*8 dd1,dd2,xi, u_interpolated
              integer nn, i, j

              
              x_odd = (/ 0.0 ,  0.25 ,  0.5 ,  0.75 ,  1.0 /)
              x_even = (/ 0.0 ,  0.25 ,  0.5 ,  0.75 ,  1.0, 1.25 /)
              upar_odd = (/ 0.020883 ,  0.021884 ,  0.036644 ,  
     &              0.025015 ,  0.031795 /)
              upar_even = (/ 0.020883 ,  0.021884 ,  0.036644 ,  
     &          0.025015 ,  0.031795, -0.013087 /)

                      
              nn = 1 !number of mapping iteration for ii component

              do i=1, 5/2
c-----------------------------------------------------------------------
c Third stage (one iteration in the z direction)
                      tabU(0) =upar_odd(2*i-1)
                      tabU(1)=upar_odd(2*i)
                      tabU(2)=upar_odd(2*i+1)

                      write(*,*) tabU
                      do j=0,10
                              xi = 0.1*j 

                              u_interpolated  =  w(nn,dd1,dd2,xi,tabU)
                              write(*,*)xi, u_interpolated
                      enddo
              end do !end par


      END SUBROUTINE simple_test

