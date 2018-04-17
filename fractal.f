
      program matr_test

      real*8 dd1,dd2,w,xi
      integer nn
      parameter(dd1=2d0**(-1.0d0/3.0d0),dd2=-2d0**(-1.0d0/3.0d0))
      real*8 tabU(0:2),u_interpolated
      integer aa, i,j,k, ii, grid_x,grid_y,grid_z,par

      real*8, DIMENSION(5) ::upar_odd, x_odd
      real*8, DIMENSION(6) ::upar_even,x_even

      
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


c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c Soubroutine w is the mapping procedure


      end program


      FUNCTION w(nn,dd1,dd2,xi,tabU) RESULT(usgs)
          integer nn
          real*8 dd1,dd2,xi,u0,w_rec
          real*8 usgs, tabU(0:2)

          usgs = w_rec(nn,dd1,dd2,xi,tabU)

          RETURN

      END FUNCTION w


      RECURSIVE FUNCTION w_rec(nn,dd1,dd2,xi,tabU) RESULT(usgs)
          implicit none
          integer nn
          real*8 dd1,dd2,xi
          real*8 q1,q2
          real*8 usgs, tabU(0:2)

          if (nn.eq.0) then
             usgs =  (tabU(2)-tabU(0))*xi + tabU(0)

          else
             if (xi.le.0.5) then
                usgs = dd1 * w_rec(nn-1,dd1,dd2,2*xi,tabU)
     &                 + q1(2*xi,dd1,tabU)
             else
                usgs = dd2 * w_rec(nn-1,dd1,dd2,2*xi-1,tabU)
     &                 + q2(2*xi-1,dd2,tabU)

             endif
          endif
      END FUNCTION w_rec



c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      pure function q1(xi,dd1,tabU)
          implicit none
          real*8, intent(in) :: xi,dd1
          real*8 , intent(in) ::tabU(0:2)
          real*8 q1

          q1=xi*(tabU(1)-tabU(0)-dd1*(tabU(2)-tabU(0)))
     &          + tabU(0)*(1-dd1)

      end function
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      pure function q2(xi,dd2,tabU)
          implicit none
          real*8, intent(in) :: xi,dd2
          real*8 , intent(in) ::tabU(0:2)
          real*8 q2

          q2=xi*(tabU(2)-tabU(1)-dd2*(tabU(2)-tabU(0)))
     &          + tabU(1)-dd2*tabU(0)

      end function
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



