      module frac_dim_helpers

      contains

      PURE FUNCTION linear_coeffs(x0,y0,x1,y1)
          implicit none
          real*8, intent(in) :: x0,y0,x1,y1
          real*8, dimension(2) :: linear_coeffs

          linear_coeffs = (/ (y1-y0)/(x1-x0), (y0*x1-y1*x0)/(x1-x0) /)

      END FUNCTION linear_coeffs

      end module frac_dim_helpers



      program linear_test

      use frac_dim_helpers
      real*8 line_02(2), line_01(2), line_12(2)

      line_02 = linear_coeffs(1.0d0,2.0d0,1.5d0,3.0d0)
      
      

      write(*,*) line_02

      end
