
      module affine_transform
      implicit none

      contains

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

      RECURSIVE FUNCTION w(nn,dd1,dd2,xi,tabU) RESULT(usgs)
          implicit none
          integer nn
          real*8 dd1,dd2,xi
          real*8 usgs, tabU(0:2)

          if (nn.eq.0) then
             usgs =  (tabU(2)-tabU(0))*xi + tabU(0)

          else
             if (xi.le.0.5) then
                usgs = dd1 * w(nn-1,dd1,dd2,2*xi,tabU) 
     &                 + q1(2*xi,dd1,tabU)
             else
                usgs = dd2 * w(nn-1,dd1,dd2,2*xi-1,tabU) 
     &                 + q2(2*xi-1,dd2,tabU)

             endif
          endif
      END FUNCTION w



      end module affine_transform

