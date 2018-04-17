      subroutine fractal(st)
c     This subroutine interpolates the fluid velocity to the particle
c     position with fractal interpolation.
c     (based on  article: Scotti, Menevau, PhysicaD 1999)
c     It makes 9 interpolations based on 27 grid points in the x 
c     direction of which 9 point are used in 3 interpolations in
c     the y directions of whitch 3 points are used in final
c     interpolation in the z direction. This gives the velocity in
c     the exact position of particle.
      
      include 'channel.cmn' 
      
c     xpar, ypar, zpar are particle position:      
c     dd1 and dd2 are scaling factors used in the interpolation and
c     are fixed the same as in the paper of Scotti and Meneveau
c     n_tab is the number of iterations of mapping procedure used in 
c     the interpolaton (different in every direction). n is fixed to
c     assure that the smallest scales of computation corresponds to
c     the Kolmogorov scale, so: 2**n ~ delta/eta; where eta is the 
c     Kolmogorov scale, delta is the resolution of the grid
c     tabU is the tablet of velocity values in the three closest
c     grid points used in the mapping procedure
c     tab_y is the tablet of 9 points obtained in the first stage of
c     interpolation which are used in the second stage
c     tab_z is similar to tab_y
c     u is the interpolated to the particle position fluid velocity 
      real*8 xpar,ypar,zpar
      real*8 dd1,dd2
      integer nn, n_tab(3), st
      parameter(dd1=2d0**(-1.0d0/3.0d0),dd2=-2d0**(-1.0d0/3.0d0))
cc      parameter(dd1=3d0**(-1.0d0/3.0d0),dd2=-3d0**(-1.0d0/3.0d0)) 
cc!fractal dimension equals 1.53
c      parameter(n_tab=(/1,2,3/)) ! czy nie powinno byc /2 1 2/?
      parameter(n_tab=(/1,1,1/)) ! czy nie powinno byc /2 1 2/?
      real*8 dyy, dxx,q, dy, dz,xi,w
      real*8 tabU(0:2), tab_y(0:2,0:2), tab_z(0:2)
      real*8 uu(3)
      integer aa, i,j,k, ii, grid_x,grid_y,grid_z,par    
      
c-----------------------------------------------------------------------
      do par=1, npar !loop over  particles

      dy = 1.0d0/(y(3)-y(1)) 
      dz = 1.0d0/(z(3)-z(1)) 

c every 2nd node
      if(st.eq.1) then
              grid_x=int(dfloat(n)/pi*acos(pos(par,1)))
              xpar = (pos(par,1)-x(grid_x))/(x(grid_x+1)-x(grid_x))
              
              q=dmod(pos(par,2)+100.0d0*width,width)
              grid_y = int(dfloat(m)*q/width)/2
              ypar=(q-y(2*grid_y))*dy
              
              q=dmod(pos(par,3)+length,length)
              grid_z = int(dfloat(o)*q/length)/2
              zpar=(q-z(2*grid_z))*dz
      end if
      
      if(st.eq.2) then
              grid_x=int(dfloat(n)/pi*acos(posnew(par,1)))
              xpar = (posnew(par,1)-x(grid_x))/(x(grid_x+1)-x(grid_x))
              
              q=dmod(posnew(par,2)+100.0d0*width,width)
              grid_y = int(dfloat(m)*q/width)/2
              ypar=(q-y(2*grid_y))*dy
              
              q=dmod(posnew(par,3)+length,length)
              grid_z = int(dfloat(o)*q/length)/2
              zpar=(q-z(2*grid_z))*dz
      end if

      do ii=1, 3 ! loop for each velocity component
      
c-----------------------------------------------------------------------
c First stage (9 iterations in the x direction):
c due to nonuniform grid points in wall-normal direction
c we will use linear interpolation there
c TODO - think how to incorporate nonuniform grids in fractal model
      tab_y = 0.0d0
      do i=0, 2
      do j=0, 2
      
              tabU(0) = upp(grid_x,2*grid_y+i,2*grid_z+j,ii)
              tabU(1) = upp(grid_x+1,2*grid_y+i,2*grid_z+j,ii)

              xi = xpar 
              tab_y(i,j) =(1-xi)*tabU(0) + xi*tabU(1)
      enddo
      enddo
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c Second stage (3 iterations in the y direction):
      tab_z = 0.0d0
      do i=0, 2
            tabU(0) = tab_y(i,0)
            tabU(1)= tab_y(i,1)
            tabU(2)= tab_y(i,2)

            xi = ypar 
            tab_z(i) = w(n_tab(2),dd1,dd2,xi,tabU)
c            tab_z(i) =(1-xi)*tabU(0) + xi*tabU(2)
      enddo
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c Third stage (one iteration in the z direction)
      tabU(0) = tab_z(0)
      tabU(1)= tab_z(1)
      tabU(2)= tab_z(2)

      xi = zpar 

      if(st.eq.1) then
          upar(par,ii) = w(n_tab(3),dd1,dd2,xi,tabU)
      end if 
      
      if(st.eq.2) then
          uparnew(par,ii) = w(n_tab(3),dd1,dd2,xi,tabU)
      end if 


      end do !end ii
      end do !end par
                  
      end subroutine fractal
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c/\/\/\/\/\\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
c Soubroutine w is the mapping procedure

      RECURSIVE FUNCTION w(nn,dd1,dd2,xi,tabU) RESULT(usgs)
          implicit none
          integer nn
          real*8 dd1,dd2,xi
          real*8 q1,q2
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


