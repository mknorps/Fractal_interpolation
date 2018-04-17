
      program matr_test

      use affine_transform

      real*8 dd1,dd2,xi
      integer nn
      parameter(dd1=2d0**(-1.0d0/3.0d0),dd2=-2d0**(-1.0d0/3.0d0))
      real*8 tabU(0:2),tab_y(0:2, 0:2), tab_z(0:2)
      real*8 upp(0:1,0:2,0:2)
      real*8 u_interpolated
      integer aa, i,j,k, ii, grid_x,grid_y,grid_z,par
      integer n_tab(3)
      parameter(n_tab=(/1,1,1/))

      real*8 xpar,ypar,zpar


      call RANDOM_NUMBER(upp)

      xpar = 0.5
      ypar = 0.5
      zpar = 0.5

      write(*,*) upp
c-----------------------------------------------------------------------
c First stage (9 iterations in the x direction):
c due to nonuniform grid points in wall-normal direction
c we will use linear interpolation there
c TODO - think how to incorporate nonuniform grids in fractal model
      tab_y = 0.0d0
      do i=0, 2
      do j=0, 2
      
              tabU(0) = upp(0,i,j)
              tabU(1) = upp(1,i,j)

              xi = xpar 
              tab_y(i,j) =(1-xi)*tabU(0) + xi*tabU(1)
              write(*,*) i,j,"xpar",tabU,tab_y(i,j)
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
            write(*,*) i,"ypar",tabU,tab_z(i)
      enddo
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c Third stage (one iteration in the z direction)
      tabU(0) = tab_z(0)
      tabU(1)= tab_z(1)
      tabU(2)= tab_z(2)

      xi = zpar 

      par_velocity = w(n_tab(3),dd1,dd2,xi,tabU)
      write(*,*) "zpar",tabU,par_velocity
      

      end program

