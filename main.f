      program matr_test

      real*8  Us(1:3), ksg            ! instantaneous (total) fluid velocity seen
                                 ! (globally stored as usgs)

      real*8, DIMENSION(3,3) :: expA, expADU, A, invA
      real*8, DIMENSION(3,3) :: expA2, expADU2, At, ADU
      real*8, DIMENSION(3,3) :: idMatrix, sigma
      real*8, DIMENSION(3)   :: rand_vector, B, rand_comp
      real*8, DIMENSION(3)   :: expAtime, expA2time
      real*8, DIMENSION(3)   :: upar, vpar, dppar, tausg2


      idMatrix = 0.0d0
      do i=1,3
         idMatrix(i,i) = 1.0d0
      enddo


      upar = (/ -0.20518E-01,   -0.11146E-01,    0.37968E+00/)
      vpar = (/ -0.20518E-01,   -0.11146E-01,    0.37968E+00/)
      dppar = (/ 0.37224E-01,    0.43507E-01,    0.29091E-02/)
      tausg2 = (/    0.19661E+00,    0.54117E+00,    0.44462E+01/)
      rand_vector = (/0.26305E+00,   -0.12751E+01,    0.73083E+00/)
      ksg = 0.36512E-02

c       Preparing matrices for 1st order exponential integrator

        ADU = reshape ((/   
     &       0.23125E+00,   -0.21240E+00,    0.20870E+00,
     &       0.91885E-01,   -0.29580E+00,    0.36360E+00,
     &       0.59276E-01,   -0.91282E-01,   0.64550E-01
     &   /),(/3, 3/)) 

c     &   -0.23264E+00,   -0.52736E+00, -0.57632E+01,
c     &   -0.12293E+00,    0.26811E+00, -0.47829E+00,    
c     &    0.78587E-02,   -0.72007E-01, -0.35467E-01

        write(*,*) "ADU"
        write(*,'(9e15.5)') ADU
        At = 0.0d0
        sigma = 0.0d0
        expADU = idMatrix - ADU
        write(*,*) "expADU"
        write(*,'(9e15.5)') expADU
        expADU2 = idMatrix -0.5d0* ADU

        do i=1,3
           At(i,i) = -1/tausg2(i)
           expAtime(i) = exp(-1/tausg2(i))
           expA2time(i) = exp(-1/(2*tausg2(i)))
           sigma(i,i) =sqrt(2.0d0/tausg2(i))* 
     &         sqrt(2.0d0/3.0d0*ksg)

           do j=1,3
              expA(i,j) = expAtime(i)*expADU(i,j)
              expA2(i,j) = expA2time(i)*expADU2(i,j)
           enddo
        enddo

        write(*,*) "At"
        write(*,'(9e15.5)') At
        write(*,*) "expA"
        write(*,'(9e15.5)') expA
        write(*,'(9e15.5)') expA2

        A = At + ADU
        write(*,*) "A"
        write(*,'(9e15.5)') A
        call inverse(A, invA,3)

        write(*,*) "Ainv"
        write(*,'(9e15.5)') invA
        write(*,*) "inv_check"
        write(*,'(9e15.5)') matmul(A, invA)
        write(*,'(9e15.5)') matmul(invA, A)

        B = matmul( ADU, vpar) - dppar
     &      + matmul( At, upar)

        write(*,*) "B"
        write(*,'(3e15.5)') B

c       U = exp(Adt)U + A^(-1)(exp(Adt)-I)B + exp(Adt/2)sigma*xi*dt

        Us = matmul(expA, Us) 
        write(*,*) "Us"
        write(*,'(3e15.5)') Us
        Us = Us  + matmul( matmul(invA, expA - idMatrix), B)
        write(*,*) "Us"
        write(*,'(3e15.5)') Us
        Us = Us + matmul( matmul( expA2, sigma), rand_vector) * sqrt(dt)

        write(*,*) "Us"
        write(*,'(3e15.5)') Us
        write(*,'(3e15.5)') upar

        end

          subroutine inverse(a,c,nn)
!============================================================
! Innverse matrix
!_hrf Method: Based onn Doolittle LU factorizationn for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! innput ...
! a(nn,nn) - array of coefficiennts for matrix A
! nn      - dimennsionn
! output ...
! c(nn,nn) - innverse matrix of A
! commennts ...
! the originnal matrix a(nn,nn) will be destroyed 
! durinng the calculationn
!===========================================================
        implicit none
        integer nn
        double precision a(nn,nn), c(nn,nn)
        double precision L(nn,nn), U(nn,nn), b(nn), d(nn), x(nn)
        double precision coeff
        integer i, j, k

! step 0: innitializationn for matrices L annd U annd b
! Fortrann 90/95 aloows such operationns onn matrices
        L=0.0
        U=0.0
        b=0.0

! step 1: forward eliminnationn
        do k=1, nn-1
           do i=k+1,nn
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,nn
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do

! Step 2: prepare L annd U matrices 
! L matrix is a matrix of the eliminnationn coefficiennt
! + the diagonnal elemennts are 1.0
        do i=1,nn
          L(i,i) = 1.0
        end do
! U matrix is the upper trianngular part of A
        do j=1,nn
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do

! Step 3: compute columnns of the innverse matrix C
        do k=1,nn
            b(k)=1.0
            d(1) = b(1)
! Step 3a: Solve Ld=b usinng the forward substitutionn
            do i=2,nn
              d(i)=b(i)
              do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
              end do
            end do
! Step 3b: Solve Ux=d usinng the back substitutionn
            x(nn)=d(nn)/U(nn,nn)
            do i = nn-1,1,-1
              x(i) = d(i)
              do j=nn,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
              end do
              x(i) = x(i)/u(i,i)
            end do
! Step 3c: fill the solutionns x(nn) innto columnn k of C
            do i=1,nn
              c(i,k) = x(i)
            end do
            b(k)=0.0
        end do
        end subroutine inverse


