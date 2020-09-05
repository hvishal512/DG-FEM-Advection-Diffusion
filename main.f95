!------------------1D ADVECTION-DIFFUSION SOLVER----------------------

!Here we attempt to solve using Discontinuous Galerkin method
!the 1D scalar conservation eqn of the form:

!\frac{\partial u)(\partial t} + c\frac{\partial u)(\partial x} -
!alpha\frac{\partial^2 u}{\partial x^2} = 0 with periodic BCs in [0,1]:
!u(0,t) = u(1,t) & u_x(0,t) = u_x(1,t)

PROGRAM SOLVER

    IMPLICIT NONE
    
    REAL, PARAMETER :: pi = 4*atan(1.0)
    
    INTEGER, PARAMETER :: N = 3         !Order of spatial discretization
    INTEGER, PARAMETER :: K = 16        !Number of elements
    REAL*8, PARAMETER :: c = 1          !Velocity
    REAL*8, PARAMETER :: alpha = 0.0001 !Diffusion coefficient
    REAL*8, PARAMETER :: dT = 0.001, saveT = 0.01, endT = 15
    REAL*8, DIMENSION(K+1) :: xnode
    INTEGER :: nsaveT  = NINT(saveT/dT), nT  = NINT(endT/dT)
    REAL*8, EXTERNAL :: dlagrange
    REAL*8 :: z
    REAL*8, DIMENSION (0:N) :: flux, dflux, ux1, ux2
    INTEGER :: i,m,t,j,l                !Iterating variables
    REAL, PARAMETER :: h = 1.0/K        !Grid-spacing
    REAL*8, DIMENSION(K,2) :: elembc    !Element coordinates
    REAL*8, DIMENSION(0:N) :: x,w       !LGL nodes and weights
    REAL*8, DIMENSION(0:N, 0:N) :: V, invM, Dr, S, Mass
    REAL*8, DIMENSION(K,0:N) :: mapping, basisweights
    REAL*8, DIMENSION(0:N,K) :: nodalweights, k1,k2,k3,k4
    
    
    REAL*8, DIMENSION(0:N,K,NINT(endT/saveT)+1):: saved
    REAL*8 :: start, finish             !To compute execution time
    CALL CPU_TIME(start)
    
    !1D Grid-spacing
    xnode(1:K+1) = [((i-1)*h,i=1,K+1)]
    !@ all the element boundaries, note there are repeats since nodes 
    !are coincident at intersections except at the domain boundaries
    elembc(:,1) = [((i-1)*h,i=1,K)]
    elembc(:,2) = [((i)*h,i=1,K)]

    !Precompute Legendre values (Vandermonde matrix) for fixed order
    !LGL quadrature points
    
    CALL lglnodes(x,w,N) 
    
    DO m = 0,N
        V(:,m) = legendre(m,x,N)*SQRT((2.0*m+1)/2.0)
    END DO
    
    !-According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
    !\int v \phi = \int u0 \phi for each element (xk-1/2 < x < xk+1/2)
    !-We will numerically compute RHS integral using quadrature points
    DO i = 1,K
        mapping(i,:) = (elembc(i,2)-elembc(i,1))*x/2 &
        + (elembc(i,1)+elembc(i,2))/2;
        DO m = 0,N
            !Apporimation of RHS integral with gauss quadrature sum
            basisweights(i,m) = SUM(w*SIN(2*pi*mapping(i,:))*V(:,m));
        END DO
    END DO

    !Modal to nodal conversion
    nodalweights = MATMUL(V,TRANSPOSE(basisweights))
    !--Now that we have u0 we can begin explicit time stepping with the 
    !semi-discrete form of the PDE.
    !Flux term for f(u)=cu: g(u-(x),g(u+(x))=u-(x) (Upwind)
    !Flux term for f(u)=-alpha*ux: g(u_x-(x),u_x+(x))=u_x+(x) (Downwind)
    
    invM = MATMUL(V,TRANSPOSE(V))   !Inverse of mass matrix           
    CALL inverse(MATMUL(V,TRANSPOSE(V)),Mass,N+1)
    !Derivative matrix of Lagrange interpolants
    DO i = 0,N
        z = x(i)
        DO j = 0,N
            Dr(i,j) = dlagrange(j,x,N,z)
        END DO
    END DO
    
    S = c*MATMUL(Mass,Dr)-alpha*MATMUL(Dr,TRANSPOSE(Dr)) !Stiffness

    saved(:,:,1) = nodalweights
    j = 1
    DO t= 1,nT
        !Element wise operations
        DO i = 1,K
        
            !Flux defn: Periodic BCs mean the upwind element for k=1
            !is k=K and vice-versa for the downwind element
            flux = 0.0
            dflux = 0.0
            IF (.not.(i == 1 .or. i == K)) THEN
                ux1 = MATMUL(Dr,nodalweights(:,i))
                ux2 = MATMUL(Dr,nodalweights(:,i+1))
                flux(0) = -c*nodalweights(N,i-1)
                flux(N) =  c*nodalweights(N,i)
                dflux(0) = -alpha*ux1(0)
                dflux(N) = alpha*ux2(0)
            ELSE IF (i == 1) THEN
                ux1 = MATMUL(Dr,nodalweights(:,i))
                ux2 = MATMUL(Dr,nodalweights(:,i+1))
                flux(0) = -c*nodalweights(N,K)
                flux(N) =  c*nodalweights(N,i)
                dflux(0) = -alpha*ux1(0)
                dflux(N) = alpha*ux2(0)
            ELSE 
                ux1 = MATMUL(Dr,nodalweights(:,i))
                ux2 = MATMUL(Dr,nodalweights(:,1))
                flux(0) = -c*nodalweights(N,i-1)
                flux(N) =  c*nodalweights(N,i)
                dflux(0) = -alpha*ux1(0)
                dflux(N) = alpha*ux2(0)
            END IF
            
            !Net flux
            flux = flux - dflux
            k1(:,i) = (2.0/h)*MATMUL(invM, MATMUL(TRANSPOSE(S),&
            & nodalweights(:,i)) - flux)
            k2(:,i) = (2.0/h)*MATMUL(invM, MATMUL(TRANSPOSE(S),&
            & nodalweights(:,i)+ dT*k1(:,i)/2.0) - flux)
            k3(:,i) = (2.0/h)*MATMUL(invM, MATMUL(TRANSPOSE(S),&
            & nodalweights(:,i)+dT*k2(:,i)/2.0) - flux)
            k4(:,i) = (2.0/h)*MATMUL(invM, MATMUL(TRANSPOSE(S),&
            & nodalweights(:,i)+dT*k3(:,i)) - flux)
        END DO
        !RK4 Update
        nodalweights = nodalweights + (dT/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        
        !Periodically save system state for plotting
        IF (1.0*t/nsaveT == NINT(1.0*t/nsaveT))  THEN  
            j= j+1;
            saved(:,:,j)=nodalweights;
        END IF
    END DO
    
    !Writing data to files for post-processing....................
    OPEN(unit=1, file='map.txt', ACTION="write", STATUS="replace")
    DO i = 1,K
        write(1, '(1000F14.7)')(mapping(i,j),j=0,N)
    END DO

    close(1)
    
    OPEN(unit=2, file='saved.txt', ACTION="write", STATUS="replace")
    DO l = 1, SIZE(saved,3)
        DO i = 0,N
            write(2, '(1000F14.7)')(saved(i,j,l) ,j=1,K)
        END DO
    
    END DO
    close(2)
    
    PRINT *,"Domain coordinates written to map.txt"
    PRINT *,"Function values written to saved.txt"
    CALL CPU_TIME(finish)
    print '(" Execution time",f6.3," seconds")',finish-start
    print *, "Use MATLAB for post-processing"
    CONTAINS
    
    !Function to calculate legendre values
    RECURSIVE FUNCTION legendre(m,x,N) RESULT(legf)

        IMPLICIT NONE
        INTEGER :: i,m,N
        REAL*8, DIMENSION(0:N), INTENT(IN) :: x
        REAL*8, DIMENSION(0:N) :: legf
        
        legf = 0
        IF (m == 0) THEN
            legf(0:N) = [(1,i=0,N)]
            RETURN
        ELSE IF (m == 1) THEN
            legf = x
            RETURN
        ELSE    
            legf = (2.0*m-1)/m*x*legendre(m-1,x,N) - (m-1.0)/m*&
            &legendre(m-2,x,N)
            RETURN
        END IF
        
    END FUNCTION legendre
        
END PROGRAM SOLVER

SUBROUTINE lglnodes(x,w,N)

  ! Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL 
  ! Vandermonde matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x).
  ! Useful for numerical integration and spectral methods. 
  !
  ! Reference on LGL nodes and weights:  
  ! C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral 
  ! Methods in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
  !
  ! Written by Greg von Winckel - 04/17/2004
  ! Contact: gregvw@chtm.unm.edu
  !
  ! Translated and modified not to output the Vandermonde matrix 
  ! by Daniel Appelo.  
  !
  IMPLICIT NONE
  INTEGER :: n,n1
  REAL(KIND=8) :: w(0:n),x(0:n),xold(0:n)
  REAL(KIND=8), PARAMETER :: pi = acos(-1.d0)
  INTEGER :: i,k
  REAL(KIND=8) :: P(1:n+1,1:n+1),eps
  ! Truncation + 1
  N1=N+1
  eps = 2.2204d-16
  
  ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  DO i = 0,n
     x(i) = -cos(pi*dble(i)/dble(N))
  END DO
  
  ! The Legendre Vandermonde Matrix
  !  P=zeros(N1,N1);
  
  ! Compute P_(N) using the recursion relation
  ! Compute its first and second derivatives and 
  ! update x using the Newton-Raphson method.
  
  xold = 2.d0
  
  DO i = 1,100 ! Ridic!   
     xold = x
     
     P(:,1) = 1.d0
     P(:,2) = x
     
     DO  k=2,n
        P(:,k+1)=( dble(2*k-1)*x*P(:,k)-dble(k-1)*P(:,k-1) )/dble(k);
     END DO
     x = xold-( x*P(:,N1)-P(:,N) )/( dble(N1)*P(:,N1) )
     if (maxval(abs(x-xold)).lt. eps ) EXIT
  END DO
  
  w=2.d0/(dble(N*N1)*P(:,N1)**2)
  
END SUBROUTINE lglnodes

FUNCTION dlagrange(j,x,N,z) result(dl)
        
    REAL*8 :: dl
    INTEGER, INTENT(IN) :: N,j
    REAL*8, INTENT(IN) :: z
    INTEGER :: l,m
    REAL*8 :: k
    REAL*8, DIMENSION(0:N), INTENT(IN) :: x
    dl = 0
    DO l = 0,N
        IF (.not.(l==j)) THEN
            k = 1.0/(x(j)-x(l))
            DO m = 0,N
                IF (.not.(m==j) .and. .not.(m==l)) THEN
                    k = 1.0*k*(z-x(m))/(x(j)-x(m))
                END IF
            END DO
            dl = dl + k
        END IF
    END DO
    
    RETURN
    
END FUNCTION dlagrange

subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed 
    ! during the calculation
    !===========================================================
    implicit none 
    integer n
    real*8 a(n,n), c(n,n)
    real*8 L(n,n), U(n,n), b(n), d(n), x(n)
    real*8 coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
end subroutine inverse
