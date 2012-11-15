      subroutine surf1 (m,n,x,y,z,iz,zx1,zxm,zy1,zyn,zxy11,
     *                  zxym1,zxy1n,zxymn,islpsw,zp,temp,
     *                  sigma,ierr)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer m,n,iz,islpsw,ierr
      real x(m),y(n),z(iz,n),zx1(n),zxm(n),zy1(m),zyn(m),
     *     zxy11,zxym1,zxy1n,zxymn,zp(m,n,3),temp(1),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the parameters necessary to
c compute an interpolatory surface passing through a rect-
c angular grid of functional values. the surface determined
c can be represented as the tensor product of splines under
c tension. the x- and y-partial derivatives around the
c boundary and the x-y-partial derivatives at the four
c corners may be specified or omitted. for actual mapping
c of points onto the surface it is necessary to call the
c function surf2.
c
c on input--
c
c   m is the number of grid lines in the x-direction, i. e.
c   lines parallel to the y-axis (m .ge. 2).
c
c   n is the number of grid lines in the y-direction, i. e.
c   lines parallel to the x-axis (n .ge. 2).
c
c   x is an array of the m x-coordinates of the grid lines
c   in the x-direction. these should be strictly increasing.
c
c   y is an array of the n y-coordinates of the grid lines
c   in the y-direction. these should be strictly increasing.
c
c   z is an array of the m * n functional values at the grid
c   points, i. e. z(i,j) contains the functional value at
c   (x(i),y(j)) for i = 1,...,m and j = 1,...,n.
c
c   iz is the row dimension of the matrix z used in the
c   calling program (iz .ge. m).
c
c   zx1 and zxm are arrays of the m x-partial derivatives
c   of the function along the x(1) and x(m) grid lines,
c   respectively. thus zx1(j) and zxm(j) contain the x-part-
c   ial derivatives at the points (x(1),y(j)) and
c   (x(m),y(j)), respectively, for j = 1,...,n. either of
c   these parameters will be ignored (and approximations
c   supplied internally) if islpsw so indicates.
c
c   zy1 and zyn are arrays of the n y-partial derivatives
c   of the function along the y(1) and y(n) grid lines,
c   respectively. thus zy1(i) and zyn(i) contain the y-part-
c   ial derivatives at the points (x(i),y(1)) and
c   (x(i),y(n)), respectively, for i = 1,...,m. either of
c   these parameters will be ignored (and estimations
c   supplied internally) if islpsw so indicates.
c
c   zxy11, zxym1, zxy1n, and zxymn are the x-y-partial
c   derivatives of the function at the four corners,
c   (x(1),y(1)), (x(m),y(1)), (x(1),y(n)), and (x(m),y(n)),
c   respectively. any of the parameters will be ignored (and
c   estimations supplied internally) if islpsw so indicates.
c
c   islpsw contains a switch indicating which boundary
c   derivative information is user-supplied and which
c   should be estimated by this subroutine. to determine
c   islpsw, let
c        i1 = 0 if zx1 is user-supplied (and = 1 otherwise),
c        i2 = 0 if zxm is user-supplied (and = 1 otherwise),
c        i3 = 0 if zy1 is user-supplied (and = 1 otherwise),
c        i4 = 0 if zyn is user-supplied (and = 1 otherwise),
c        i5 = 0 if zxy11 is user-supplied
c                                       (and = 1 otherwise),
c        i6 = 0 if zxym1 is user-supplied
c                                       (and = 1 otherwise),
c        i7 = 0 if zxy1n is user-supplied
c                                       (and = 1 otherwise),
c        i8 = 0 if zxymn is user-supplied
c                                       (and = 1 otherwise),
c   then islpsw = i1 + 2*i2 + 4*i3 + 8*i4 + 16*i5 + 32*i6
c                   + 64*i7 + 128*i8
c   thus islpsw = 0 indicates all derivative information is
c   user-supplied and islpsw = 255 indicates no derivative
c   information is user-supplied. any value between these
c   limits is valid.
c
c   zp is an array of at least 3*m*n locations.
c
c   temp is an array of at least n+n+m locations which is
c   used for scratch storage.
c
c and
c
c   sigma contains the tension factor. this value indicates
c   the curviness desired. if abs(sigma) is nearly zero
c   (e. g. .001) the resulting surface is approximately the
c   tensor product of cubic splines. if abs(sigma) is large
c   (e. g. 50.) the resulting surface is approximately
c   bi-linear. if sigma equals zero tensor products of
c   cubic splines result. a standard value for sigma is
c   approximately 1. in absolute value.
c
c on output--
c
c   zp contains the values of the xx-, yy-, and xxyy-partial
c   derivatives of the surface at the given nodes.
c
c   ierr contains an error flag,
c        = 0 for normal return,
c        = 1 if n is less than 2 or m is less than 2,
c        = 2 if the x-values or y-values are not strictly
c            increasing.
c
c and
c
c   m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1,
c   zxy1n, zxymn, islpsw, and sigma are unaltered.
c
c this subroutine references package modules ceez, terms,
c and snhcsh.
c
c-----------------------------------------------------------
c
      mm1 = m-1
      mp1 = m+1
      nm1 = n-1
      np1 = n+1
      npm = n+m
      ierr = 0
      if (n .le. 1 .or. m .le. 1) go to 46
      if (y(n) .le. y(1)) go to 47
c
c denormalize tension factor in y-direction
c
      sigmay = abs(sigma)*float(n-1)/(y(n)-y(1))
c
c obtain y-partial derivatives along y = y(1)
c
      if ((islpsw/8)*2 .ne. (islpsw/4)) go to 2
      do 1 i = 1,m
    1   zp(i,1,1) = zy1(i)
      go to 5
    2 dely1 = y(2)-y(1)
      dely2 = dely1+dely1
      if (n .gt. 2) dely2 = y(3)-y(1)
      if (dely1 .le. 0. .or. dely2 .le. dely1) go to 47
      call ceez (dely1,dely2,sigmay,c1,c2,c3,n)
      do 3 i = 1,m
    3   zp(i,1,1) = c1*z(i,1)+c2*z(i,2)
      if (n .eq. 2) go to 5
      do i = 1,m
         zp(i,1,1) = zp(i,1,1)+c3*z(i,3)
         zy1(i)=zp(i,1,1)
      enddo
c
c obtain y-partial derivatives along y = y(n)
c
    5 if ((islpsw/16)*2 .ne. (islpsw/8)) go to 7
      do 6 i = 1,m
        npi = n+i
    6   temp(npi) = zyn(i)
      go to 10
    7 delyn = y(n)-y(nm1)
      delynm = delyn+delyn
      if (n .gt. 2) delynm = y(n)-y(n-2)
      if (delyn .le. 0. .or. delynm .le. delyn) go to 47
      call ceez (-delyn,-delynm,sigmay,c1,c2,c3,n)
      do 8 i = 1,m
        npi = n+i
    8   temp(npi) = c1*z(i,n)+c2*z(i,nm1)
      if (n .eq. 2) go to 10
      do 9 i = 1,m
        npi = n+i
    9   temp(npi) = temp(npi)+c3*z(i,n-2)
   10 if (x(m) .le. x(1)) go to 47
c
c denormalize tension factor in x-direction
c
      sigmax = abs(sigma)*float(m-1)/(x(m)-x(1))
c
c obtain x-partial derivatives along x = x(1)
c
      if ((islpsw/2)*2 .ne. islpsw) go to 12
      do 11 j = 1,n
   11   zp(1,j,2) = zx1(j)
      if ((islpsw/32)*2 .eq. (islpsw/16) .and.
     *    (islpsw/128)*2  .eq. (islpsw/64)) go to 15
   12 delx1 = x(2)-x(1)
      delx2 = delx1+delx1
      if (m .gt. 2) delx2 = x(3)-x(1)
      if (delx1 .le. 0. .or. delx2 .le. delx1) go to 47
      call ceez (delx1,delx2,sigmax,c1,c2,c3,m)
      if ((islpsw/2)*2 .eq. islpsw) go to 15
      do 13 j = 1,n
   13   zp(1,j,2) = c1*z(1,j)+c2*z(2,j)
      if (m .eq. 2) go to 15
      do j = 1,n
         zp(1,j,2) = zp(1,j,2)+c3*z(3,j)
         zx1(j)=zp(1,j,2)
      enddo
c
c obtain x-y-partial derivative at (x(1),y(1))
c
   15 if ((islpsw/32)*2 .ne. (islpsw/16)) go to 16
      zp(1,1,3) = zxy11
      go to 17
   16 zp(1,1,3) = c1*zp(1,1,1)+c2*zp(2,1,1)
      if (m .gt. 2) zp(1,1,3) = zp(1,1,3)+c3*zp(3,1,1)
      zxy11=zp(1,1,3)
c
c obtain x-y-partial derivative at (x(1),y(n))
c
   17 if ((islpsw/128)*2 .ne. (islpsw/64)) go to 18
      zxy1ns = zxy1n
      go to 19
   18 zxy1ns = c1*temp(n+1)+c2*temp(n+2)
      if (m .gt. 2) zxy1ns = zxy1ns+c3*temp(n+3)
c
c obtain x-partial derivative along x = x(m)
c
   19 if ((islpsw/4)*2 .ne. (islpsw/2)) go to 21
      do 20 j = 1,n
        npmpj = npm+j
   20   temp(npmpj) = zxm(j)
      if ((islpsw/64)*2 .eq. (islpsw/32) .and.
     *    (islpsw/256)*2 .eq. (islpsw/128)) go to 24
   21 delxm = x(m)-x(mm1)
      delxmm = delxm+delxm
      if (m .gt. 2) delxmm = x(m)-x(m-2)
      if (delxm .le. 0. .or. delxmm .le. delxm) go to 47
      call ceez (-delxm,-delxmm,sigmax,c1,c2,c3,m)
      if ((islpsw/4)*2 .eq. (islpsw/2)) go to 24
      do 22 j = 1,n
        npmpj = npm+j
   22   temp(npmpj) = c1*z(m,j)+c2*z(mm1,j)
      if (m .eq. 2) go to 24
      do 23 j = 1,n
        npmpj = npm+j
   23   temp(npmpj) = temp(npmpj)+c3*z(m-2,j)
c
c obtain x-y-partial derivative at (x(m),y(1))
c
   24 if ((islpsw/64)*2 .ne. (islpsw/32)) go to 25
      zp(m,1,3) = zxym1
      go to 26
   25 zp(m,1,3) = c1*zp(m,1,1)+c2*zp(mm1,1,1)
      if (m .gt. 2) zp(m,1,3) = zp(m,1,3)+c3*zp(m-2,1,1)
c
c obtain x-y-partial derivative at (x(m),y(n))
c
   26 if ((islpsw/256)*2 .ne. (islpsw/128)) go to 27
      zxymns = zxymn
      go to 28
   27 zxymns = c1*temp(npm)+c2*temp(npm-1)
      if (m .gt. 2) zxymns = zxymns+c3*temp(npm-2)
c
c set up right hand sides and tridiagonal system for y-grid
c perform forward elimination
c
   28 del1 = y(2)-y(1)
      if (del1 .le. 0.) go to 47
      deli = 1./del1
      do 29 i = 1,m
   29   zp(i,2,1) = deli*(z(i,2)-z(i,1))
      zp(1,2,3) = deli*(zp(1,2,2)-zp(1,1,2))
      zp(m,2,3) = deli*(temp(npm+2)-temp(npm+1))
      call terms (diag1,sdiag1,sigmay,del1)
      diagi = 1./diag1
      do 30 i = 1,m
   30   zp(i,1,1) = diagi*(zp(i,2,1)-zp(i,1,1))
      zp(1,1,3) = diagi*(zp(1,2,3)-zp(1,1,3))
      zp(m,1,3) = diagi*(zp(m,2,3)-zp(m,1,3))
      temp(1) = diagi*sdiag1
      if (n .eq. 2) go to 34
      do 33 j = 2,nm1
        jm1 = j-1
        jp1 = j+1
        npmpj = npm+j
        del2 = y(jp1)-y(j)
        if (del2 .le. 0.) go to 47
        deli = 1./del2
        do 31 i = 1,m
   31     zp(i,jp1,1) = deli*(z(i,jp1)-z(i,j))
        zp(1,jp1,3) = deli*(zp(1,jp1,2)-zp(1,j,2))
        zp(m,jp1,3) = deli*(temp(npmpj+1)-temp(npmpj))
        call terms (diag2,sdiag2,sigmay,del2)
        diagin = 1./(diag1+diag2-sdiag1*temp(jm1))
        do 32 i = 1,m
   32     zp(i,j,1) = diagin*(zp(i,jp1,1)-zp(i,j,1)-
     *                        sdiag1*zp(i,jm1,1))
        zp(1,j,3) = diagin*(zp(1,jp1,3)-zp(1,j,3)-
     *                      sdiag1*zp(1,jm1,3))
        zp(m,j,3) = diagin*(zp(m,jp1,3)-zp(m,j,3)-
     *                      sdiag1*zp(m,jm1,3))
        temp(j) = diagin*sdiag2
        diag1 = diag2
   33   sdiag1 = sdiag2
   34 diagin = 1./(diag1-sdiag1*temp(nm1))
      do 35 i = 1,m
        npi = n+i
   35   zp(i,n,1) = diagin*(temp(npi)-zp(i,n,1)-
     *                      sdiag1*zp(i,nm1,1))
      zp(1,n,3) = diagin*(zxy1ns-zp(1,n,3)-
     *                    sdiag1*zp(1,nm1,3))
      temp(n) = diagin*(zxymns-zp(m,n,3)-
     *                  sdiag1*zp(m,nm1,3))
c
c perform back substitution
c
      do 37 j = 2,n
        jbak = np1-j
        jbakp1 = jbak+1
        t = temp(jbak)
        do 36 i = 1,m
   36     zp(i,jbak,1) = zp(i,jbak,1)-t*zp(i,jbakp1,1)
        zp(1,jbak,3) = zp(1,jbak,3)-t*zp(1,jbakp1,3)
   37   temp(jbak) = zp(m,jbak,3)-t*temp(jbakp1)
c
c set up right hand sides and tridiagonal system for x-grid
c perform forward elimination
c
      del1 = x(2)-x(1)
      if (del1 .le. 0.) go to 47
      deli = 1./del1
      do 38 j = 1,n
        zp(2,j,2) = deli*(z(2,j)-z(1,j))
   38   zp(2,j,3) = deli*(zp(2,j,1)-zp(1,j,1))
      call terms (diag1,sdiag1,sigmax,del1)
      diagi = 1./diag1
      do 39 j = 1,n
        zp(1,j,2) = diagi*(zp(2,j,2)-zp(1,j,2))
   39   zp(1,j,3) = diagi*(zp(2,j,3)-zp(1,j,3))
      temp(n+1) = diagi*sdiag1
      if (m  .eq. 2) go to 43
      do 42 i = 2,mm1
        im1 = i-1
        ip1 = i+1
        npi = n+i
        del2 = x(ip1)-x(i)
        if (del2 .le. 0.) go to 47
        deli = 1./del2
        do 40 j = 1,n
          zp(ip1,j,2) = deli*(z(ip1,j)-z(i,j))
   40     zp(ip1,j,3) = deli*(zp(ip1,j,1)-zp(i,j,1))
        call terms (diag2,sdiag2,sigmax,del2)
        diagin = 1./(diag1+diag2-sdiag1*temp(npi-1))
        do 41 j = 1,n
          zp(i,j,2) = diagin*(zp(ip1,j,2)-zp(i,j,2)-
     *                        sdiag1*zp(im1,j,2))
   41     zp(i,j,3) = diagin*(zp(ip1,j,3)-zp(i,j,3)-
     *                        sdiag1*zp(im1,j,3))
        temp(npi) = diagin*sdiag2
        diag1 = diag2
   42   sdiag1 = sdiag2
   43 diagin = 1./(diag1-sdiag1*temp(npm-1))
      do 44 j = 1,n
        npmpj = npm+j
        zp(m,j,2) = diagin*(temp(npmpj)-zp(m,j,2)-
     *                      sdiag1*zp(mm1,j,2))
   44   zp(m,j,3) = diagin*(temp(j)-zp(m,j,3)-
     *                      sdiag1*zp(mm1,j,3))
c
c perform back substitution
c
      do 45 i = 2,m
        ibak = mp1-i
        ibakp1 = ibak+1
        npibak = n+ibak
        t = temp(npibak)
        do 45 j = 1,n
          zp(ibak,j,2) = zp(ibak,j,2)-t*zp(ibakp1,j,2)
   45     zp(ibak,j,3) = zp(ibak,j,3)-t*zp(ibakp1,j,3)
      return
c
c too few points
c
   46 ierr = 1
      return
c
c points not strictly increasing
c
   47 ierr = 2
      return
      end



      function surf2 (xx,yy,m,n,x,y,z,iz,zp,sigma)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer m,n,iz
      real xx,yy,x(m),y(n),z(iz,n),zp(m,n,3),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function interpolates a surface at a given coordinate
c pair using a bi-spline under tension. the subroutine surf1
c should be called earlier to determine certain necessary
c parameters.
c
c on input--
c
c   xx and yy contain the x- and y-coordinates of the point
c   to be mapped onto the interpolating surface.
c
c   m and n contain the number of grid lines in the x- and
c   y-directions, respectively, of the rectangular grid
c   which specified the surface.
c
c   x and y are arrays containing the x- and y-grid values,
c   respectively, each in increasing order.
c
c   z is a matrix containing the m * n functional values
c   corresponding to the grid values (i. e. z(i,j) is the
c   surface value at the point (x(i),y(j)) for i = 1,...,m
c   and j = 1,...,n).
c
c   iz contains the row dimension of the array z as declared
c   in the calling program.
c
c   zp is an array of 3*m*n locations stored with the
c   various surface derivative information determined by
c   surf1.
c
c and
c
c   sigma contains the tension factor (its sign is ignored).
c
c the parameters m, n, x, y, z, iz, zp, and sigma should be
c input unaltered from the output of surf1.
c
c on output--
c
c   surf2 contains the interpolated surface value.
c
c none of the input parameters are altered.
c
c this function references package modules intrvl and
c snhcsh.
c
c-----------------------------------------------------------
c
c inline one dimensional cubic spline interpolation
c
      hermz (f1,f2,fp1,fp2) = (f2*del1+f1*del2)/dels-del1*
     *                        del2*(fp2*(del1+dels)+
     *                              fp1*(del2+dels))/
     *                        (6.*dels)
c
c inline one dimensional spline under tension interpolation
c
      hermnz (f1,f2,fp1,fp2,sigmap) = (f2*del1+f1*del2)/dels
     *          +(fp2*del1*(sinhm1-sinhms)
     *           +fp1*del2*(sinhm2-sinhms)
     *          )/(sigmap*sigmap*dels*(1.+sinhms))
c
c denormalize tension factor in x and y direction
c
      sigmax = abs(sigma)*float(m-1)/(x(m)-x(1))
      sigmay = abs(sigma)*float(n-1)/(y(n)-y(1))
c
c determine y interval
c
      jm1 = intrvl (yy,y,n)
      j = jm1+1
c
c determine x interval
c
      im1 = intrvl (xx,x,m)
      i = im1+1
      del1 = yy-y(jm1)
      del2 = y(j)-yy
      dels = y(j)-y(jm1)
      if (sigmay .ne. 0.) go to 1
c
c perform four interpolations in y-direction
c
      zim1 = hermz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),
     *                                  zp(i-1,j,1))
      zi = hermz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1))
      zxxim1 = hermz(zp(i-1,j-1,2),zp(i-1,j,2),
     *                zp(i-1,j-1,3),zp(i-1,j,3))
      zxxi = hermz(zp(i,j-1,2),zp(i,j,2),
     *              zp(i,j-1,3),zp(i,j,3))
      go to 2
    1 call snhcsh (sinhm1,dummy,sigmay*del1,-1)
      call snhcsh (sinhm2,dummy,sigmay*del2,-1)
      call snhcsh (sinhms,dummy,sigmay*dels,-1)
      zim1 = hermnz(z(i-1,j-1),z(i-1,j),zp(i-1,j-1,1),
     *               zp(i-1,j,1),sigmay)
      zi = hermnz(z(i,j-1),z(i,j),zp(i,j-1,1),zp(i,j,1),
     *             sigmay)
      zxxim1 = hermnz(zp(i-1,j-1,2),zp(i-1,j,2),
     *                 zp(i-1,j-1,3),zp(i-1,j,3),sigmay)
      zxxi = hermnz(zp(i,j-1,2),zp(i,j,2),
     *               zp(i,j-1,3),zp(i,j,3),sigmay)
c
c perform final interpolation in x-direction
c
    2 del1 = xx-x(im1)
      del2 = x(i)-xx
      dels = x(i)-x(im1)
      if (sigmax .ne. 0.) go to 3
      surf2 = hermz(zim1,zi,zxxim1,zxxi)
      return
    3 call snhcsh (sinhm1,dummy,sigmax*del1,-1)
      call snhcsh (sinhm2,dummy,sigmax*del2,-1)
      call snhcsh (sinhms,dummy,sigmax*dels,-1)
      surf2 = hermnz(zim1,zi,zxxim1,zxxi,sigmax)
      return
      end



      subroutine ceez (del1,del2,sigma,c1,c2,c3,n)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      real del1,del2,sigma,c1,c2,c3
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the coefficients c1, c2, and c3
c used to determine endpoint slopes. specifically, if
c function values y1, y2, and y3 are given at points x1, x2,
c and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3
c is the value of the derivative at x1 of a spline under
c tension (with tension factor sigma) passing through the
c three points and having third derivative equal to zero at
c x1. optionally, only two values, c1 and c2 are determined.
c
c on input--
c
c   del1 is x2-x1 (.gt. 0.).
c
c   del2 is x3-x1 (.gt. 0.). if n .eq. 2, this parameter is
c   ignored.
c
c   sigma is the tension factor.
c
c and
c
c   n is a switch indicating the number of coefficients to
c   be returned. if n .eq. 2 only two coefficients are
c   returned. otherwise all three are returned.
c
c on output--
c
c   c1, c2, and c3 contain the coefficients.
c
c none of the input parameters are altered.
c
c this subroutine references package module snhcsh.
c
c-----------------------------------------------------------
c
      if (n .eq. 2) go to 2
      if (sigma .ne. 0.) go to 1
      del = del2-del1
c
c tension .eq. 0.
c
      c1 = -(del1+del2)/(del1*del2)
      c2 = del2/(del1*del)
      c3 = -del1/(del2*del)
      return
c
c tension .ne. 0.
c
    1 call snhcsh (dummy,coshm1,sigma*del1,1)
      call snhcsh (dummy,coshm2,sigma*del2,1)
      delp = sigma*(del2+del1)/2.
      delm = sigma*(del2-del1)/2.
      call snhcsh (sinhmp,dummy,delp,-1)
      call snhcsh (sinhmm,dummy,delm,-1)
      denom = coshm1*(del2-del1)-2.*del1*delp*delm*
     *        (1.+sinhmp)*(1.+sinhmm)
      c1 = 2.*delp*delm*(1.+sinhmp)*(1.+sinhmm)/denom
      c2 = -coshm2/denom
      c3 = coshm1/denom
      return
c
c two coefficients
c
    2 c1 = -1./del1
      c2 = -c1
      return
      end



      subroutine curvss (n,x,y,d,isw,s,eps,ys,ysp,sigma,td,
     *                   tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,
     *                   ierr)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer n,isw,ierr
      real x(n),y(n),d(n),s,eps,ys(n),ysp(n),sigma,td(n),
     *     tsd1(n),hd(n),hsd1(n),hsd2(n),rd(n),rsd1(n),
     *     rsd2(n),v(n)
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the parameters necessary to
c compute a smoothing spline under tension. for a given
c increasing sequence of abscissae (x(i)), i = 1,..., n and
c associated ordinates (y(i)), i = 1,..., n, the function
c determined minimizes the summation from i = 1 to n-1 of
c the square of the second derivative of f plus sigma
c squared times the difference of the first derivative of f
c and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
c functions f with two continuous derivatives such that the
c summation of the square of (f(x(i))-y(i))/d(i) is less
c than or equal to a given constant s, where (d(i)), i = 1,
c ..., n are a given set of observation weights. the
c function determined is a spline under tension with third
c derivative discontinuities at (x(i)), i = 2,..., n-1. for
c actual computation of points on the curve it is necessary
c to call the function curv2.
c
c on input--
c
c   n is the number of values to be smoothed (n.ge.2).
c
c   x is an array of the n increasing abscissae of the
c   values to be smoothed.
c
c   y is an array of the n ordinates of the values to be
c   smoothed, (i. e. y(k) is the functional value
c   corresponding to x(k) ).
c
c   d is a parameter containing the observation weights.
c   this may either be an array of length n or a scalar
c   (interpreted as a constant). the value of d
c   corresponding to the observation (x(k),y(k)) should
c   be an approximation to the standard deviation of error.
c
c   isw contains a switch indicating whether the parameter
c   d is to be considered a vector or a scalar,
c          = 0 if d is an array of length n,
c          = 1 if d is a scalar.
c
c   s contains the value controlling the smoothing. this
c   must be non-negative. for s equal to zero, the
c   subroutine does interpolation, larger values lead to
c   smoother funtions. if parameter d contains standard
c   deviation estimates, a reasonable value for s is
c   float(n).
c
c   eps contains a tolerance on the relative precision to
c   which s is to be interpreted. this must be greater than
c   or equal to zero and less than equal or equal to one. a
c   reasonable value for eps is sqrt(2./float(n)).
c
c   ys is an array of length at least n.
c
c   ysp is an array of length at least n.
c
c   sigma contains the tension factor. this value indicates
c   the degree to which the first derivative part of the
c   smoothing functional is emphasized. if sigma is nearly
c   zero (e. g. .001) the resulting curve is approximately a
c   cubic spline. if sigma is large (e. g. 50.) the
c   resulting curve is nearly a polygonal line. if sigma
c   equals zero a cubic spline results. a standard value for
c   sigma is approximately 1.
c
c and
c
c   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, and v are
c   arrays of length at least n which are used for scratch
c   storage.
c
c on output--
c
c   ys contains the smoothed ordinate values.
c
c   ysp contains the values of the second derivative of the
c   smoothed curve at the given nodes.
c
c   ierr contains an error flag,
c        = 0 for normal return,
c        = 1 if n is less than 2,
c        = 2 if s is negative,
c        = 3 if eps is negative or greater than one,
c        = 4 if x-values are not strictly increasing,
c        = 5 if a d-value is non-positive.
c
c and
c
c   n, x, y, d, isw, s, eps, and sigma are unaltered.
c
c this subroutine references package modules terms and
c snhcsh.
c
c-----------------------------------------------------------
c
      if (n .lt. 2) go to 16
      if (s .lt. 0.) go to 17
      if (eps .lt. 0. .or. eps .gt. 1.) go to 18
      ierr = 0
      p = 0.
      v(1) = 0.
      v(n) = 0.
      ysp(1) = 0.
      ysp(n) = 0.
      if (n .eq. 2) go to 14
      rsd1(1) = 0.
      rd(1) = 0.
      rsd2(n) = 0.
      rdim1 = 0.
      yspim2 = 0.
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c form t matrix and second differences of y into ys
c
      nm1 = n-1
      nm3 = n-3
      delxi1 = 1.
      delyi1 = 0.
      dim1 = 0.
      do 1 i = 1,nm1
        delxi = x(i+1)-x(i)
        if (delxi .le. 0.) go to 19
        delyi = (y(i+1)-y(i))/delxi
        ys(i) = delyi-delyi1
        call terms (di,tsd1(i+1),sigmap,delxi)
        td(i) = di+dim1
        hd(i) = -(1./delxi+1./delxi1)
        hsd1(i+1) = 1./delxi
        delxi1 = delxi
        delyi1 = delyi
    1   dim1 = di
c
c calculate lower and upper tolerances
c
      sl = s*(1.-eps)
      su = s*(1.+eps)
      if (isw .eq. 1) go to 3
c
c form h matrix - d array
c
      if (d(1) .le. 0. .or. d(2) .le. 0.) go to 20
      betapp = 0.
      betap = 0.
      alphap = 0.
      do 2 i = 2,nm1
        alpha = hd(i)*d(i)*d(i)
        if (d(i+1) .le. 0.) go to 20
        beta = hsd1(i+1)*d(i+1)*d(i+1)
        hd(i) = (hsd1(i)*d(i-1))**2+alpha*hd(i)
     *                             +beta*hsd1(i+1)
        hsd2(i) = hsd1(i)*betapp
        hsd1(i) = hsd1(i)*(alpha+alphap)
        alphap = alpha
        betapp = betap
    2   betap = beta
      go to 5
c
c form h matrix - d constant
c
    3 if (d(1) .le. 0.) go to 20
      sl = d(1)*d(1)*sl
      su = d(1)*d(1)*su
      hsd1p = 0.
      hdim1 = 0.
      do 4 i = 2,nm1
        hdi = hd(i)
        hd(i) = hsd1(i)*hsd1(i)+hdi*hdi+hsd1(i+1)*hsd1(i+1)
        hsd2(i) = hsd1(i)*hsd1p
        hsd1p = hsd1(i)
        hsd1(i) = hsd1p*(hdi+hdim1)
    4   hdim1 = hdi
c
c top of iteration
c cholesky factorization of p*t+h into r
c
    5 do 6 i = 2,nm1
        rsd2i = hsd2(i)
        rsd1i = p*tsd1(i)+hsd1(i)-rsd2i*rsd1(i-1)
        rsd2(i) = rsd2i*rdim1
        rdim1 = rd(i-1)
        rsd1(i) = rsd1i*rdim1
        rd(i) = 1./(p*td(i)+hd(i)-rsd1i*rsd1(i)
     *                           -rsd2i*rsd2(i))
        ysp(i) = ys(i)-rsd1(i)*ysp(i-1)-rsd2(i)*yspim2
    6   yspim2 = ysp(i-1)
c
c back solve of r(transpose)* r * ysp = ys
c
      ysp(nm1) = rd(nm1)*ysp(nm1)
      if (n .eq. 3) go to 8
      do 7 ibak = 1,nm3
        i = nm1-ibak
    7   ysp(i) = rd(i)*ysp(i)-rsd1(i+1)*ysp(i+1)
     *                       -rsd2(i+2)*ysp(i+2)
    8 sum = 0.
      delyi1 = 0.
      if (isw .eq. 1) go to 10
c
c calculation of residual norm
c  - d array
c
      do 9 i = 1,nm1
        delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
        v(i) = (delyi-delyi1)*d(i)*d(i)
        sum = sum+v(i)*(delyi-delyi1)
    9   delyi1 = delyi
      v(n) = -delyi1*d(n)*d(n)
      go to 12
c
c calculation of residual norm
c  - d constant
c
   10 do 11 i = 1,nm1
        delyi = (ysp(i+1)-ysp(i))/(x(i+1)-x(i))
        v(i) = delyi-delyi1
        sum = sum+v(i)*(delyi-delyi1)
   11   delyi1 = delyi
      v(n) = -delyi1
   12 sum = sum-v(n)*delyi1
c
c test for convergence
c
      if (sum .le. su) go to 14
c
c calculation of newton correction
c
      f = 0.
      g = 0.
      wim2 = 0.
      wim1 = 0.
      do 13 i = 2,nm1
        tui = tsd1(i)*ysp(i-1)+td(i)*ysp(i)
     *                        +tsd1(i+1)*ysp(i+1)
        wi = tui-rsd1(i)*wim1-rsd2(i)*wim2
        f = f+tui*ysp(i)
        g = g+wi*wi*rd(i)
        wim2 = wim1
   13   wim1 = wi
      h = f-p*g
      if (h .le. 0.) go to 14
c
c update p - newton step
c
      step = (sum-sqrt(sum*sl))/h
      if (sl .ne. 0.) step = step*sqrt(sum/sl)
      p = p+step
      go to 5
c
c store smoothed y-values and second derivatives
c
   14 do 15 i = 1,n
        ys(i) = y(i)-v(i)
   15   ysp(i) = p*ysp(i)
      return
c
c n less than 2
c
   16 ierr = 1
      return
c
c s negative
c
   17 ierr = 2
      return
c
c eps negative or greater than 1
c
   18 ierr = 3
      return
c
c x-values not strictly increasing
c
   19 ierr = 4
      return
c
c weight non-positive
c
   20 ierr = 5
      return
      end
      function intrvl (t,x,n)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer n
      real t,x(n)
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function determines the index of the interval
c (determined by a given increasing sequence) in which
c a given value lies.
c
c on input--
c
c   t is the given value.
c
c   x is a vector of strictly increasing values.
c
c and
c
c   n is the length of x (n .ge. 2).
c
c on output--
c
c   intrvl returns an integer i such that
c
c          i =  1       if         e   t .lt. x(2)  ,
c          i =  n-1     if x(n-1) .le. t            ,
c          otherwise       x(i)  .le. t .le. x(i+1),
c
c none of the input parameters are altered.
c
c-----------------------------------------------------------
c
      save i
      data i /1/
c
      tt = t
c
c check for illegal i
c
      if (i .ge. n) i = n/2
c
c check old interval and extremes
c
      if (tt .lt. x(i)) then
        if (tt .le. x(2)) then
          i = 1
          intrvl = 1
          return
        else
          il = 2
          ih = i
        end if
      else if (tt .le. x(i+1)) then
        intrvl = i
        return
      else if (tt .ge. x(n-1)) then
        i = n-1
        intrvl = n-1
        return
      else
        il = i+1
        ih = n-1
      end if
c
c binary search loop
c
    1 i = (il+ih)/2
      if (tt .lt. x(i)) then
         ih = i
      else if (tt .gt. x(i+1)) then
         il = i+1
      else
         intrvl = i
         return
      end if
      go to 1
      end



      function intrvp (t,x,n,p,tp)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer n
      real t,x(n),p,tp
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function determines the index of the interval
c (determined by a given increasing sequence) in which a
c given value lies, after translating the value to within
c the correct period.  it also returns this translated value.
c
c on input--
c
c   t is the given value.
c
c   x is a vector of strictly increasing values.
c
c   n is the length of x (n .ge. 2).
c
c and
c
c   p contains the period.
c
c on output--
c
c   tp contains a translated value of t (i. e. x(1) .le. tp,
c   tp .lt. x(1)+p, and tp = t + k*p for some integer k).
c
c   intrvl returns an integer i such that
c
c          i = 1       if             tp .lt. x(2)  ,
c          i = n       if   x(n) .le. tp            ,
c          otherwise       x(i)  .le. tp .lt. x(i+1),
c
c none of the input parameters are altered.
c
c-----------------------------------------------------------
c
      save i
      data i /1/
c
      nper = (t-x(1))/p
      tp = t-float(nper)*p
      if (tp .lt. x(1)) tp = tp+p
      tt = tp
c
c check for illegal i
c
      if (i .ge. n) i = n/2
c
c check old interval and extremes
c
      if (tt .lt. x(i)) then
        if (tt .le. x(2)) then
          i = 1
          intrvp = 1
          return
        else
          il = 2
          ih = i
        end if
      else if (tt .le. x(i+1)) then
        intrvp = i
        return
      else if (tt .ge. x(n)) then
        i = n
        intrvp = n
        return
      else
        il = i+1
        ih = n
      end if
c
c binary search loop
c
    1 i = (il+ih)/2
      if (tt .lt. x(i)) then
         ih = i
      else if (tt .gt. x(i+1)) then
         il = i+1
      else
         intrvp = i
         return
      end if
      go to 1
      end



c$$$      subroutine snhcsh (sinhm,coshm,x,isw)
c$$$c
c$$$      integer isw
c$$$      real sinhm,coshm,x
c$$$c
c$$$c                                 coded by alan kaylor cline
c$$$c                           from fitpack -- january 26, 1987
c$$$c                        a curve and surface fitting package
c$$$c                      a product of pleasant valley software
c$$$c                  8603 altus cove, austin, texas 78759, usa
c$$$c
c$$$c this subroutine returns approximations to
c$$$c       sinhm(x) = sinh(x)/x-1
c$$$c       coshm(x) = cosh(x)-1
c$$$c and
c$$$c       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
c$$$c with relative error less than 1.0e-6
c$$$c
c$$$c on input--
c$$$c
c$$$c   x contains the value of the independent variable.
c$$$c
c$$$c   isw indicates the function desired
c$$$c           = -1 if only sinhm is desired,
c$$$c           =  0 if both sinhm and coshm are desired,
c$$$c           =  1 if only coshm is desired,
c$$$c           =  2 if only coshmm is desired,
c$$$c           =  3 if both sinhm and coshmm are desired.
c$$$c
c$$$c on output--
c$$$c
c$$$c   sinhm contains the value of sinhm(x) if isw .le. 0 or
c$$$c   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
c$$$c   2).
c$$$c
c$$$c   coshm contains the value of coshm(x) if isw .eq. 0 or
c$$$c   isw .eq. 1 and contains the value of coshmm(x) if isw
c$$$c   .ge. 2 (coshm is unaltered if isw .eq. -1).
c$$$c
c$$$c and
c$$$c
c$$$c   x and isw are unaltered.
c$$$c
c$$$c-----------------------------------------------------------
c$$$c
c$$$      data sp13/.3029390e-5/,
c$$$     *     sp12/.1975135e-3/,
c$$$     *     sp11/.8334261e-2/,
c$$$     *     sp10/.1666665e0/
c$$$      data sp24/.3693467e-7/,
c$$$     *     sp23/.2459974e-5/,
c$$$     *     sp22/.2018107e-3/,
c$$$     *     sp21/.8315072e-2/,
c$$$     *     sp20/.1667035e0/
c$$$      data sp33/.6666558e-5/,
c$$$     *     sp32/.6646307e-3/,
c$$$     *     sp31/.4001477e-1/,
c$$$     *     sq32/.2037930e-3/,
c$$$     *     sq31/-.6372739e-1/,
c$$$     *     sq30/.6017497e1/
c$$$      data sp43/.2311816e-4/,
c$$$     *     sp42/.2729702e-3/,
c$$$     *     sp41/.9868757e-1/,
c$$$     *     sq42/.1776637e-3/,
c$$$     *     sq41/-.7549779e-1/,
c$$$     *     sq40/.9110034e1/
c$$$      data cp4/.2982628e-6/,
c$$$     *     cp3/.2472673e-4/,
c$$$     *     cp2/.1388967e-2/,
c$$$     *     cp1/.4166665e-1/,
c$$$     *     cp0/.5000000e0/
c$$$c
c$$$      ax = abs(x)
c$$$      if (isw .ge. 0) go to 5
c$$$c
c$$$c sinhm approximation
c$$$c
c$$$      if (ax .gt. 4.45) go to 2
c$$$      xs = ax*ax
c$$$      if (ax .gt. 2.3) go to 1
c$$$c
c$$$c sinhm approximation on (0.,2.3)
c$$$c
c$$$      sinhm = xs*(((sp13*xs+sp12)*xs+sp11)*xs+sp10)
c$$$      return
c$$$c
c$$$c sinhm approximation on (2.3,4.45)
c$$$c
c$$$    1 sinhm = xs*((((sp24*xs+sp23)*xs+sp22)*xs+sp21)
c$$$     .               *xs+sp20)
c$$$      return
c$$$    2 if (ax .gt. 7.65) go to 3
c$$$c
c$$$c sinhm approximation on (4.45,7.65)
c$$$c
c$$$      xs = ax*ax
c$$$      sinhm = xs*(((sp33*xs+sp32)*xs+sp31)*xs+1.)/
c$$$     .             ((sq32*xs+sq31)*xs+sq30)
c$$$      return
c$$$    3 if (ax .gt. 10.1) go to 4
c$$$c
c$$$c sinhm approximation on (7.65,10.1)
c$$$c
c$$$      xs = ax*ax
c$$$      sinhm = xs*(((sp43*xs+sp42)*xs+sp41)*xs+1.)/
c$$$     .             ((sq42*xs+sq41)*xs+sq40)
c$$$      return
c$$$c
c$$$c sinhm approximation above 10.1
c$$$c
c$$$    4 sinhm = exp(ax)/(ax+ax)-1.
c$$$      return
c$$$c
c$$$c coshm and (possibly) sinhm approximation
c$$$c
c$$$    5 if (isw .ge. 2) go to 7
c$$$      if (ax .gt. 2.3) go to 6
c$$$      xs = ax*ax
c$$$      coshm = xs*((((cp4*xs+cp3)*xs+cp2)*xs+cp1)*xs+cp0)
c$$$      if (isw .eq. 0) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)
c$$$     .                              *xs+sp10)
c$$$      return
c$$$    6 expx = exp(ax)
c$$$      coshm = (expx+1./expx)/2.-1.
c$$$      if (isw .eq. 0) sinhm = (expx-1./expx)/(ax+ax)-1.
c$$$      return
c$$$c
c$$$c coshmm and (possibly) sinhm approximation
c$$$c
c$$$    7 xs = ax*ax
c$$$      if (ax .gt. 2.3) go to 8
c$$$      coshm = xs*(((cp4*xs+cp3)*xs+cp2)*xs+cp1)
c$$$      if (isw .eq. 3) sinhm = xs*(((sp13*xs+sp12)*xs+sp11)
c$$$     .                              *xs+sp10)
c$$$      return
c$$$    8 expx = exp(ax)
c$$$      coshm = ((expx+1./expx-xs)/2.-1.)/xs
c$$$      if (isw .eq. 3) sinhm = (expx-1./expx)/(ax+ax)-1.
c$$$      return
c$$$      end



      subroutine snhcsh (sinhm,coshm,x,isw)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      integer isw
      real sinhm,coshm,x
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine returns approximations to
c       sinhm(x) = sinh(x)/x-1
c       coshm(x) = cosh(x)-1
c and
c       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
c with relative error less than 4.0e-14.
c
c on input--
c
c   x contains the value of the independent variable.
c
c   isw indicates the function desired
c           = -1 if only sinhm is desired,
c           =  0 if both sinhm and coshm are desired,
c           =  1 if only coshm is desired,
c           =  2 if only coshmm is desired,
c           =  3 if both sinhm and coshmm are desired.
c
c on output--
c
c   sinhm contains the value of sinhm(x) if isw .le. 0 or
c   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
c   2).
c
c   coshm contains the value of coshm(x) if isw .eq. 0 or
c   isw .eq. 1 and contains the value of coshmm(x) if isw
c   .ge. 2 (coshm is unaltered if isw .eq. -1).
c
c and
c
c   x and isw are unaltered.
c
c-----------------------------------------------------------
c
      data sp14/.227581660976348e-7/,
     *     sp13/.612189863171694e-5/,
     *     sp12/.715314759211209e-3/,
     *     sp11/.398088289992973e-1/,
     *     sq12/.206382701413725e-3/,
     *     sq11/-.611470260009508e-1/,
     *     sq10/.599999999999986e+1/
      data sp25/.129094158037272e-9/,
     *     sp24/.473731823101666e-7/,
     *     sp23/.849213455598455e-5/,
     *     sp22/.833264803327242e-3/,
     *     sp21/.425024142813226e-1/,
     *     sq22/.106008515744821e-3/,
     *     sq21/-.449855169512505e-1/,
     *     sq20/.600000000268619e+1/
      data sp35/.155193945864942e-9/,
     *     sp34/.511529451668737e-7/,
     *     sp33/.884775635776784e-5/,
     *     sp32/.850447617691392e-3/,
     *     sp31/.428888148791777e-1/,
     *     sq32/.933128831061610e-4/,
     *     sq31/-.426677570538507e-1/,
     *     sq30/.600000145086489e+1/
      data sp45/.188070632058331e-9/,
     *     sp44/.545792817714192e-7/,
     *     sp43/.920119535795222e-5/,
     *     sp42/.866559391672985e-3/,
     *     sp41/.432535234960858e-1/,
     *     sq42/.824891748820670e-4/,
     *     sq41/-.404938841672262e-1/,
     *     sq40/.600005006283834e+1/
      data cp5/.552200614584744e-9/,
     *     cp4/.181666923620944e-6/,
     *     cp3/.270540125846525e-4/,
     *     cp2/.206270719503934e-2/,
     *     cp1/.744437205569040e-1/,
     *     cq2/.514609638642689e-4/,
     *     cq1/-.177792255528382e-1/,
     *     cq0/.200000000000000e+1/
      data zp4/.664418805876835e-8/,
     *     zp3/.218274535686385e-5/,
     *     zp2/.324851059327161e-3/,
     *     zp1/.244515150174258e-1/,
     *     zq2/.616165782306621e-3/,
     *     zq1/-.213163639579425e0/,
     *     zq0/.240000000000000e+2/
c
      ax = abs(x)
      if (isw .ge. 0) go to 5
c
c sinhm approximation
c
      if (ax .gt. 3.9) go to 2
      xs = ax*ax
      if (ax .gt. 2.2) go to 1
c
c sinhm approximation on (0.,2.2)
c
      sinhm = xs*((((sp14*xs+sp13)*xs+sp12)*xs+sp11)*xs+1.)/
     .             ((sq12*xs+sq11)*xs+sq10)
      return
c
c sinhm approximation on (2.2,3.9)
c
    1 sinhm = xs*(((((sp25*xs+sp24)*xs+sp23)*xs+sp22)*xs+sp21)
     .        *xs+1.)/((sq22*xs+sq21)*xs+sq20)
      return
    2 if (ax .gt. 5.1) go to 3
c
c sinhm approximation on (3.9,5.1)
c
      xs = ax*ax
      sinhm = xs*(((((sp35*xs+sp34)*xs+sp33)*xs+sp32)*xs+sp31)
     .        *xs+1.)/((sq32*xs+sq31)*xs+sq30)
      return
    3 if (ax .gt. 6.1) go to 4
c
c sinhm approximation on (5.1,6.1)
c
      xs = ax*ax
      sinhm = xs*(((((sp45*xs+sp44)*xs+sp43)*xs+sp42)*xs+sp41)
     .        *xs+1.)/((sq42*xs+sq41)*xs+sq40)
      return
c
c sinhm approximation above 6.1
c
    4 expx = exp(ax)
      sinhm = (expx-1./expx)/(ax+ax)-1.
      return
c
c coshm and (possibly) sinhm approximation
c
    5 if (isw .ge. 2) go to 7
      if (ax .gt. 2.2) go to 6
      xs = ax*ax
      coshm = xs*(((((cp5*xs+cp4)*xs+cp3)*xs+cp2)*xs+cp1)
     .        *xs+1.)/((cq2*xs+cq1)*xs+cq0)
      if (isw .eq. 0) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)
     .          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)
      return
    6 expx = exp(ax)
      coshm = (expx+1./expx)/2.-1.
      if (isw .eq. 0) sinhm = (expx-1./expx)/(ax+ax)-1.
      return
c
c coshmm and (possibly) sinhm approximation
c
    7 xs = ax*ax
      if (ax .gt. 2.2) go to 8
      coshm = xs*((((zp4*xs+zp3)*xs+zp2)*xs+zp1)*xs+1.)/
     .             ((zq2*xs+zq1)*xs+zq0)
      if (isw .eq. 3) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)
     .          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)
      return
    8 expx = exp(ax)
      coshm = ((expx+1./expx-xs)/2.-1.)/xs
      if (isw .eq. 3) sinhm = (expx-1./expx)/(ax+ax)-1.
      return
      end



      subroutine terms (diag,sdiag,sigma,del)
c
      IMPLICIT REAL(A-H,O-Z), INTEGER(I-N)
      real diag,sdiag,sigma,del
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine computes the diagonal and superdiagonal
c terms of the tridiagonal linear system associated with
c spline under tension interpolation.
c
c on input--
c
c   sigma contains the tension factor.
c
c and
c
c   del contains the step size.
c
c on output--
c
c                sigma*del*cosh(sigma*del) - sinh(sigma*del)
c   diag = del*--------------------------------------------.
c                     (sigma*del)**2 * sinh(sigma*del)
c
c                   sinh(sigma*del) - sigma*del
c   sdiag = del*----------------------------------.
c                (sigma*del)**2 * sinh(sigma*del)
c
c and
c
c   sigma and del are unaltered.
c
c this subroutine references package module snhcsh.
c
c-----------------------------------------------------------
c
      if (sigma .ne. 0.) go to 1
      diag = del/3.
      sdiag = del/6.
      return
    1 sigdel = sigma*del
      call snhcsh (sinhm,coshm,sigdel,0)
      denom = sigma*sigdel*(1.+sinhm)
      diag = (coshm-sinhm)/denom
      sdiag = sinhm/denom
      return
      end
