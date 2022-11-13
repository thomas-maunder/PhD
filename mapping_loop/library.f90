SUBROUTINE mapping(nr, nth, nph, yywt, nx, ny, nz, nnuc, x0, y0, &
  z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew, &
  muc, mul, mur, phic, phil, phir, tol)

  USE, INTRINSIC :: IEEE_ARITHMETIC

  IMPLICIT NONE

  real*8, parameter :: x12 = 1.d0/12.d0

  INTEGER, INTENT(in)  ::  nr, nth, nph, nx, ny, nz, nnuc
  REAL*8, INTENT(in)  ::  x0, y0, z0, tol
  REAL*8, INTENT(in)  ::  yywt(0:nth-1,0:nph-1)
  REAL*8, INTENT(in)  ::  dphi, dx, dy, dz
  REAL*8, INTENT(in)  ::  r3c(0:nr-1), r3l(0:nr-1), r3r(0:nr-1)
  REAL*8, INTENT(in)  ::  muc(0:nth-1), mul(0:nth-1), mur(0:nth-1)
  REAL*8, INTENT(in)  ::  phic(0:nph-1), phil(0:nph-1),  phir(0:nph-1)
  REAL*8, INTENT(in)  ::  rho(0:nr-1, 0:nth-1, 0:nph-1)
  REAL*8, INTENT(in)  ::  xnu(0:nr-1, 0:nth-1, 0:nph-1, 0:nnuc-1)
  REAL*8, INTENT(inout)  ::  rhonew(0:nx-1, 0:ny-1, 0:nz-1)
  REAL*8, INTENT(inout)  ::  xnunew(0:nx-1, 0:ny-1, 0:nz-1, 0:nnuc-1)
  REAL*8, PARAMETER  ::  pi = 4.d0 * ATAN(1.d0)
  ! REAL*8 :: mapped_mass = 0.0
  ! REAL*8 :: input_mass = 0.0
  INTEGER  ::  kk, km1, kp1, jj, jm1, jp1, ii, im1, ip1, nsub, nsubmax, flipyy
  REAL*8  ::  sr, sl, s1rho, s2rho, s3rho, r3c0, muc0, phic0, dv, dv0, &
    rho0, rho0i, dphii, dxi, dyi, dzi, xoffs, yoffs, zoffs
  REAL*8  ::  slxnu(0:nnuc-1), srxnu(0:nnuc-1), s1xnu(0:nnuc-1), &
    s2xnu(0:nnuc-1), s3xnu(0:nnuc-1), xnu0(0:nnuc-1), &
    delta_xnu0(0:nnuc-1)

  print*,'Using tolerance',tol

  ! print*,nx,ny,nz
  print*,dx,dy,dz

  flipyy = 0

  dphii = 1.d0 / dphi

  dxi = 1.d0/dx
  dyi = 1.d0/dy
  dzi = 1.d0/dz

  xoffs = - x0 * dxi
  yoffs = - y0 * dyi
  zoffs = - z0 * dzi

   DO kk = 0,nph-1
     print*,'kk = ',kk

     km1 = max(kk-1,0)
     kp1 = min(kk+1,nph-1)

     theta_loop: DO jj = 0,nth-1
      ! print*,'jj = ',jj

       jm1 = max(jj-1,0)
       jp1 = min(jj+1,nth-1)

       IF (nph > 1) THEN
         IF (jj == nth/2-1) THEN
           jp1 = jj
         ELSE IF (jj == nth/2) THEN
           jm1 = jj
         ELSE IF (jj >= nth/2) THEN
           flipyy = 1
         ELSE
           flipyy = 0
         END IF
       END IF

       IF (yywt(jj,kk) == 0) CYCLE theta_loop

       DO ii = 0,nr-1
        ! print*,'ii = ',ii

          !IF (nph==1) THEN
          !  dv = abs(Xr-Xl) * abs(Yr-Yl) * dphi
          !ELSE
            dv = abs(r3r(ii) - r3l(ii)) * abs(mur(jj)-mul(jj)) * abs(phir(kk)-phil(kk)) * yywt(jj,kk)
          !END IF
          ! input_mass = input_mass + rho(ii,jj,kk) * dv

         im1 = max(ii-1,0)
         ip1 = min(ii+1,nr-1)

         if ((ii==0).or.(ii==nr-1)) then
            s1rho = 0.d0
            s1xnu = 0.d0
         else
            sr = (rho(ip1, jj, kk) - rho(ii, jj, kk)) /(r3c(ip1) - r3c(ii))
            sl = (rho(ii, jj, kk) - rho(im1, jj, kk)) /(r3c(ii) - r3c(im1))
            s1rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)

            srxnu = (xnu(ip1, jj, kk, :) - xnu(ii, jj, kk, :)) / (r3c(ip1) - r3c(ii))
            slxnu = (xnu(ii, jj, kk, :) - xnu(im1, jj, kk, :)) / (r3c(ii) - r3c(im1))
            s1xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)
         endif

         if ((jj==0).or.(jj==nth-1)) then
            s2rho = 0.d0
            s2xnu = 0.d0
         else
            sr = (rho(ii, jp1, kk) - rho(ii, jj, kk)) /(muc(jp1) - muc(jj))
            sl = (rho(ii, jj, kk) - rho(ii, jm1, kk)) /(muc(jj) - muc(jm1))
            s2rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)

            srxnu = (xnu(ii, jp1, kk, :) - xnu(ii, jj, kk, :)) / (muc(jp1) - muc(jj))
            slxnu = (xnu(ii, jj, kk, :) - xnu(ii, jm1, kk, :)) / (muc(jj) - muc(jm1))
            s2xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)
         endif

         sr = (rho(ii, jj, kp1) - rho(ii, jj, kk)) * dphii
         sl = (rho(ii, jj, kk) - rho(ii, jj, km1)) * dphii
         s3rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)

         srxnu = (xnu(ii, jj, kp1, :) - xnu(ii, jj, kk, :)) * dphii
         slxnu = (xnu(ii, jj, kk, :) - xnu(ii, jj, km1, :)) * dphii
         s3xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)

         rho0 = rho(ii, jj, kk)
         xnu0 = xnu(ii, jj, kk, :)

         rho0i = 1.d0 / rho0

         s1xnu = s1xnu / (1.d0+ abs(s1rho * (r3r(ii) - r3l(ii)) * rho0i) * x12)
         s2xnu = s2xnu / (1.d0+ abs(s2rho * (mur(jj) - mul(jj)) * rho0i) * x12)
         s3xnu = s3xnu / (1.d0+ abs(s3rho * (phir(kk) - phil(kk)) * rho0i) * x12)


	 !s1rho = 0.d0
	 !s2rho = 0.d0
	 !s3rho = 0.d0
	 !s1xnu = 0.d0
	 !s2xnu = 0.d0
	 !s3xnu = 0.d0

         delta_xnu0 = (s1rho * s1xnu * (r3r(ii) - r3l(ii))**2 &
                       + s2rho * s2xnu * abs(mur(jj) - mul(jj))**2 &
                       + s3rho * s3xnu * abs(phir(kk) - phil(kk))**2) * x12
         xnu0 = xnu0 - delta_xnu0 * rho0i

         r3c0 = r3c(ii)
         muc0 = muc(jj)
         phic0 = phic(kk)

         nsub = 1
         nsubmax = 1

         call map_cube(r3l(ii), r3r(ii), mul(jj), mur(jj), phil(kk), phir(kk))

       END DO
     END DO theta_loop
  END DO

  ! PRINT*,'maped_mass',mapped_mass,mapped_mass/1.99d33
  ! PRINT*,'input_mass',input_mass,input_mass/1.99d33


 CONTAINS
   ELEMENTAL REAL*8 FUNCTION heaviside(x1)
     IMPLICIT NONE

     REAL*8, INTENT(in)  ::  x1

     IF (x1 < 0) THEN
       heaviside = 0.d0
     ELSE
       heaviside = 1.d0
     END IF

     RETURN

   END FUNCTION heaviside

   FUNCTION to_cart(R3, Mu, Phi, flipyy) result(triple)
     IMPLICIT NONE

     real*8, parameter :: p13 = 1.d0/3.d0

     REAL*8, INTENT(in)  ::  R3, Mu, Phi
     INTEGER, INTENT(in)  ::  flipyy
     REAL*8, DIMENSION(3)  ::  triple
     REAL*8  ::  x, y, z, tempx, tempy, tempz, R, Theta, s
     REAL*8  ::  n1, n2, n3, cx, cy, cz, tempx2, tempy2, tempz2
     REAL*8, DIMENSION(3)  :: cross_prod

     R = (3.d0*R3)**p13
  	 Theta = acos(Mu)
     s = sin(Theta)

     IF (ny>1 .and. nz >1) THEN !3d
       tempx = R * s * cos(Phi)
       tempy = R * s * sin(Phi)
       tempz = R * Mu
       IF (flipyy == 0) THEN
         x = tempx
         y = tempy
         z = tempz
       ELSE
         x = -tempx
         y = tempz
         z = tempy
       END IF
    !  ELSEIF (ny>1 .and. nz>1 .and. nr>1 .and. nth==1 .and. nph>1) THEN !3d to 2d
    !    n1 = 1.d0
    !    n2 = 0.d0
    !    n3 = 0.d0
    !    tempx = R * s * cos(Phi)
    !    tempy = R * s * sin(Phi)
    !    tempz = R * cos(Theta)
    !    IF (flipyy == 0) THEN
    !      tempx2 = tempx
    !      tempy2 = tempy
    !      tempz2 = tempz
    !    ELSE
    !      tempx2 = -tempx
    !      tempy2 = tempz
    !      tempz2 = tempy
    !    END IF
       cx = tempy2 * n3 - tempz2 * n2
       cy = tempz2 * n1 - tempx2 * n3
       cz = tempx2 * n2 - tempy2 * n1
       cross_prod(:) = (/ cx, cy, cz /)
       x = NORM2(cross_prod)
       y = 0.d0
       z = DOT_PRODUCT((/ tempx2, tempy2, tempz2 /), (/ n1, n2, n3 /))
     ELSEIF (ny==1 .and. nz >1) THEN !2d
       x = R * s
       y = 0.
       z = R * Mu
     ELSE !1d
       x = R
       y = 0.0
       z = 0.0
     ENDIF

    !  n1 = 1.d0
    !  n2 = 0.d0
    !  n3 = 0.d0
    !  tempx = R * s * cos(Phi)
    !  tempy = R * s * sin(Phi)
    !  tempz = R * cos(Theta)
    !  cx = tempy * n3 - tempz * n2
    !  cy = tempz * n1 - tempx * n3
    !  cz = tempx * n2 - tempy * n1
    !  cross_prod(:) = (/ cx, cy, cz /)
    !  x = NORM2(cross_prod)
    !  y = 0.d0
    !  z = DOT_PRODUCT((/ tempx, tempy, tempz /), (/ n1, n2, n3 /))
    
     

     triple(:) = (/ x, y, z /)

   END FUNCTION to_cart

   ! FUNCTION target_indices(vec, x0, y0, z0, dx, dy, dz, m, n, o)
   FUNCTION target_indices(vec) RESULT(triple)
     IMPLICIT NONE

     REAL*8, DIMENSION(3), INTENT(in)  ::  vec
     integer, DIMENSION(3)  ::  triple
     REAL*8  ::  x, y, z

     x = vec(1)
     y = vec(2)
     z = vec(3)

     triple(1) = floor(x*dxi + xoffs)
     triple(2) = floor(y*dyi + yoffs)
     triple(3) = floor(z*dzi + zoffs)
     IF (ny.eq.1) triple(2)=0
     IF (nz.eq.1) triple(3)=0

   END FUNCTION target_indices


   RECURSIVE SUBROUTINE map_cube(Xl, Xr, Yl, Yr, Zl, Zr)
     USE, INTRINSIC :: IEEE_ARITHMETIC

     IMPLICIT NONE

     REAL*8, intent(in)  ::  Xl, Xr, Yl, Yr, Zl, Zr

     integer, DIMENSION(3,8)  ::  c
     integer, DIMENSION(3)  ::  res
     REAL*8  ::  Xc, Yc, Zc, temp
     INTEGER  ::  m0, m1, n0, n1, o0, o1, mc, nc, oc, i

     c(:,1) = target_indices(to_cart(Xl, Yl, Zl, flipyy))
     c(:,2) = target_indices(to_cart(Xl, Yl, Zr, flipyy))
     c(:,3) = target_indices(to_cart(Xl, Yr, Zl, flipyy))
     c(:,4) = target_indices(to_cart(Xl, Yr, Zr, flipyy))
     c(:,5) = target_indices(to_cart(Xr, Yl, Zl, flipyy))
     c(:,6) = target_indices(to_cart(Xr, Yl, Zr, flipyy))
     c(:,7) = target_indices(to_cart(Xr, Yr, Zl, flipyy))
     c(:,8) = target_indices(to_cart(Xr, Yr, Zr, flipyy))

!    print*,c1,c2,c3,c4,c5,c6,c7,c8

     m0 = minval(c(1,:))
     m1 = maxval(c(1,:))
     n0 = minval(c(2,:))
     n1 = maxval(c(2,:))
     o0 = minval(c(3,:))
     o1 = maxval(c(3,:))

!    print*,'mno',m0,m1,n0,n1,o0,o1

     Xc = 0.5d0*(Xr+Xl)
     Yc = 0.5d0*(Yr+Yl)
     Zc = 0.5d0*(Zr+Zl)
     res = target_indices(to_cart(Xc, Yc, Zc, flipyy))
     mc = res(1)
     nc = res(2)
     oc = res(3)

!    print*,Xc,Yc,Zc,mc,nc,oc

     IF (nph==1) THEN
       dv = abs(Xr-Xl) * abs(Yr-Yl) * dphi
     ELSE
       dv = abs(Xr-Xl) * abs(Yr-Yl) * abs(Zr-Zl) * yywt(jj,kk)
     END IF
!    print*,dv
     IF (ny==1) THEN
       IF (nz==1) THEN
         dv0 = (4./3.)*pi*(3*mc**2+3*mc+1.0)*dx**3
       ELSE
         dv0 = pi*(2*mc+1)*dx**2*dz
       END IF
     ELSE
       ! dv0 = float(dx*dy*dz)
       dv0 = dx*dy*dz
     END IF
!    print*,dv0
     IF ((max(m0,m1)>=nx) .OR. (max(n0,n1)>=ny) .OR. (max(o0,o1)>=nz)) THEN
     !   print*,'outside ARTIS grid (a)',m0,m1,n0,n1,o0,o1,(r3c0*3)**(1.0/3.0), &
     !    muc0,phic0
       RETURN
     END IF
     IF ((min(m0,m1)<0) .OR. (min(n0,n1)<0) .OR. (min(o0,o1)<0)) THEN
     !   print*,'outside ARTIS grid (b)',m0,m1,n0,n1,o0,o1,(r3c0*3)**(1.0/3.0), &
     !    muc0,phic0
       RETURN
     END IF

     IF ( ((m0==m1) .AND. (n0==n1) .AND. (o0==o1)) .OR. (dv < tol*dv0) ) THEN
        temp = (rho0+s1rho*(Xc-r3c0)+s2rho*(Yc-muc0) + s3rho*(Zc-phic0)) * dv/dv0
        rhonew(mc,nc,oc) = rhonew(mc,nc,oc) + temp
        xnunew(mc,nc,oc,:) = xnunew(mc,nc,oc,:) + temp * &
             (xnu0(:)+s1xnu(:)*(Xc-r3c0)+s2xnu(:)*(Yc-muc0)+s3xnu(:)*(Zc-phic0))
        ! print*,'map',(r3c0*3)**(1.0/3.0),muc0,phic0,mc,nc,oc,dx*(mc+0.5),dv,dv0
        ! mapped_mass = mapped_mass + temp * dv0
       RETURN
    endif

    nsub = nsub + 1
    nsubmax = max(nsub,nsubmax)
    call map_cube(Xl, Xc, Yl, Yc, Zl, Zc)
    call map_cube(Xl, Xc, Yc, Yr, Zl, Zc)
    call map_cube(Xc, Xr, Yl, Yc, Zl, Zc)
    call map_cube(Xc, Xr, Yc, Yr, Zl, Zc)
    IF (nph > 1) THEN
       call map_cube(Xl, Xc, Yl, Yc, Zc, Zr)
       call map_cube(Xl, Xc, Yc, Yr, Zc, Zr)
       call map_cube(Xc, Xr, Yl, Yc, Zc, Zr)
       call map_cube(Xc, Xr, Yc, Yr, Zc, Zr)
    END IF
    nsub = nsub - 1


   END SUBROUTINE map_cube
END SUBROUTINE mapping
