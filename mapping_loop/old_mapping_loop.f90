! PROGRAM test
! 	IMPLICIT NONE
! 	print*, "6.8968486193938672 0. -4.24297576045549740e-2"
! END PROGRAM test

!SUBROUTINE mapping(nr, nth, nph, yywt_, nx, ny, nz, nnuc, x0, y0, &
! 	z0, rho_, xnu_, r3c_, r3l_, r3r_, dx, dy, dz, dphi, rhonew_, xnunew_, &
! 	muc_, mul_, mur_, phic_, phil_, phir_)
SUBROUTINE mapping(nr, nth, nph, yywt, nx, ny, nz, nnuc, x0, y0, &
	z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew, &
	muc, mul, mur, phic, phil, phir)
	IMPLICIT NONE

	!INTEGER, INTENT(in)	:: nr, nth, nph
	INTEGER, INTENT(in)	::	nr, nth, nph, nx, ny, nz, nnuc
	REAL*8, INTENT(in)	::	x0, y0, z0
	REAL*8, INTENT(in)	::	yywt(0:nth-1,0:nph-1)
	!REAL*8, INTENT(in), TARGET	::	r3c_(:), r3l_(:), r3r_(:)
	!REAL*8, INTENT(in), TARGET	::	muc_(:), mul_(:), mur_(:)
	!REAL*8, INTENT(in), TARGET	::	phic_(:), phil_(:),	phir_(:)
	!REAL*8, INTENT(in), TARGET	::	rho_(:, :, :)
	!REAL*8, INTENT(in), TARGET	::	xnu_(:,:,:,:)
	!REAL*8, INTENT(out), TARGET	::	rhonew_(nx,ny,nz)
	!REAL*8, INTENT(out), TARGET	::	xnunew_(nx,ny,nz,nnuc)
	INTEGER	::	kk, km1, kp1, jj, jm1, jp1, ii, im1, ip1, nsub, nsubmax, flipyy
	REAL*8	::	sr, sl, s1rho, s2rho, s3rho, r3c0, muc0, phic0, tol, dv, dv0, &
		rho0
	REAL*8	::	slxnu(0:nnuc-1), srxnu(0:nnuc-1), s1xnu(0:nnuc-1), &
		s2xnu(0:nnuc-1), s3xnu(0:nnuc-1), xnu0(0:nnuc-1), &
		delta_xnu0(0:nnuc-1)
	REAL*8, INTENT(in)	::	dphi, dx, dy, dz
!	REAL*8	::	r3c(0:nr-1)
	REAL*8, INTENT(in)	::	r3c(0:nr-1), r3l(0:nr-1), r3r(0:nr-1)
	REAL*8, INTENT(in)	::	muc(0:nth-1), mul(0:nth-1), mur(0:nth-1)
	REAL*8, INTENT(in)	::	phic(0:nph-1), phil(0:nph-1),	phir(0:nph-1)
	REAL*8, INTENT(in)	::	rho(0:nr-1, 0:nth-1, 0:nph-1)
	REAL*8, INTENT(in)	::	xnu(0:nr-1, 0:nth-1, 0:nph-1, 0:nnuc-1)
	REAL*8, INTENT(inout)	::	rhonew(0:nx-1, 0:ny-1, 0:nz-1)
	REAL*8, INTENT(inout)	::	xnunew(0:nr-1, 0:nth-1, 0:nph-1, 0:nnuc-1)
	! INTEGER	::	kk, km1, kp1, jj, jm1, jp1, ii, im1, ip1, nsub, nsubmax, flipyy
	! REAL*8	::	sr, sl, s1rho, s2rho, s3rho, r3c0, muc0, phic0, tol, dv, dv0, &
	! 	rho0
	! REAL*8	::	slxnu(0:nnuc-1), srxnu(0:nnuc-1), s1xnu(0:nnuc-1), &
	!		s2xnu(0:nnuc-1), s3xnu(0:nnuc-1), xnu0(0:nnuc-1),	delta_xnu0(0:nnuc-1)
	REAL*8, PARAMETER	::	pi = 4.d0 * ATAN(1.d0)
	!REAL*8, POINTER	::	yywt(:, :)
	!REAL*8, POINTER	::	r3c(:), r3l(:), r3r(:)
	!REAL*8, POINTER	::	muc(:), mul(:), mur(:)
	!REAL*8, POINTER	::	phic(:), phil(:), phir(:)
	!REAL*8, POINTER	::	rho(:,:,:), xnu(:,:,:,:)
	!REAL*8, POINTER	::	rhonew(:,:,:), xnunew(:,:,:,:)


	!yywt => yywt_(2:2,2:2)
	!r3c => r3c_(2:2)
	!r3l => r3l_(2:2)
	!r3r => r3r_(2:2)
	!muc => muc_(2:2)
	!mul => mul_(2:2)
	!mur => mur_(2:2)
	!phic => phic_(2:2)
	!phil => phil_(2:2)
	!phir => phir_(2:2)
	!rho => rho_(2:2,2:2,2:2)
	!xnu => xnu_(2:2,2:2,2:2,2:2)
	!rhonew => rhonew_(2:2,2:2,2:2)
	!xnunew => xnunew_(2:2,2:2,2:2,2:2)


 	flipyy = 0

 	DO kk = 0,nph-1
	print*,'kk = ',kk

 		km1 = max(kk-1,0)
 		kp1 = min(kk+1,nph-1)

 		theta_loop: DO jj = 0,nth-1
			print*,'jj = ',jj

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

 			! IF (yywt(jj,kk) == 0) CYCLE theta_loop

 			DO ii = 0,nr-1
				print*,'ii = ',ii

 				im1 = max(ii-1,0)
 				ip1 = min(ii+1,nr-1)

 				sr = (rho(ip1, jj, kk) - rho(ii, jj, kk)) /(r3c(ip1) - r3c(ii) + 1.d-80)
 				sl = (rho(ii, jj, kk) - rho(im1, jj, kk)) /(r3c(ii) - r3c(im1) + 1.d-80)
 				s1rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)
!				print*,rho(ip1,jj,kk), rho(ii,jj,kk),r3c(ip1),r3c(ii)
!				print*,rho(ii,jj,kk)-rho(im1,jj,kk),r3c(ii)-r3c(im1)
!				print*,sr,sl,s1rho

				sr = (rho(ii, jp1, kk) - rho(ii, jj, kk)) /(muc(jp1) - muc(jj) + 1.d-80)
 				sl = (rho(ii, jj, kk) - rho(ii, jm1, kk)) /(muc(jj) - muc(jm1) + 1.d-80)
 				s2rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)
				! print*,sr,sl,s2rho
				print*,'rho(ii,jj,kk)', rho(ii,jj,kk)
				print*,'rho(ii,jm1,kk)', rho(ii, jm1, kk)
				print*,'rho(ii,jj,kk)-rho(ii,jm1,kk)'
				print*,'sr ',sr
				print*,'sl ',sl
				print*,'s2rho ',s2rho

 				sr = (rho(ii, jj, kp1) - rho(ii, jj, kk)) / dphi
 				sl = (rho(ii, jj, kk) - rho(ii, jj, km1)) / dphi
 				s3rho = sign(min(abs(sr), abs(sl)), sr) * heaviside(sl*sr)
!				print*,sr,sl,s3rho

 				srxnu = (xnu(ip1, jj, kk, :) - xnu(ii, jj, kk, :)) / (r3c(ip1) - r3c(ii) + 1.d-80)
 				slxnu = (xnu(ii, jj, kk, :) - xnu(im1, jj, kk, :)) / (r3c(ii) - r3c(im1) + 1.d-80)
 				s1xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)
!				print*,srxnu,slxnu,s1xnu

 				srxnu = (xnu(ii, jp1, kk, :) - xnu(ii, jj, kk, :)) / (muc(jp1) - muc(jj) + 1.d-80)
 				slxnu = (xnu(ii, jj, kk, :) - xnu(ii, jm1, kk, :)) / (muc(jj) - muc(jm1) + 1.d-80)
 				s2xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)
!				print*,srxnu,slxnu,s2xnu

 				srxnu = (xnu(ii, jj, kp1, :) - xnu(ii, jj, kk, :)) / dphi
 				slxnu = (xnu(ii, jj, kk, :) - xnu(ii, jj, km1, :)) / dphi
 				s3xnu = sign(min(abs(srxnu), abs(slxnu)), srxnu) * heaviside(slxnu*srxnu)
!				print*,srxnu,slxnu,s3xnu

 				rho0 = rho(ii, jj, kk)
 				xnu0 = xnu(ii, jj, kk, :)

 				s1xnu = s1xnu / (1. + abs(s1rho * (r3r(ii) - r3l(ii)) / rho0) / 12.)
 				s2xnu = s2xnu / (1. + abs(s2rho * (mur(jj) - mul(jj)) / rho0) / 12.)
 				s3xnu = s3xnu / (1. + abs(s3rho * (phir(kk) - phil(kk)) / rho0) / 12.)
!				print*,s1xnu,s2xnu,s3xnu
 				delta_xnu0 = (s1rho * s1xnu * (r3r(ii) - r3l(ii))**2 &
 											+ s2rho * s2xnu * abs(mur(jj) - mul(jj))**2 &
 											+ s3rho * s3xnu * abs(phir(kk) - phil(kk))**2) / 12.
 				xnu0 = xnu0 - delta_xnu0 / rho0

 				r3c0 = r3c(ii)
 				muc0 = muc(jj)
 				phic0 = phic(kk)

 				nsub = 1
 				nsubmax = 1

				print*,'s1rho ',s1rho
				print*,'s2rho ',s2rho
				print*,'s3rho ',s3rho
				print*,'Calling map_cube'
 				call map_cube(r3l(ii), r3r(ii), mul(jj), mur(jj), phil(kk), phir(kk), tol)
				print*,'Called map_cube'

 			END DO
 		END DO theta_loop
 	END DO

 CONTAINS
 	ELEMENTAL REAL*8 FUNCTION heaviside(x1)
 		IMPLICIT NONE

 		REAL*8, INTENT(in)	::	x1

 		IF (x1 < 0) THEN
 			heaviside = 0
 		ELSE
 			heaviside = 1
 		END IF

 		RETURN

 	END FUNCTION heaviside

 	FUNCTION to_cart(R3, Mu, Phi, flipyy) result(triple)
 		IMPLICIT NONE

 		REAL*8, INTENT(in)	::	R3, Mu, Phi
 		INTEGER, INTENT(in)	::	flipyy
 		REAL*8, DIMENSION(3)	::	triple
 		REAL*8	::	x, y, z, tempx, tempy, tempz, R, Theta

 		R = (3*R3)**(1./3.)
		IF (ny == 1) THEN
			Mu = 0
		END IF
 		Theta = acos(Mu)
		IF (nz == 1) THEN
			Phi = 0
		END IF
 		tempx = R * sin(Theta) * cos(Phi)
 		tempy = R * sin(Theta) * sin(Phi)
 		tempz = R * cos(Theta)
 		IF (flipyy == 0) THEN
 			x = tempx
 			y = tempy
 			z = tempz
 		ELSE
 			x = -tempx
 			y = tempz
 			z = tempy
 		END IF

 		triple(:) = (/ x, y, z /)

 	END FUNCTION to_cart

 	! FUNCTION target_indices(vec, x0, y0, z0, dx, dy, dz, m, n, o)
 	FUNCTION target_indices(vec) RESULT(triple)
 		IMPLICIT NONE

 		REAL*8, DIMENSION(3), INTENT(in)	::	vec
 		REAL*8, DIMENSION(3)	::	triple
 		REAL*8	::	x, y, z

 		x = vec(1)
 		y = vec(2)
 		z = vec(3)

 		triple(1) = int(floor((x-x0)/dx))
 		triple(2) = int(floor((y-y0)/dy))
 		triple(3) = int(floor((z-z0)/dz))

 	END FUNCTION target_indices

 	! RECURSIVE FUNCTION map_cube(Xl, Xr, Yl, Yr, Zl, Zr, tol)
 	RECURSIVE SUBROUTINE map_cube(Xl, Xr, Yl, Yr, Zl, Zr, tol)
 		IMPLICIT NONE

 		REAL*8, DIMENSION(3)	::	c1, c2, c3, c4, c5, c6, c7, c8, res
 		REAL*8	::	Xl, Xr, Yl, Yr, Zl, Zr, tol
 		REAL*8	::	m0, m1, n0, n1, o0, o1, Xc, Yc, Zc, mc, nc, oc

 		c1 = target_indices(to_cart(Xl, Yl, Zl, flipyy))
 		c2 = target_indices(to_cart(Xl, Yl, Zr, flipyy))
 		c3 = target_indices(to_cart(Xl, Yr, Zl, flipyy))
 		c4 = target_indices(to_cart(Xl, Yr, Zr, flipyy))
 		c5 = target_indices(to_cart(Xr, Yl, Zl, flipyy))
 		c6 = target_indices(to_cart(Xr, Yl, Zr, flipyy))
 		c7 = target_indices(to_cart(Xr, Yr, Zl, flipyy))
 		c8 = target_indices(to_cart(Xr, Yr, Zr, flipyy))

!		print*,c1,c2,c3,c4,c5,c6,c7,c8

 		m0 = min(c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1))
 		m1 = max(c1(1), c2(1), c3(1), c4(1), c5(1), c6(1), c7(1), c8(1))
 		n0 = min(c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2))
 		n1 = max(c1(2), c2(2), c3(2), c4(2), c5(2), c6(2), c7(2), c8(2))
 		o0 = min(c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3))
 		o1 = max(c1(3), c2(3), c3(3), c4(3), c5(3), c6(3), c7(3), c8(3))

!		print*,'mno',m0,m1,n0,n1,o0,o1

 		Xc = 0.5*(Xr+Xl)
 		Yc = 0.5*(Yr+Yl)
 		Zc = 0.5*(Zr+Zl)
 		res = target_indices(to_cart(Xc, Yc, Zc, flipyy))
 		mc = res(1)
 		nc = res(2)
 		oc = res(3)

		print*,Xc,Yc,Zc,mc,nc,oc

 		IF (nph==1) THEN
 			dv = abs(Xr-Xl) * abs(Yr-Yl) * dphi
 		ELSE
 			dv = abs(Xr-Xl) * abs(Yr-Yl) * abs(Zr-Zl) * yywt(jj,kk)
 		END IF
!		print*,dv
 		IF (ny==1) THEN
			IF (nz==1) THEN
				dv0 = (4./3.)*pi*(3*mc**2+3*mc+1)*dx**3
			ELSE
	 			dv0 = pi*(2*mc+1)*dx**2*dz
			END IF
 		ELSE
 			! dv0 = float(dx*dy*dz)
 			dv0 = dx*dy*dz
 		END IF
		print*,dv0
 		IF ((max(m0,m1)>=nx) .OR. (max(n0,n1)>=ny) .OR. (max(o0,o1)>=nz)) THEN
 			RETURN
 		END IF
 		IF ((min(m0,m1)<0) .OR. (min(n0,n1)<0) .OR. (min(o0,o1)<0)) THEN
 			RETURN
 		END IF

 		IF ( ((m0==m1) .AND. (n0==n1) .AND. (o0==o1)) .OR. (dv < tol*dv0) ) THEN
			! print*,rho0
			! print*,Xc-r3c0
			! print*,Yc-muc0
			! print*,Zc-phic0
			print*,'s1rho ',s1rho
			print*,'s2rho ',s2rho
			print*,'s3rho ',s3rho
 			rhonew(mc,nc,oc) = rhonew(mc,nc,oc) + (rho0 + s1rho*(Xc-r3c0) &
 				+ s2rho*(Yc-muc0) + s3rho*(Zc-phic0)) * dv/dv0
			print*,'rhonew(mc,nc,oc) ',rhonew(mc,nc,oc)
 			xnunew(mc,nc,oc,:) = xnunew(mc,nc,oc,:) + (rho0 + s1rho*(Xc-r3c0) &
 				+ s2rho*(Yc-muc0) + s3rho*(Zc-phic0))*(xnu0 + s1xnu*(Xc-r3c0) &
 				+ s2xnu*(Yc-muc0) + s3xnu*(Zc-phic0)) * dv/dv0
!			print*,xnunew
 			RETURN

!		print*,'rhonew = ',rhonew
!		print*,'xnunew = ',xnunew

 		ELSE
 			nsub = nsub + 1
 			nsubmax = max(nsub,nsubmax)
 			call map_cube(Xl, Xc, Yl, Yc, Zl, Zc, tol)
 			call map_cube(Xl, Xc, Yc, Yr, Zl, Zc, tol)
 			call map_cube(Xc, Xr, Yl, Yc, Zl, Zc, tol)
 			call map_cube(Xc, Xr, Yc, Yr, Zl, Zc, tol)
 			IF (nph > 1) THEN
 				call map_cube(Xl, Xc, Yl, Yc, Zc, Zr, tol)
 				call map_cube(Xl, Xc, Yc, Yr, Zc, Zr, tol)
 				call map_cube(Xc, Xr, Yl, Yc, Zc, Zr, tol)
 				call map_cube(Xc, Xr, Yc, Yr, Zc, Zr, tol)
 			END IF
 		END IF

 			nsub = nsub - 1

 	END SUBROUTINE map_cube
END SUBROUTINE mapping
