SUBROUTINE mapping_(nr, nth, nph, yywt, nx, ny, nz, nnuc, x0, y0, &
 	z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew, &
 	muc, mul, mur, phic, phil, phir, tol)
	IMPLICIT NONE

	INTEGER, INTENT(in)	::	nr, nth, nph, nx, ny, nz, nnuc
	REAL*8, INTENT(in)	::	x0, y0, z0
	REAL*8, INTENT(in)	::	yywt(:,:)
	REAL*8, INTENT(in)	::	r3c(:), r3l(:), r3r(:)
	REAL*8, INTENT(in)	::	muc(:), mul(:), mur(:)
	REAL*8, INTENT(in)	::	phic(:), phil(:),	phir(:)
	REAL*8, INTENT(in)	::	rho(:,:,:)
	REAL*8, INTENT(in)	::	xnu(:,:,:,:)
	REAL*8, INTENT(inout)	::	rhonew(:,:,:)
	REAL*8, INTENT(inout)	::	xnunew(:,:,:,:)
	REAL*8, INTENT(in)	::	dphi, dx, dy, dz, tol

 call  mapping(nr, nth, nph, yywt, nx, ny, nz, nnuc, x0, y0, &
 	z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew, &
 	muc, mul, mur, phic, phil, phir, tol)

end SUBROUTINE mapping_
