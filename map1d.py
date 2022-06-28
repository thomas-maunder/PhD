#25.18

import numpy as np
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

# Use this function to map 2D->2D or 3D->3D
def to_cart (R3, Mu, Phi, ny, nz):
    R = np.cbrt(3*R3)
    if ny == 1:
        Mu = 0
    Theta = np.arccos(Mu)
    if nz == 1:
        Phi = 0
    x = R * np.sin (Theta) * np.cos (Phi)
    y = R * np.sin (Theta) * np.sin (Phi)
    z = R * np.cos (Theta)
#    print(R, Theta, Phi,x,y,z)
    return x,y,z

#Use this_function to map 3D->2D
def to_cart_32D (R, Theta, Phi):
    n1, n2, n3 = (1,0,0)
    xi = R * np.sin (Theta) * np.cos (Phi)
    yi = R * np.sin (Theta) * np.sin (Phi)
    zi = R * np.cos (Theta)

    x = np.linalg.norm(np.cross((xi,yi,zi),(n1,n2,n3)))
    y = 0.0
    z = np.dot(a,b)

    return x,y,z

def map2d(model, rl, rr, thl, thr, phil, phir, x0, x1, z0, z1, nx, ny, nz, tol):

    rho = model.den()
    xnu = model.xnu()

    nr  = np.size(rl)
    nth = np.size(thl)
    nph = np.size(phil)
    nnuc = np.size(xnu[0,0,0,:])

    xigmap = 0.0

# calculate explosion properties
    dv_r = 1./3. * (model.xzr()**3-model.xzl()**3)
    dv_theta = abs(np.cos(model.yzl())-np.cos(model.yzr()))
    dv_phi = model.zzr() - model.zzl()
    dv = dv_r[:,np.newaxis,np.newaxis]*dv_theta[np.newaxis,:,np.newaxis]*dv_phi[np.newaxis,np.newaxis,:]
    xig=np.sum(xnu[:,:,:,14:nnuc-1],axis=3)
    print('Ejecta mass',np.sum(dv*rho*(1.0-xnu[:,:,:,2])/1.99e33),'M_sun')
    print('Explosion energy',np.sum(dv*rho*model.ene()),'erg')
    print('Nickel/IG mass',np.sum(dv*rho*xig/1.99e33),'M_sum')
    print('He mass',np.sum(dv*rho*xnu[:,:,:,4]/1.99e33),'M_sun')
    print('C mass',np.sum(dv*rho*xnu[:,:,:,5]/1.99e33),'M_sun')
    print('O mass',np.sum(dv*rho*xnu[:,:,:,7]/1.99e33),'M_sun')
    print('Ne mass',np.sum(dv*rho*xnu[:,:,:,8]/1.99e33),'M_sun')
    print('Mg mass',np.sum(dv*rho*xnu[:,:,:,9]/1.99e33),'M_sun')
    print('Si mass',np.sum(dv*rho*xnu[:,:,:,10]/1.99e33),'M_sun')
    print('S mass',np.sum(dv*rho*xnu[:,:,:,11]/1.99e33),'M_sun')
    print('Ar mass',np.sum(dv*rho*xnu[:,:,:,12]/1.99e33),'M_sun')
    print('Ca mass',np.sum(dv*rho*xnu[:,:,:,13]/1.99e33),'M_sun')
    print('Ti mass',np.sum(dv*rho*xnu[:,:,:,14]/1.99e33),'M_sun')

    dm = np.sum (rho[:,:,:]*dv,axis=(1,2))
    xhe = np.sum(rho[:,:,:]*xnu[:,:,:,4]*dv,axis=(1,2)) / dm
    xc  = np.sum(rho[:,:,:]*xnu[:,:,:,5]*dv,axis=(1,2)) / dm
    xo  = np.sum(rho[:,:,:]*xnu[:,:,:,7]*dv,axis=(1,2)) / dm
    xne = np.sum(rho[:,:,:]*xnu[:,:,:,8]*dv,axis=(1,2)) / dm
    xmg = np.sum(rho[:,:,:]*xnu[:,:,:,9]*dv,axis=(1,2)) / dm
    xsi = np.sum(rho[:,:,:]*xnu[:,:,:,10]*dv,axis=(1,2)) / dm
    xs  = np.sum(rho[:,:,:]*xnu[:,:,:,11]*dv,axis=(1,2)) / dm
    xar = np.sum(rho[:,:,:]*xnu[:,:,:,12]*dv,axis=(1,2)) / dm
    xca = np.sum(rho[:,:,:]*xnu[:,:,:,13]*dv,axis=(1,2)) / dm
    xti = np.sum(rho[:,:,:]*xnu[:,:,:,14]*dv,axis=(1,2)) / dm
    xigav = np.sum(rho[:,:,:]*xig          *dv,axis=(1,2)) / dm

    vmid = np.sum(rho[:,:,:]*model.vex()      *dv,axis=(1,2)) / dm

    plt.close('all')
    plt.clf()

    f, ax = plt.subplots(1, sharex=True, figsize=(6,7))
    f.subplots_adjust(hspace=0.0,left=0.18,bottom=0.1,top=0.96,right=0.9)

    qc1,=ax.plot(model.xzn()/1e5,xhe,label='He')
    qc2,=ax.plot(model.xzn()/1e5,xc ,label='C')
    qc3,=ax.plot(model.xzn()/1e5,xo ,label='O')
    qc4,=ax.plot(model.xzn()/1e5,xne,label='Ne')
    qc5,=ax.plot(model.xzn()/1e5,xmg,label='Mg')
    qc6,=ax.plot(model.xzn()/1e5,xsi,label='Si')
    qc7,=ax.plot(model.xzn()/1e5,xca,label='Ca')
    qc8,=ax.plot(model.xzn()/1e5,xigav,label='Fe group')

    ax.legend(handles=[qc1,qc2,qc3,qc4,qc5,qc6,qc7,qc8],fontsize=11,handlelength=1.4,columnspacing=0.4,loc='upper left'\
,ncol=4)
#    ax.set_xlabel(r'$v_r\ [\mathrm{km}\, \mathrm{s}^{-1}]$')
    ax.set_xlabel(r'$r\ [\mathrm{km}]$')
    ax.set_ylabel(r'$\Delta M_i\ [M_\odot]$')

    plt.savefig('composition.pdf')
    plt.close('all')

#    for i in range(nnuc):
#        xnu[:,:,:,i]=xnu[:,:,:,i]*rho[:,:,:]

    mul = np.cos(thl)
    mur = np.cos(thr)

    r3l = rl ** 3 / 3.0
    r3r = rr ** 3 / 3.0

    #cell centre coordinates in PROMETHEUS
    r3c  = 0.5 * (r3l+r3r)
    muc  = 0.5 * (mul+mur)
    thc  = 0.5 * (thl+thr)
    phic = 0.5 * (phil+phir)
    dphi = 2.0*np.pi/nph

    if (nph == 1):
        phil[:] = -0.0
        phir[:] =  0.0
        phic[:] =  0.0

    #cell spcaing in ARTIS
    dx = (x1-x0) / float(nx)
    dy = (x1-x0) / float(ny) # same as dx if we're in 3D, otherwise this should not be used anyway
    dz = (z1-z0) / float(nz)

    y0 = x0

    rhonew = np.zeros([nx,ny,nz])
    xnunew = np.zeros([nx,ny,nz,nnuc])
    xnunew [:,:,:,4]= 1e-50

    # Map each cell
    for kk in tqdm(range (nph), desc='Phi loop'):

        km1 = max(kk-1,0)
        kp1 = min(kk+1,nph-1)

        for jj in tqdm(range (nth), desc='Theta loop', leave=False):

            jm1 = max(jj-1,0)
            jp1 = min(jj+1,nth-1)

            for ii in tqdm(range (nr), desc='Radial loop', leave=False):

                im1 = max(ii-1,0)
                ip1 = min(ii+1,nr-1)

# minmod slopes in each direction
                sr = (rho[ip1,jj ,kk ]-rho[ii , jj,kk ]) / (r3c[ip1] - r3c [ii] + 1e-80)
                sl = (rho[ii ,jj ,kk ]-rho[im1, jj,kk ]) / (r3c[ii ] - r3c [im1] + 1e-80)
                s1rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (rho[ii ,jp1,kk ]-rho[ii ,jj ,kk]) / (muc[jp1] - muc [jj] + 1e-80)
                sl = (rho[ii ,jj ,kk ]-rho[ii ,jm1,kk]) / (muc[jj ] - muc [jm1] + 1e-80)
                s2rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (rho[ii ,jj ,kp1]-rho[ii , jj,kk ]) / dphi
                sl = (rho[ii ,jj ,kk ]-rho[ii , jj,km1]) / dphi
                s3rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ip1,jj ,kk ,:]-xnu[ii , jj,kk ,:]) / (r3c[ip1] - r3c [ii] + 1e-80)
                sl = (xnu[ii ,jj ,kk ,:]-xnu[im1, jj,kk ,:]) / (r3c[ii ] - r3c [im1] + 1e-80)
                s1xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ii ,jp1,kk ,:]-xnu[ii ,jj ,kk,:]) / (muc[jp1] - muc [jj]  + 1e-80)
                sl = (xnu[ii ,jj ,kk ,:]-xnu[ii ,jm1,kk,:]) / (muc[jj ] - muc [jm1] + 1e-80)
                s2xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ii ,jj ,kp1,:]-xnu[ii , jj,kk ,:]) / dphi
                sl = (xnu[ii ,jj ,kk ,:]-xnu[ii , jj,km1,:]) / dphi
                s3xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                rho0 = rho[ii,jj,kk]
                xnu0 = xnu[ii,jj,kk,:]

# additional slope limiting for mass fraction
                s1xnu = s1xnu / (1.0 + np.abs(s1rho * (r3r[ii] - r3l[ii]) / rho0) / 12.0)
                s2xnu = s2xnu / (1.0 + np.abs(s2rho * (mur[jj] - mul[jj]) / rho0) / 12.0)
                s3xnu = 0*s3xnu / (1.0 + np.abs(s3rho * (phir[kk] - phil[kk]) / rho0) / 12.0)

                delta_xnu0 = (s1rho * s1xnu * (r3r[ii] - r3l[ii]) ** 2 \
                    + s2rho * s2xnu * np.abs(mur[jj] - mul[jj]) ** 2 \
                    + s3rho * s3xnu * np.abs(phir[kk] - phil[kk]) ** 2) / 12.0
                xnu0 = xnu0 - delta_xnu0 / rho0
#                sys.exit()

                r3c0  = r3c[ii]
                muc0  = muc[jj]
                phic0 = phic[kk]

#debug
#                s1rho = 0.0
#                s2rho = 0.0
#                s3rho = 0.0

                def target_indices (vec):
                    x, y, z = vec
                    nonlocal x0, y0, z0, dx, dy, dz
                    m = int(np.floor ((x-x0) / dx))
                    n = int(np.floor ((y-x0) / dy))
                    o = int(np.floor ((z-z0) / dz))
                    if (nz==1):
                        n = 0
                        o = 0
                    return m, n, o

                def map_cube(Xl, Xr, Yl, Yr, Zl, Zr, tol):

                    nonlocal rho,rhonew,xnu,xnunew,rho0,xnu0,s1rho,s2rho,s3rho,s1xnu,s2xnu,s3xnu,nx,ny,nz,nsub,nsubmax,dx,dy,dz,r3c0,muc0,phic0,xigmap

#                    print('map_cube',Xl,Xr,Yl,Yr,Zl,Zr)
# maybe doing all of the eight corners is overkill...
                    c1 = target_indices (to_cart(Xl, Yl, Zl, ny, nz))
                    c2 = target_indices (to_cart(Xl, Yl, Zr, ny, nz))
                    c3 = target_indices (to_cart(Xl, Yr, Zl, ny, nz))
                    c4 = target_indices (to_cart(Xl, Yr, Zr, ny, nz))
                    c5 = target_indices (to_cart(Xr, Yl, Zl, ny, nz))
                    c6 = target_indices (to_cart(Xr, Yl, Zr, ny, nz))
                    c7 = target_indices (to_cart(Xr, Yr, Zl, ny, nz))
                    c8 = target_indices (to_cart(Xr, Yr, Zr, ny, nz))

# Check whether we're overlapping with more than one target cell
                    m0 = min(c1[0], c2[0], c3[0], c4[0], c5[0], c6[0], c7[0], c8[0])
                    m1 = max(c1[0], c2[0], c3[0], c4[0], c5[0], c6[0], c7[0], c8[0])
                    n0 = min(c1[1], c2[1], c3[1], c4[1], c5[1], c6[1], c7[1], c8[1])
                    n1 = max(c1[1], c2[1], c3[1], c4[1], c5[1], c6[1], c7[1], c8[1])
                    o0 = min(c1[2], c2[2], c3[2], c4[2], c5[2], c6[2], c7[2], c8[2])
                    o1 = max(c1[2], c2[2], c3[2], c4[2], c5[2], c6[2], c7[2], c8[2])

                    Xc = 0.5 * (Xr+Xl)
                    Yc = 0.5 * (Yr+Yl)
                    Zc = 0.5 * (Zr+Zl)
                    mc, nc, oc = target_indices (to_cart(Xc, Yc, Zc, ny, nz))

                    if (nph==1): dv = abs(Xr-Xl) * abs(Yr-Yl) * dphi
                    if (nph>1): dv = abs(Xr-Xl) * abs(Yr-Yl) * abs(Zr-Zl)
                    if (ny==1):
                        if (nz==1):
                            dv0 = (4/3)*np.pi*(3*mc**2+3*mc+1)*dx**3
                        else:
                            dv0 = np.pi*(2*mc+1)*dx**2*dz
                    if (ny>1): dv0 = float(dx*dy*dz)

#                    if (nsub==1):
#                        print ('map_cube2',m0,m1,n0,n1,o0,o1,mc,nc,oc,dv,tol*dv0,nsub)
#                        print (np.cbrt(Xc*3),Yc,Zc,abs(Xr-Xl),abs(Yr-Yl),nph)
#                    print ('map_cube2',m0,m1,n0,n1,o0,o1,Xr-Xl,Yr-Yl,Zr-Zl)
                    if ((max(m0,m1)>=nx) | (max(n0,n1)>=ny) | (max(o0,o1)>=nz)): return
                    if ((min(m0,m1)<0) | (min(n0,n1)<0) | (min(o0,o1)<0)): return

                    if ((m0==m1) & (n0==n1) & (o0==o1)) | (dv < tol*dv0):
                        #map entire cube to one destination cell
#                        print ('converged',mc,nc,oc,(Xc*3)**0.33333,Yc,Zc)
#                        if ((rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0))<0):
#                            print('rho<0',rho0,s1rho,s2rho,s3rho,Xc,Yc,Zc)
#                            sys.exit()
                        rhonew[mc,nc,oc]   += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * dv/dv0
# 1st order
#                        xnunew[mc,nc,oc,:] += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * xnu0 * dv/dv0
# 2nd order
                        xnunew[mc,nc,oc,:] += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * \
                            (xnu0 + s1xnu * (Xc-r3c0) + s2xnu * (Yc-muc0) + s3xnu * (Zc-phic0)) * \
                            dv/dv0
#                        xigmap +=  (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * np.sum(xnu0[14:nnuc-1]) * dv
                        return
                    else:
                        #divide into eight subcubes and map them
                        nsub = nsub + 1
                        nsubmax = max(nsub,nsubmax)
                        map_cube(Xl, Xc, Yl, Yc, Zl, Zc, tol)
                        map_cube(Xl, Xc, Yc, Yr, Zl, Zc, tol)
                        map_cube(Xc, Xr, Yl, Yc, Zl, Zc, tol)
                        map_cube(Xc, Xr, Yc, Yr, Zl, Zc, tol)
                        # print(nph)
                        if (nph > 1):
                            map_cube(Xl, Xc, Yl, Yc, Zc, Zr, tol)
                            map_cube(Xl, Xc, Yc, Yr, Zc, Zr, tol)
                            map_cube(Xc, Xr, Yl, Yc, Zc, Zr, tol)
                            map_cube(Xc, Xr, Yc, Yr, Zc, Zr, tol)

                        nsub = nsub - 1

                nsub  = 1
                nsubmax = 1
                map_cube(r3l[ii], r3r[ii], mul[jj], mur[jj], phil[kk], phir[kk], tol)

                # print ('Mapped zone',ii,jj,kk,nsubmax)

    rhonew = np.maximum(rhonew, 1e-50)
    for i in range(nnuc):
        xnunew[:,:,:,i]=xnunew[:,:,:,i]/rhonew[:,:,:]

    x = np.arange(x0,x1,dx)+0.5*dx
    z = np.arange(z0,z1,dz)+0.5*dz
    r = np.sqrt(np.outer(x**2,z**0) + np.outer(x**0,z**2))

    xignew=np.sum(xnunew[:,:,:,14:nnuc-1],axis=3)

    plt.plot(x,np.log10(rhonew[:,0,:]))
    plt.plot(model.xzn(), np.log10(model.den()[:,0,0]))
    plt.plot(model.xzn(), np.log10(model.den()[:,112,0]))
    plt.plot(model.xzn(), np.log10(model.den()[:,48,0]))
    plt.plot(model.xzn(), np.log10(model.den()[:,32,0]))
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$\rho\ [\mathrm{g/cm}]$')
    plt.savefig('density.pdf')

    plt.plot(x,np.log10(rhonew[:,0,:]), label='1D Map')
    plt.plot(model.xzn(), np.log10(model.den()[:,0,0]), label='2D')
    plt.plot(model.xzn(), np.log10(model.den()[:,112,0]))
    plt.plot(model.xzn(), np.log10(model.den()[:,48,0]))
    plt.plot(model.xzn(), np.log10(model.den()[:,32,0]))
    plt.xlim(min(model.xzn()), max(model.xzn()))
    plt.ylim(min(np.log10(model.den()[:,0,0])), max(np.log10(model.den()[:,0,0])))
    plt.legend()
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$\rho\ [\mathrm{g/cm}]$')
    plt.savefig('density_zoomed.pdf')

    plt.plot(x,xignew[:,0,:])
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$\Delta M_{\mathrm{Iron group}}\ [M_\odot]$')
    plt.savefig('iron_group.pdf')

    plt.plot(x,xnunew[:,0,:,7])
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$\Delta M_O\ [M_\odot]$')
    plt.savefig('oxygen.pdf')

    plt.plot(x,xnunew[:,0,:,4])
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$\Delta M_{He}\ [M_\odot]$')
    plt.savefig('helium.pdf')

# Density and Detailed composition for ARTIS
# Assuming slightly proton-rich ejecta

    xartis = np.zeros([nx,ny,nz,1+30])

    xartis[:,:,:,0] = rhonew [:,:,:]

    xartis[:,:,:,1]  = xnunew [:,:,:,0] + xnunew [:,:,:,1] # n, p -> H
    xartis[:,:,:,2]  = xnunew [:,:,:,4] # He
    xartis[:,:,:,3]  = 0.0 # Li
    xartis[:,:,:,4]  = 0.0 # Be
    xartis[:,:,:,5]  = 0.0 # B
    xartis[:,:,:,6]  = xnunew [:,:,:,5] # C
    xartis[:,:,:,7]  = xnunew [:,:,:,6] # N
    xartis[:,:,:,8]  = xnunew [:,:,:,7] # O
    xartis[:,:,:,10] = xnunew [:,:,:,8] # Ne
    xartis[:,:,:,12] = xnunew [:,:,:,9] # Mg
    xartis[:,:,:,14] = xnunew [:,:,:,10] # Si
    xartis[:,:,:,16] = xnunew [:,:,:,11] # S
    xartis[:,:,:,18] = xnunew [:,:,:,12] # Ar
    xartis[:,:,:,20] = xnunew [:,:,:,13] + xignew * 3e-3 # Ca
    xartis[:,:,:,22] = xignew * 1.5e-3 # Ti
    xartis[:,:,:,24] = xignew * 2.5e-3 # Cr
    xartis[:,:,:,26] = xignew * 2.0e-3 # Fe52

    xdec = 0.5**(model.time()/(6.077*86400.0))
    xcr48 = xignew * 2.5e-3 # Cr
    xfe52 = xignew * 2.0e-3 # Fe52
    xco56 = xignew * 0.99392 * (1.0-xdec)
    xni56 = xignew * 0.99392 * xdec

    xartis[:,:,:,27] = xco56 # Co56
    xartis[:,:,:,28] = xni56 + xignew *  8e-5 # Ni56, Ni58

# Normalise mass fractions
    xsum = np.sum(xartis[:,:,:,1:],3)
    for i in range (30):
        xartis[:,:,:,i+1] = xartis[:,:,:,i+1] / xsum
    xcr48 = xcr48 / xsum
    xfe52 = xfe52 / xsum
    xco56 = xco56 / xsum
    xni56 = xni56 / xsum
    xignew = xignew / xsum

    print('Writing input files for ARTIS.')

    f = open('model.txt', 'tw')
    g = open('abundances.txt', 'tw')


    ij = 0
    f.write ("%d %d \n" % (nx, nz))
    time = model.time()
    f.write ("%16.8e \n" % (time/86400.))
    f.write ("%16.8e \n" % (max(x1/time,z1/time)))
    for jj in range(nz):
        for ii in range(nx):
            ij += 1
            f.write ("%d %16.8e %16.8e %16.8e \n %16.8e %16.8e %16.8e %16.8e %16.8e \n" % (ij, x[ii], z[jj], xartis[ii,0,jj,0], xignew[ii,0,jj], xni56[ii,0,jj], xco56[ii,0,jj], xfe52[ii,0,jj], xcr48[ii,0,jj]))
            out = " ".join(("%16.8e " %dat) for dat in tuple(xartis[ii,0,jj,1:]))
            g.write(("%d " % ij) + out + " \n")

    f.closed
    g.closed


# Check mass of mapped model
    dv0 = 2*np.pi*np.outer(x,z**0) *dx*dz
    print("Mass of mapped model:",np.sum(xartis[:,0,:,0]*dv0)/1.99e33,"M_sun")
    print("Iron group:          ",np.sum(xartis[:,0,:,0]*xignew[:,0,:]*dv0)/1.99e33,"M_sun")
#    print("Iron group:          ",np.sum(xartis[:,0,:,0]*xignew[:,0,:]*dv0)/1.99e33,xigmap/1.99e33,"M_sun")


    print('Done!')
    return rhonew,xnunew,x,z,r
