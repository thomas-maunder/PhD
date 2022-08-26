#! /bin/env python3

import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from . import mapping


# Use this function to map 2D->2D or 3D->3D
def to_cart (R3, Mu, Phi, flipyy = 0):
    R = np.cbrt(3*R3)
    if ny == 1:
        Mu = 0
    Theta = np.arccos(Mu)
    if nz == 1:
        Phi = 0
    x = R * np.sin (Theta) * np.cos (Phi)
    y = R * np.sin (Theta) * np.sin (Phi)
    z = R * np.cos (Theta)
    if flipyy == 0:
#        print('flipyy',0)
        return  x,y,z
    else:
#        print('flipyy',1)
        return -x,z,y

#Use this_function to map 3D->2D
def to_cart_32D (R, Theta, Phi):
    n1, n2, n3 = (1,0,0)
    xi = R * np.sin (Theta) * np.cos (Phi)
    yi = R * np.sin (Theta) * np.sin (Phi)
    zi = R * np.cos (Theta)

    x = np.linalg.norm(np.cross((xi,yi,zi),(n1,n2,n3)))
    y = 0.0
    z = np.dot(a,b)

    return x,y,znth, nph, yywt, nx, ny, nz, nnuc, x0, y0, z0
    # return x,y,znth, nph, yywt, nx, ny, nz, nnuc, x0, y0, &
	# z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew, &
	# muc, mul, mur, phic, phil, phir)

# @jit
def map3d(model, rl, rr, thl, thr, phil, phir, x0, x1, z0, z1, nx, ny, nz, tol):

    print('Starting Mapping...')

    rho = model.den()
    xnu = model.xnu()

    # print(np.shape(rho))
    # exit()

    thl = np.concatenate((thl, thl))
    thr = np.concatenate((thr, thr))

    nr  = np.size(rl)
    nth = np.size(thl)
    nph = np.size(phil)
    nnuc = np.size(xnu[0,0,0,:])

    yzl = np.concatenate((model.yzl(), model.yzl()))
    yzr = np.concatenate((model.yzr(), model.yzr()))
    yzn = np.concatenate((model.yzn(), model.yzn()))

    yywt = np.concatenate((model.yinyang_weight(), model.yinyang_weight()))

    xigmap = 0.0

    solmassi = 1. / 1.99e33

# calculate explosion properties
    dv_r = 1./3. * (model.xzr()**3-model.xzl()**3)
    dv_theta = abs(np.cos(yzl)-np.cos(yzr))
    dv_phi = model.zzr() - model.zzl()
    # CHECK YINYANG GRID AND WEIGHT
    dv = dv_r[:,np.newaxis,np.newaxis]*dv_theta[np.newaxis,:,np.newaxis]*dv_phi[np.newaxis,np.newaxis,:]*yywt[np.newaxis,:,:]
    # print(np.shape(dv))
    # exit()
    xig=np.sum(xnu[:,:,:,14:nnuc-1],axis=3)
    initial_mass = np.sum(dv*rho*(1.0-xnu[:,:,:,2])*solmassi)
    initial_IG_mass = np.sum(dv*rho*xig*solmassi)
    print('Ejecta mass',np.sum(dv*rho*(1.0-xnu[:,:,:,2])*solmassi),'M_sun')
    print('Explosion energy',np.sum(dv*rho*model.ene()),'erg')
    print('Nickel/IG mass',np.sum(dv*rho*xig*solmassi),'M_sum')
    print('He mass',np.sum(dv*rho*xnu[:,:,:,4]*solmassi),'M_sun')
    print('C mass',np.sum(dv*rho*xnu[:,:,:,5]*solmassi),'M_sun')
    print('O mass',np.sum(dv*rho*xnu[:,:,:,7]*solmassi),'M_sun')
    print('Ne mass',np.sum(dv*rho*xnu[:,:,:,8]*solmassi),'M_sun')
    print('Mg mass',np.sum(dv*rho*xnu[:,:,:,9]*solmassi),'M_sun')
    print('Si mass',np.sum(dv*rho*xnu[:,:,:,10]*solmassi),'M_sun')
    print('S mass',np.sum(dv*rho*xnu[:,:,:,11]*solmassi),'M_sun')
    print('Ar mass',np.sum(dv*rho*xnu[:,:,:,12]*solmassi),'M_sun')
    print('Ca mass',np.sum(dv*rho*xnu[:,:,:,13]*solmassi),'M_sun')
    print('Ti mass',np.sum(dv*rho*xnu[:,:,:,14]*solmassi),'M_sun')

    dmi = 1. / np.sum (rho[:,:,:]*dv,axis=(1,2))
    xhe = np.sum(rho[:,:,:]*xnu[:,:,:,4]*dv,axis=(1,2)) * dmi
    xc  = np.sum(rho[:,:,:]*xnu[:,:,:,5]*dv,axis=(1,2)) * dmi
    xo  = np.sum(rho[:,:,:]*xnu[:,:,:,7]*dv,axis=(1,2)) * dmi
    xne = np.sum(rho[:,:,:]*xnu[:,:,:,8]*dv,axis=(1,2)) * dmi
    xmg = np.sum(rho[:,:,:]*xnu[:,:,:,9]*dv,axis=(1,2)) * dmi
    xsi = np.sum(rho[:,:,:]*xnu[:,:,:,10]*dv,axis=(1,2)) * dmi
    xs  = np.sum(rho[:,:,:]*xnu[:,:,:,11]*dv,axis=(1,2)) * dmi
    xar = np.sum(rho[:,:,:]*xnu[:,:,:,12]*dv,axis=(1,2)) * dmi
    xca = np.sum(rho[:,:,:]*xnu[:,:,:,13]*dv,axis=(1,2)) * dmi
    xti = np.sum(rho[:,:,:]*xnu[:,:,:,14]*dv,axis=(1,2)) * dmi
    xigav = np.sum(rho[:,:,:]*xig          *dv,axis=(1,2)) * dmi

    vmid = np.sum(rho[:,:,:]*model.vex()      *dv,axis=(1,2)) * dmi

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
    dphi = 2.0 * np.pi / nph

    if (nph == 1):
        phil[:] = -0.0
        phir[:] =  0.0
        phic[:] =  0.0

    #cell spcaing in ARTIS
    y0 = x0
    y1 = x1

    dx = (x1-x0) / float(nx)
    if ny > 1:
        dy = (y1-y0) / float(ny) # same as dx if we're in 3D, otherwise this should not be used anyway
        if nz > 1:
            dz = (z1-z0) / float(nz)
        else:
            dz = 1
    else:
        dy = 1
        dz = 1


    rhonew = np.zeros([nx,ny,nz], order='F')
    xnunew = np.zeros([nx,ny,nz,nnuc], order='F')

    flipyy = 0

    # FORTRAN FUNCTION IS CALLED HERE

    print('Starting FORTRAN code...')
    mapping(nr, nth, nph, yywt, nx, ny, nz, nnuc, x0, y0, z0, rho, xnu, r3c, r3l, r3r, dx, dy, dz, dphi, rhonew, xnunew,  muc, mul, mur, phic, phil, phir, tol)
    print('Exited FORTRAN code')

#    rhonew = np.maximum(rhonew, 1e-50)
    rhonew = np.maximum(rhonew, 1e-5 * np.min(rho))

 #   for i in range(nnuc):
    xnunew[:,:,:,:]=xnunew[:,:,:,:]/rhonew[:,:,:,np.newaxis]

    x = np.arange(x0,x1,dx)+0.5*dx
    if ny > 1:
        y = np.arange(y0,y1,dy)+0.5*dy
    if nz > 1:
        z = np.arange(z0,z1,dz)+0.5*dz
    else:
        z = 0

    print(len(x))
    # r = np.sqrt(np.outer(x**2,z**0) + np.outer(x**0,z**2))
    r = x
    # r = np.sqrt(x**2 + y**2 + z**2)
    # GENERALISE R TO 3D

    xignew=np.sum(xnunew[:,:,:,14:nnuc-1],axis=3)

    # print('Saving plots...')
    # plt.contourf(z,x,np.log10(rhonew[:,ny//2,:]),50, cmap='inferno')
    # plt.xlabel(r'$x\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.savefig('density.pdf')
    #
    # plt.contourf(y,x,np.log10(rhonew[:,:,ny//2]),50, cmap='inferno')
    # plt.xlabel(r'$x\ [\mathrm{cm}]$')
    # plt.ylabel(r'$y\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    #
    # plt.contourf(z,y,np.log10(rhonew[ny//2,:,:]),50, cmap='inferno')
    # plt.xlabel(r'$y\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    #
    # plt.pcolor(z,y,np.log10(rhonew[ny//2,:,:]), cmap='inferno')
    # plt.xlabel(r'$y\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    #
    # plt.contourf(z,x,xignew[:,ny//2,:],50)
    # plt.xlabel(r'$x\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.savefig('iron_group.pdf')
    #
    # plt.contourf(z,x,xnunew[:,ny//2,:,7],50)
    # plt.xlabel(r'$x\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.savefig('oxygen.pdf')
    #
    # plt.contourf(z,x,xnunew[:,ny//2,:,4],50)
    # plt.xlabel(r'$x\ [\mathrm{cm}]$')
    # plt.ylabel(r'$z\ [\mathrm{cm}]$')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.savefig('helium.pdf')
    # print('Plots saved')

# Density and Detailed composition for ARTIS
# Assuming slightly proton-rich ejecta

    # print(nx, ny, nz)

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
    xsumi = 1. / np.maximum(1.e-10, np.sum(xartis[:,:,:,1:],3))
    xartis[:,:,:,1:31] = xartis[:,:,:,1:31] * xsumi[:,:,:,np.newaxis]
    xcr48 = xcr48 * xsumi
    xfe52 = xfe52 * xsumi
    xco56 = xco56 * xsumi
    xni56 = xni56 * xsumi
    xignew = xignew * xsumi

#     print('Writing input files for ARTIS.')
#
#     f = open('model.txt', 'tw')
#     g = open('abundances.txt', 'tw')
#
#
# # CHECK WITH STUART AND FINN HOW THEY SET UP 3D INPUT DATA
#     ij = 0
#     f.write ("%d \n" % (nx*ny*nz))
#     time = model.time()
#     f.write ("%16.8e \n" % (time/86400.))
#     f.write ("%16.8e \n" % (max(x1/time,z1/time)))
#     for kk in range(nz):
#         for jj in range(ny):
#             for ii in range(nx):
#                 ij += 1
#                 f.write ("%d %16.8e %16.8e %16.8e %16.8e \n %16.8e %16.8e %16.8e %16.8e %16.8e \n" % (ij, x[ii], y[jj], z[kk], xartis[ii,jj,kk,0], xignew[ii,jj,kk], xni56[ii,jj,kk], xco56[ii,jj,kk], xfe52[ii,jj,kk], xcr48[ii,jj,kk]))
#                 out = " ".join(("%16.8e " %dat) for dat in tuple(xartis[ii,jj,kk,1:]))
#                 g.write(("%d " % ij) + out + " \n")
#
#     f.closed
#     g.closed

    print('Writing input files for ARTIS.')

    f = open('model_1d_100.txt', 'w')
    # h = open('model.txt', 'tw')
    g = open('abundances.txt', 'w')


# CHECK WITH STUART AND FINN HOW THEY SET UP 3D INPUT DATA
    ij = 0
    f.write ("%d \n" % (nx*ny*nz))
    time = model.time()
    f.write ("%16.8e \n" % (time/86400.))
    # f.write ("%16.8e \n" % (max(x1/time,z1/time)))
    # print(x)
    # print(time)
    for kk in range(nz):
        for jj in range(ny):
            for ii in range(nx):
                ij += 1
                velocity = (x[ii]/1.e5)/time
                # print(velocity)
                f.write ("%d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n" % (ij, velocity, np.log10(xartis[ii,jj,kk,0]), xignew[ii,jj,kk], xni56[ii,jj,kk], xco56[ii,jj,kk], xfe52[ii,jj,kk], xcr48[ii,jj,kk]))
                # f.write("%d %16.8e %16.8e %16.8e %16.8e \n" % (ij, velocity, np.log10(xartis[ii,jj,kk,0]), xignew[ii,jj,kk], xni56[ii,jj,kk]))
                # print(ij, velocity, np.log10(xartis[ii,jj,kk,0]), xni56[ii,jj,kk], xignew[ii,jj,kk])
                # h.write ("%d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n" % (ij, x[ii], xartis[ii,jj,kk,0], xignew[ii,jj,kk], xni56[ii,jj,kk], xco56[ii,jj,kk], xfe52[ii,jj,kk], xcr48[ii,jj,kk]))
                out = " ".join(("%16.8e " %dat) for dat in tuple(xartis[ii,jj,kk,1:]))
                g.write(("%d " % ij) + out + " \n")

    f.closed
    g.closed
    # h.closed
#
#
#
#     print('Writing HDF file...')
#     with h5py.File('s3.5_mapped.h5', 'w') as hf:
#        hf.create_dataset('model', data=xartis)
#     print('Wrote HDF file')

# # Check mass of mapped model
    dv0 = dx*dy*dz
    if ny ==1 and nz == 1:
        mc = np.arange(0, nx)
        dv0 = (4./3.)*np.pi*(3*mc**2+3*mc+1)*dx**3
    print("Mass of mapped model:",np.sum(xartis[:,:,:,0]*dv0[:,np.newaxis,np.newaxis])*solmassi,"M_sun")
    print("Iron group:          ",np.sum(xartis[:,:,:,0]*xignew[:,:,:]*dv0[:,np.newaxis,np.newaxis])*solmassi,"M_sun")
    print("Mapped mass/Inital mass:", np.sum(xartis[:,:,:,0]*dv0[:,np.newaxis,np.newaxis])*solmassi/initial_mass)
    print("Mapped IG mass/Inital IG mass:", np.sum(xartis[:,:,:,0]*xignew[:,:,:]*dv0[:,np.newaxis,np.newaxis])*solmassi/initial_IG_mass)
    # print(xartis[:,:,:,0])

    # if np.argwhere(np.isnan(xartis[:,:,:,0])) != []:
        # print('Problems in Zones ', np.argwhere(np.isnan(xartis[:,:,:,0])))
    print('Done!')
    # breakpoint()
    # return rhonew,xnunew,x,y,z,r
    # return np.sum(xartis[:,:,:,0]*dv0)*solmassi
