c Fortran Library for Skeleton 3D Electrostatic PIC Code
c written by Viktor K. Decyk, UCLA
c-----------------------------------------------------------------------
      subroutine DISTR3(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,idimp, 
     1nop,nx,ny,nz,ipbc)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny/nz = system length in x/y/z direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, npz, idimp, nop, nx, ny, nz, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
c local data
      integer j, k, l, k1, l1, npxy, npxyz
      real edgelx, edgely, edgelz, at1, at2, at3, at4, at5
      real sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
      npxyz = npxy*npz
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      at3 = real(nz)/real(npz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
         at3 = real(nz-2)/real(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 0.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      endif
c uniform density profile
      do 30 l = 1, npz
      l1 = npxy*(l - 1)
      at5 = edgelz + at3*(real(l) - .5)
      do 20 k = 1, npy
      k1 = npx*(k - 1) + l1
      at4 = edgely + at2*(real(k) - .5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - .5)
      part(2,j+k1) = at4
      part(3,j+k1) = at5
   10 continue
   20 continue
   30 continue
c maxwellian velocity distribution
      do 40 j = 1, npxyz
      part(4,j) = vtx*ranorm()
      part(5,j) = vty*ranorm()
      part(6,j) = vtz*ranorm()
   40 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 50 j = 1, npxyz
      dsum1 = dsum1 + part(4,j)
      dsum2 = dsum2 + part(5,j)
      dsum3 = dsum3 + part(6,j)
   50 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1.0/real(npxyz)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 60 j = 1, npxyz
      part(4,j) = part(4,j) - sum1
      part(5,j) = part(5,j) - sum2
      part(6,j) = part(6,j) - sum3
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GPUSH3L(part,fxyz,qbm,dt,ek,idimp,nop,nx,ny,nz,nxv,nyv,
     1nzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space
c scalar version using guard cells
c 94 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of field array, must be >= nx+1
c nyv = third dimension of field array, must be >= ny+1
c nzv = fourth dimension of field array, must be >= nz+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, nx, ny, nz, nxv, nyv, nzv, ipbc
      real qbm, dt, ek
      real part, fxyz
      dimension part(idimp,nop), fxyz(3,nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxp, dyp, dzp, amx, amy, amz, dx1, dx, dy, dz
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c find acceleration
      dx = amz*(amx*fxyz(1,nn,mm,ll) + amy*fxyz(1,np,mm,ll)             
     1        + dyp*fxyz(1,nn,mp,ll) + dx1*fxyz(1,np,mp,ll))            
     2   + dzp*(amx*fxyz(1,nn,mm,lp) + amy*fxyz(1,np,mm,lp)             
     3        + dyp*fxyz(1,nn,mp,lp) + dx1*fxyz(1,np,mp,lp))            
      dy = amz*(amx*fxyz(2,nn,mm,ll) + amy*fxyz(2,np,mm,ll)             
     1        + dyp*fxyz(2,nn,mp,ll) + dx1*fxyz(2,np,mp,ll))            
     2   + dzp*(amx*fxyz(2,nn,mm,lp) + amy*fxyz(2,np,mm,lp)             
     3        + dyp*fxyz(2,nn,mp,lp) + dx1*fxyz(2,np,mp,lp))            
      dz = amz*(amx*fxyz(3,nn,mm,ll) + amy*fxyz(3,np,mm,ll)             
     1        + dyp*fxyz(3,nn,mp,ll) + dx1*fxyz(3,np,mp,ll))            
     2   + dzp*(amx*fxyz(3,nn,mm,lp) + amy*fxyz(3,np,mm,lp)             
     3        + dyp*fxyz(3,nn,mp,lp) + dx1*fxyz(3,np,mp,lp))            
c new velocity
      dx = part(4,j) + qtm*dx
      dy = part(5,j) + qtm*dy
      dz = part(6,j) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j))**2 + (dy + part(5,j))**2           
     1            + (dz + part(6,j))**2
      part(4,j) = dx
      part(5,j) = dy
      part(6,j) = dz
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
      dz = part(3,j) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j)
            part(6,j) = -part(6,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
      part(3,j) = dz
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine GPOST3L(part,q,qm,nop,idimp,nxv,nyv,nzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 6
c nxv = first dimension of charge array, must be >= nx+1
c nyv = second dimension of charge array, must be >= ny+1
c nzv = third dimension of charge array, must be >= nz+1
      implicit none
      integer nop, idimp, nxv, nyv, nzv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxv,nyv,nzv)
c local data
      integer j, nn, mm, ll, np, mp, lp
      real dx1, dxp, dyp, dzp, amx, amy, amz
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      ll = part(3,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      dzp = part(3,j) - real(ll)
      nn = nn + 1
      amx = qm - dxp
      amy = 1.0 - dyp
      np = nn + 1
      mm = mm + 1
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1.0 - dzp
      ll = ll + 1
      amy = dxp*amy
      lp = ll + 1
c deposit charge
      q(nn,mm,ll) = q(nn,mm,ll) + amx*amz
      q(np,mm,ll) = q(np,mm,ll) + amy*amz
      q(nn,mp,ll) = q(nn,mp,ll) + dyp*amz
      q(np,mp,ll) = q(np,mp,ll) + dx1*amz
      q(nn,mm,lp) = q(nn,mm,lp) + amx*dzp
      q(np,mm,lp) = q(np,mm,lp) + amy*dzp
      q(nn,mp,lp) = q(nn,mp,lp) + dyp*dzp
      q(np,mp,lp) = q(np,mp,lp) + dx1*dzp
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DSORTP3YZL(parta,partb,npic,idimp,nop,ny1,nyz1)
c this subroutine sorts particles by y,z grid
c linear interpolation
c part = particle array
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 6
c nop = number of particles
c ny1 = system length in y direction + 1
c nyz1 = ny1*nz1, where nz1 = system length in z direction + 1
      implicit none
      integer npic, idimp, nop, ny1, nyz1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(nyz1)
c local data
      integer i, j, k, m, l, isum, ist, ip
c clear counter array
      do 10 k = 1, nyz1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j) + 1.0
      l = parta(3,j)
      l = m + ny1*l
      npic(l) = npic(l) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyz1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(2,j) + 1.0
      l = parta(3,j)
      l = m + ny1*l
      ip = npic(l) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(l) = ip
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine CGUARD3L(fxyz,nx,ny,nz,nxe,nye,nze)
c replicate extended periodic vector field fxyz
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real fxyz
      integer nx, ny, nz, nxe, nye, nze
      dimension fxyz(3,nxe,nye,nze)
c local data
      integer j, k, l
      do 30 l = 1, nz
      do 10 k = 1, ny
      fxyz(1,nx+1,k,l) = fxyz(1,1,k,l)
      fxyz(2,nx+1,k,l) = fxyz(2,1,k,l)
      fxyz(3,nx+1,k,l) = fxyz(3,1,k,l)
   10 continue
      do 20 j = 1, nx
      fxyz(1,j,ny+1,l) = fxyz(1,j,1,l)
      fxyz(2,j,ny+1,l) = fxyz(2,j,1,l)
      fxyz(3,j,ny+1,l) = fxyz(3,j,1,l)
   20 continue
      fxyz(1,nx+1,ny+1,l) = fxyz(1,1,1,l)
      fxyz(2,nx+1,ny+1,l) = fxyz(2,1,1,l)
      fxyz(3,nx+1,ny+1,l) = fxyz(3,1,1,l)
   30 continue
      do 50 k = 1, ny
      do 40 j = 1, nx
      fxyz(1,j,k,nz+1) = fxyz(1,j,k,1)
      fxyz(2,j,k,nz+1) = fxyz(2,j,k,1)
      fxyz(3,j,k,nz+1) = fxyz(3,j,k,1)
   40 continue
      fxyz(1,nx+1,k,nz+1) = fxyz(1,1,k,1)
      fxyz(2,nx+1,k,nz+1) = fxyz(2,1,k,1)
      fxyz(3,nx+1,k,nz+1) = fxyz(3,1,k,1)
   50 continue
      do 60 j = 1, nx
      fxyz(1,j,ny+1,nz+1) = fxyz(1,j,1,1)
      fxyz(2,j,ny+1,nz+1) = fxyz(2,j,1,1)
      fxyz(3,j,ny+1,nz+1) = fxyz(3,j,1,1)
   60 continue
      fxyz(1,nx+1,ny+1,nz+1) = fxyz(1,1,1,1)
      fxyz(2,nx+1,ny+1,nz+1) = fxyz(2,1,1,1)
      fxyz(3,nx+1,ny+1,nz+1) = fxyz(3,1,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine AGUARD3L(q,nx,ny,nz,nxe,nye,nze)
c accumulate extended periodic scalar field q
c linear interpolation
c nx/ny/nz = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nze = third dimension of field arrays, must be >= nz+1
      implicit none
      real q
      integer nx, ny, nz, nxe, nye, nze
      dimension q(nxe,nye,nze)
      integer j, k, l
c accumulate edges of extended field
      do 30 l = 1, nz
      do 10 k = 1, ny
      q(1,k,l) = q(1,k,l) + q(nx+1,k,l)
      q(nx+1,k,l) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1,l) = q(j,1,l) + q(j,ny+1,l)
      q(j,ny+1,l) = 0.0
   20 continue
      q(1,1,l) = q(1,1,l) + q(nx+1,ny+1,l)
      q(nx+1,ny+1,l) = 0.0
   30 continue
      do 50 k = 1, ny
      do 40 j = 1, nx
      q(j,k,1) = q(j,k,1) + q(j,k,nz+1)
      q(j,k,nz+1) = 0.0
   40 continue
      q(1,k,1) = q(1,k,1) + q(nx+1,k,nz+1)
      q(nx+1,k,nz+1) = 0.0
   50 continue
      do 60 j = 1, nx
      q(j,1,1) = q(j,1,1) + q(j,ny+1,nz+1)
      q(j,ny+1,nz+1) = 0.0
   60 continue
      q(1,1,1) = q(1,1,1) + q(nx+1,ny+1,nz+1)
      q(nx+1,ny+1,nz+1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine POIS33(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxvh,
     1nyv,nzv,nxhd,nyhd,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, output: ffc
c input: isign,ax,ay,az,affp,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c for isign = -1, output: fxyz, we
c input: q,ffc,isign,nx,ny,nz,nxvh,nyv,nzv,nxhd,nyhd,nzhd
c approximate flop count is:
c 59*nxc*nyc*nzc + 26*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1, nzc = nz/2 - 1
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c q(j,k,l) = complex charge density for fourier mode (j-1,k-1,l-1)
c fxyz(1,j,k,l) = x component of complex force/charge
c fxyz(2,j,k,l) = y component of complex force/charge
c fxyz(3,j,k,l) = z component of complex force/charge
c all for fourier mode (j-1,k-1,l-1)
c aimag(ffc(j,k,l)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1,l-1)
c real(ffc(j,k,l)) = potential green's function g
c for fourier mode (j-1,k-1,l-1)
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nzv = third dimension of field arrays, must be >= nz
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
c nzhd = third dimension of form factor array, must be >= nzh
      implicit none
      integer isign, nx, ny, nz, nxvh, nyv, nzv, nxhd, nyhd, nzhd
      real ax, ay, az, affp, we
      complex q, fxyz, ffc
      dimension q(nxvh,nyv,nzv), fxyz(3,nxvh,nyv,nzv)
      dimension ffc(nxhd,nyhd,nzhd)
c local data
      integer nxh, nyh, nzh, ny2, nz2, j, k, l, k1, l1
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      nzh = max(1,nz/2)
      ny2 = ny + 2
      nz2 = nz + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dnz = 6.28318530717959/real(nz)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 40
c prepare form factor array
      do 30 l = 1, nzh
      dkz = dnz*real(l - 1)
      at1 = dkz*dkz
      at2 = (dkz*az)**2
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = (dky*ay)**2 + at2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at5 = dkx*dkx + at3
      at6 = exp(-.5*((dkx*ax)**2 + at4))
      if (at5.eq.0.) then
         ffc(j,k,l) = cmplx(affp,1.0)
      else
         ffc(j,k,l) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      do 90 l = 2, nzh
      l1 = nz2 - l
      dkz = dnz*real(l - 1)
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 50 j = 2, nxh
      at1 = real(ffc(j,k,l))*aimag(ffc(j,k,l))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(j,k,l)),-real(q(j,k,l)))
      zt2 = cmplx(aimag(q(j,k1,l)),-real(q(j,k1,l)))
      fxyz(1,j,k,l) = at2*zt1
      fxyz(2,j,k,l) = at3*zt1
      fxyz(3,j,k,l) = at4*zt1
      fxyz(1,j,k1,l) = at2*zt2
      fxyz(2,j,k1,l) = -at3*zt2
      fxyz(3,j,k1,l) = at4*zt2
      zt1 = cmplx(aimag(q(j,k,l1)),-real(q(j,k,l1)))
      zt2 = cmplx(aimag(q(j,k1,l1)),-real(q(j,k1,l1)))
      fxyz(1,j,k,l1) = at2*zt1
      fxyz(2,j,k,l1) = at3*zt1
      fxyz(3,j,k,l1) = -at4*zt1
      fxyz(1,j,k1,l1) = at2*zt2
      fxyz(2,j,k1,l1) = -at3*zt2
      fxyz(3,j,k1,l1) = -at4*zt2
      wp = wp + at1*(q(j,k,l)*conjg(q(j,k,l))                           
     1   + q(j,k1,l)*conjg(q(j,k1,l)) + q(j,k,l1)*conjg(q(j,k,l1))      
     2   + q(j,k1,l1)*conjg(q(j,k1,l1)))
   50 continue
   60 continue
c mode numbers kx = 0, nx/2
      do 70 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k,l))*aimag(ffc(1,k,l))
      at3 = dny*real(k - 1)*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(1,k,l)),-real(q(1,k,l)))
      zt2 = cmplx(aimag(q(1,k,l1)),-real(q(1,k,l1)))
      fxyz(1,1,k,l) = zero
      fxyz(2,1,k,l) = at3*zt1
      fxyz(3,1,k,l) = at4*zt1
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = at3*zt2
      fxyz(3,1,k,l1) = -at4*zt2
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      wp = wp + at1*(q(1,k,l)*conjg(q(1,k,l))                           
     1   + q(1,k,l1)*conjg(q(1,k,l1)))
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 80 j = 2, nxh
      at1 = real(ffc(j,1,l))*aimag(ffc(j,1,l))
      at2 = dnx*real(j - 1)*at1
      at4 = dkz*at1
      zt1 = cmplx(aimag(q(j,1,l)),-real(q(j,1,l)))
      zt2 = cmplx(aimag(q(j,1,l1)),-real(q(j,1,l1)))
      fxyz(1,j,1,l) = at2*zt1
      fxyz(2,j,1,l) = zero
      fxyz(3,j,1,l) = at4*zt1
      fxyz(1,j,k1,l) = zero
      fxyz(2,j,k1,l) = zero
      fxyz(3,j,k1,l) = zero
      fxyz(1,j,1,l1) = at2*zt2
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = -at4*zt2
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      wp = wp + at1*(q(j,1,l)*conjg(q(j,1,l))                           
     1   + q(j,1,l1)*conjg(q(j,1,l1)))
   80 continue
c mode numbers kx = 0, nx/2
      at1 = real(ffc(1,1,l))*aimag(ffc(1,1,l))
      at4 = dkz*at1
      fxyz(1,1,1,l) = zero
      fxyz(2,1,1,l) = zero
      fxyz(3,1,1,l) = at4*cmplx(aimag(q(1,1,l)),-real(q(1,1,l)))
      fxyz(1,1,k1,l) = zero
      fxyz(2,1,k1,l) = zero
      fxyz(3,1,k1,l) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      wp = wp + at1*(q(1,1,l)*conjg(q(1,1,l)))
   90 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 110 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 100 j = 2, nxh
      at1 = real(ffc(j,k,1))*aimag(ffc(j,k,1))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k,1)),-real(q(j,k,1)))
      zt2 = cmplx(aimag(q(j,k1,1)),-real(q(j,k1,1)))
      fxyz(1,j,k,1) = at2*zt1
      fxyz(2,j,k,1) = at3*zt1
      fxyz(3,j,k,1) = zero
      fxyz(1,j,k1,1) = at2*zt2
      fxyz(2,j,k1,1) = -at3*zt2
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,k,l1) = zero
      fxyz(2,j,k,l1) = zero
      fxyz(3,j,k,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      wp = wp + at1*(q(j,k,1)*conjg(q(j,k,1))
     1   + q(j,k1,1)*conjg(q(j,k1,1)))
  100 continue
  110 continue
c mode numbers kx = 0, nx/2
      do 120 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k,1))*aimag(ffc(1,k,1))
      at3 = dny*real(k - 1)*at1
      fxyz(1,1,k,1) = zero
      fxyz(2,1,k,1) = at3*cmplx(aimag(q(1,k,1)),-real(q(1,k,1)))
      fxyz(3,1,k,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,k,l1) = zero
      fxyz(2,1,k,l1) = zero
      fxyz(3,1,k,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      wp = wp + at1*(q(1,k,1)*conjg(q(1,k,1)))
  120 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 130 j = 2, nxh
      at1 = real(ffc(j,1,1))*aimag(ffc(j,1,1))
      at2 = dnx*real(j - 1)*at1
      fxyz(1,j,1,1) = at2*cmplx(aimag(q(j,1,1)),-real(q(j,1,1)))
      fxyz(2,j,1,1) = zero
      fxyz(3,j,1,1) = zero
      fxyz(1,j,k1,1) = zero
      fxyz(2,j,k1,1) = zero
      fxyz(3,j,k1,1) = zero
      fxyz(1,j,1,l1) = zero
      fxyz(2,j,1,l1) = zero
      fxyz(3,j,1,l1) = zero
      fxyz(1,j,k1,l1) = zero
      fxyz(2,j,k1,l1) = zero
      fxyz(3,j,k1,l1) = zero
      wp = wp + at1*(q(j,1,1)*conjg(q(j,1,1)))
  130 continue
      fxyz(1,1,1,1) = zero
      fxyz(2,1,1,1) = zero
      fxyz(3,1,1,1) = zero
      fxyz(1,1,k1,1) = zero
      fxyz(2,1,k1,1) = zero
      fxyz(3,1,k1,1) = zero
      fxyz(1,1,1,l1) = zero
      fxyz(2,1,1,l1) = zero
      fxyz(3,1,1,l1) = zero
      fxyz(1,1,k1,l1) = zero
      fxyz(2,1,k1,l1) = zero
      fxyz(3,1,k1,l1) = zero
      we = real(nx)*real(ny)*real(nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
c this subroutine calculates tables needed by a three dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, indz, nxhyzd, nxyzhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, indz, nxhyzd, nxyzhd
      integer mixup
      complex sct
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh
      integer j, k, lb, ll, jb, it
      real dnxyz, arg
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 10 k = 1, ndx1yz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/real(nxyz)
      do 30 j = 1, nxyzh
      arg = dnxyz*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3RX(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd, 
     1nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WFFT3R3(f,isign,mixup,sct,indx,indy,indz,nxhd,nyd,nzd, 
     1nxhyzd,nxyzhd)
c wrapper function for 3 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, indz, nxhd, nyd, nzd, nxhyzd, nxyzhd
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer ny, nz, nyi, nzi
      data nyi, nzi /1,1/
c calculate range of indices
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform xy fft
         call FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
c perform z fft
         call FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform z fft
         call FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,ny,nxhd,nyd,n
     1zd,nxhyzd,nxyzhd)
c perform xy fft
         call FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nz,nxhd,nyd,
     1nzd,nxhyzd,nxyzhd)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RXY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd, 
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of z,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(n,m,i) = (1/nx*ny*nz)*sum(f(j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)*
c       exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k,l) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nzt, nrx, nry
      integer i, j, k, l, n, j1, j2, k1, k2, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 180
c inverse fourier transform
      do 170 n = nzi, nzt
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 i = 1, ny
      t1 = f(j1,i,n)
      f(j1,i,n) = f(j,i,n)
      f(j,i,n) = t1
   10 continue
   20 continue
c first transform in x
      nrx = nxyz/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = ani*(t1 + t2)
      f(nxh2-j,k,n) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = 1, ny
      f(nxhh+1,k,n) = ani*conjg(f(nxhh+1,k,n))
      f(1,k,n) = ani*cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                     real(f(1,k,n)) - aimag(f(1,k,n)))
   90 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      do 100 i = 1, nxh
      t1 = f(i,k1,n)
      f(i,k1,n) = f(i,k,n)
      f(i,k,n) = t1
  100 continue
  110 continue
c then transform in y
      nry = nxyz/ny
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 120 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  120 continue
  130 continue
  140 continue
  150 continue
c unscramble modes kx = 0, nx/2
      do 160 k = 2, nyh
      t1 = f(1,ny2-k,n)
      f(1,ny2-k,n) = 0.5*cmplx(aimag(f(1,k,n) + t1),real(f(1,k,n) - t1))
      f(1,k,n) = 0.5*cmplx(real(f(1,k,n) + t1),aimag(f(1,k,n) - t1))
  160 continue
  170 continue
      return
c forward fourier transform
  180 do 350 n = nzi, nzt
c scramble modes kx = 0, nx/2
      do 190 k = 2, nyh
      t1 = cmplx(aimag(f(1,ny2-k,n)),real(f(1,ny2-k,n)))
      f(1,ny2-k,n) = conjg(f(1,k,n) - t1)
      f(1,k,n) = f(1,k,n) + t1
  190 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 210 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 210
      do 200 i = 1, nxh
      t1 = f(i,k1,n)
      f(i,k1,n) = f(i,k,n)
      f(i,k,n) = t1
  200 continue
  210 continue
c then transform in y
      nry = nxyz/ny
      do 250 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = 1, nxh
      t2 = t1*f(i,j2,n)
      f(i,j2,n) = f(i,j1,n) - t2
      f(i,j1,n) = f(i,j1,n) + t2
  220 continue
  230 continue
  240 continue
  250 continue
c scramble coefficients
      kmr = nxyz/nx
      do 270 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 260 k = 1, ny
      t2 = conjg(f(nxh2-j,k,n))
      t1 = f(j,k,n) + t2
      t2 = (f(j,k,n) - t2)*t3
      f(j,k,n) = t1 + t2
      f(nxh2-j,k,n) = conjg(t1 - t2)
  260 continue
  270 continue
      do 280 k = 1, ny
      f(nxhh+1,k,n) = 2.0*conjg(f(nxhh+1,k,n))
      f(1,k,n) = cmplx(real(f(1,k,n)) + aimag(f(1,k,n)),
     1                 real(f(1,k,n)) - aimag(f(1,k,n)))
  280 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 300 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 300
      do 290 i = 1, ny
      t1 = f(j1,i,n)
      f(j1,i,n) = f(j,i,n)
      f(j,i,n) = t1
  290 continue
  300 continue
c finally transform in x
      nrx = nxyz/nxh
      do 340 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 330 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 320 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 310 i = 1, ny
      t2 = t1*f(j2,i,n)
      f(j2,i,n) = f(j1,i,n) - t2
      f(j1,i,n) = f(j1,i,n) + t2
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3RXZ(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd, 
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(j,k,l) = sum(f(j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, a forward fourier transform is performed
c f(n,m,i) = sum(f(n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f
c nyd,nzd = second and third dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(j,k,l)= real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1,1,1)) = real part of mode nx/2,0,0
c imag(f(1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz
      integer i, j, k, l, n, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 30 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 30
      do 20 n = nyi, nyt
      do 10 i = 1, nxh
      t1 = f(i,n,l1)
      f(i,n,l1) = f(i,n,l)
      f(i,n,l) = t1
   10 continue
   20 continue
   30 continue
c finally transform in z
      nrz = nxyz/nz
      do 80 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 n = nyi, nyt
      do 40 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 90 n = 2, nzh
      if (nyi.eq.1) then
         t1 = f(1,1,nz2-n)
         f(1,1,nz2-n) = 0.5*cmplx(aimag(f(1,1,n) + t1),
     1                            real(f(1,1,n) - t1))
         f(1,1,n) = 0.5*cmplx(real(f(1,1,n) + t1),aimag(f(1,1,n) - t1))
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = f(1,nyh+1,nz2-n)
         f(1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(1,nyh+1,n) + t1),
     1                                real(f(1,nyh+1,n) - t1))
         f(1,nyh+1,n) = 0.5*cmplx(real(f(1,nyh+1,n) + t1),
     1                            aimag(f(1,nyh+1,n) - t1))
      endif
   90 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  100 do 110 n = 2, nzh
      if (nyi.eq.1) then
         t1 = cmplx(aimag(f(1,1,nz2-n)),real(f(1,1,nz2-n)))
         f(1,1,nz2-n) = conjg(f(1,1,n) - t1)
         f(1,1,n) = f(1,1,n) + t1
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         t1 = cmplx(aimag(f(1,nyh+1,nz2-n)),real(f(1,nyh+1,nz2-n)))
         f(1,nyh+1,nz2-n) = conjg(f(1,nyh+1,n) - t1)
         f(1,nyh+1,n) = f(1,nyh+1,n) + t1
      endif
  110 continue
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 140 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 140
      do 130 n = nyi, nyt
      do 120 i = 1, nxh
      t1 = f(i,n,l1)
      f(i,n,l1) = f(i,n,l)
      f(i,n,l) = t1
  120 continue
  130 continue
  140 continue
c first transform in z
      nrz = nxyz/nz
      do 190 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 n = nyi, nyt
      do 150 i = 1, nxh
      t2 = t1*f(i,n,j2)
      f(i,n,j2) = f(i,n,j1) - t2
      f(i,n,j1) = f(i,n,j1) + t2
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3R3XY(f,isign,mixup,sct,indx,indy,indz,nzi,nzp,nxhd,
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the x-y part of 3 three dimensional complex
c to real fast fourier transforms and their inverses, for a subset of z,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:3,n,m,i) = (1/nx*ny*nz)*sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*n*j/nx)
c       *exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,j,k,l) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nzi = initial z index used
c nzp = number of z indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l) = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nzi, nzp, nxhd, nyd, nzd
      integer nxhyzd,nxyzhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nzt, nrx, nry
      integer i, j, k, l, n, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nzt = nzi + nzp - 1
      if (isign.gt.0) go to 230
c inverse fourier transform
      do 220 n = nzi, nzt
c swap complex components
      do 20 i = 1, ny
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(real(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = aimag(f(2,j,i,n))
      f(2,j,i,n) = cmplx(aimag(f(1,j,i,n)),at1)
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 i = 1, ny
      t1 = f(1,j1,i,n)
      t2 = f(2,j1,i,n)
      t3 = f(3,j1,i,n)
      f(1,j1,i,n) = f(1,j,i,n)
      f(2,j1,i,n) = f(2,j,i,n)
      f(3,j1,i,n) = f(3,j,i,n)
      f(1,j,i,n) = t1
      f(2,j,i,n) = t2
      f(3,j,i,n) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxyz/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(real(nx)*real(ny)*real(nz))
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = 1, ny
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = ani*(t1 + t2)
      f(jj,nxh2-j,k,n) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.0*ani
      do 130 k = 1, ny
      do 120 jj = 1, 3
      f(jj,nxhh+1,k,n) = ani*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = ani*cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                        real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  120 continue
  130 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 150 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 150
      do 140 i = 1, nxh
      t1 = f(1,i,k1,n)
      t2 = f(2,i,k1,n)
      t3 = f(3,i,k1,n)
      f(1,i,k1,n) = f(1,i,k,n)
      f(2,i,k1,n) = f(2,i,k,n)
      f(3,i,k1,n) = f(3,i,k,n)
      f(1,i,k,n) = t1
      f(2,i,k,n) = t2
      f(3,i,k,n) = t3
  140 continue
  150 continue
c then transform in y
      nry = nxyz/ny
      do 190 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 160 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  160 continue
  170 continue
  180 continue
  190 continue
c unscramble modes kx = 0, nx/2
      do 210 k = 2, nyh
      do 200 jj = 1, 3
      t1 = f(jj,1,ny2-k,n)
      f(jj,1,ny2-k,n) = 0.5*cmplx(aimag(f(jj,1,k,n) + t1),
     1                            real(f(jj,1,k,n) - t1))
      f(jj,1,k,n) = 0.5*cmplx(real(f(jj,1,k,n) + t1),
     1                        aimag(f(jj,1,k,n) - t1))
  200 continue
  210 continue
  220 continue
      return
c forward fourier transform
  230 do 450 n = nzi, nzt
c scramble modes kx = 0, nx/2
      do 250 k = 2, nyh
      do 240 jj = 1, 3
      t1 = cmplx(aimag(f(jj,1,ny2-k,n)),real(f(jj,1,ny2-k,n)))
      f(jj,1,ny2-k,n) = conjg(f(jj,1,k,n) - t1)
      f(jj,1,k,n) = f(jj,1,k,n) + t1
  240 continue
  250 continue
c bit-reverse array elements in y
      nry = nxhyz/ny
      do 270 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 270
      do 260 i = 1, nxh
      t1 = f(1,i,k1,n)
      t2 = f(2,i,k1,n)
      t3 = f(3,i,k1,n)
      f(1,i,k1,n) = f(1,i,k,n)
      f(2,i,k1,n) = f(2,i,k,n)
      f(3,i,k1,n) = f(3,i,k,n)
      f(1,i,k,n) = t1
      f(2,i,k,n) = t2
      f(3,i,k,n) = t3
  260 continue
  270 continue
c then transform in y
      nry = nxyz/ny
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 280 i = 1, nxh
      t2 = t1*f(1,i,j2,n)
      t3 = t1*f(2,i,j2,n)
      t4 = t1*f(3,i,j2,n)
      f(1,i,j2,n) = f(1,i,j1,n) - t2
      f(2,i,j2,n) = f(2,i,j1,n) - t3
      f(3,i,j2,n) = f(3,i,j1,n) - t4
      f(1,i,j1,n) = f(1,i,j1,n) + t2
      f(2,i,j1,n) = f(2,i,j1,n) + t3
      f(3,i,j1,n) = f(3,i,j1,n) + t4
  280 continue
  290 continue
  300 continue
  310 continue
c scramble coefficients
      kmr = nxyz/nx
      do 340 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 330 k = 1, ny
      do 320 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k,n))
      t1 = f(jj,j,k,n) + t2
      t2 = (f(jj,j,k,n) - t2)*t3
      f(jj,j,k,n) = t1 + t2
      f(jj,nxh2-j,k,n) = conjg(t1 - t2)
  320 continue
  330 continue
  340 continue
      do 360 k = 1, ny
      do 350 jj = 1, 3
      f(jj,nxhh+1,k,n) = 2.0*conjg(f(jj,nxhh+1,k,n))
      f(jj,1,k,n) = cmplx(real(f(jj,1,k,n)) + aimag(f(jj,1,k,n)),
     1                    real(f(jj,1,k,n)) - aimag(f(jj,1,k,n)))
  350 continue
  360 continue
c bit-reverse array elements in x
      nrx = nxhyz/nxh
      do 380 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 380
      do 370 i = 1, ny
      t1 = f(1,j1,i,n)
      t2 = f(2,j1,i,n)
      t3 = f(3,j1,i,n)
      f(1,j1,i,n) = f(1,j,i,n)
      f(2,j1,i,n) = f(2,j,i,n)
      f(3,j1,i,n) = f(3,j,i,n)
      f(1,j,i,n) = t1
      f(2,j,i,n) = t2
      f(3,j,i,n) = t3
  370 continue
  380 continue
c finally transform in x
      nrx = nxyz/nxh
      do 420 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 410 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 400 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 390 i = 1, ny
      t2 = t1*f(1,j2,i,n)
      t3 = t1*f(2,j2,i,n)
      t4 = t1*f(3,j2,i,n)
      f(1,j2,i,n) = f(1,j1,i,n) - t2
      f(2,j2,i,n) = f(2,j1,i,n) - t3
      f(3,j2,i,n) = f(3,j1,i,n) - t4
      f(1,j1,i,n) = f(1,j1,i,n) + t2
      f(2,j1,i,n) = f(2,j1,i,n) + t3
      f(3,j1,i,n) = f(3,j1,i,n) + t4
  390 continue
  400 continue
  410 continue
  420 continue
c swap complex components
      do 440 i = 1, ny
      do 430 j = 1, nxh
      at1 = real(f(3,j,i,n))
      f(3,j,i,n) = cmplx(aimag(f(2,j,i,n)),aimag(f(3,j,i,n)))
      at2 = real(f(2,j,i,n))
      f(2,j,i,n) = cmplx(at1,aimag(f(1,j,i,n)))
      f(1,j,i,n) = cmplx(real(f(1,j,i,n)),at2)
  430 continue
  440 continue
  450 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FFT3R3Z(f,isign,mixup,sct,indx,indy,indz,nyi,nyp,nxhd, 
     1nyd,nzd,nxhyzd,nxyzhd)
c this subroutine performs the z part of 3 three dimensional complex to
c real fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny*nz
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, three inverse fourier transforms is performed
c f(1:3,j,k,l) = sum(f(1:3,j,k,i)*exp(-sqrt(-1)*2pi*l*i/nz))
c if isign = 1, three forward fourier transforms are performed
c f(1:3,n,m,i) = sum(f(1:3,n,m,l)*exp(sqrt(-1)*2pi*l*i/nz))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd,nzd = third and fourth dimensions of f
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = maximum of (nx,ny,nz)/2
c fourier coefficients are stored as follows:
c f(1:3,j,k,l) = real, imaginary part of mode j-1,k-1,l-1
c where 1 <= j <= nx/2, 1 <= k <= ny, 1 <= l <= nz, except for
c f(1:3,1,k,l), = real, imaginary part of mode nx/2,k-1,l-1,
c where ny/2+2 <= k <= ny and 1 <= l <= nz, and
c f(1:3,1,1,l) = real, imaginary part of mode nx/2,0,l-1,
c f(1:3,1,ny/2+1,l) = real, imaginary part mode nx/2,ny/2,l-1,
c where nz/2+2 <= l <= nz, and
c imag(f(1:3,1,1,1)) = real part of mode nx/2,0,0
c imag(f(1:3,1,ny/2+1,1)) = real part of mode nx/2,ny/2,0
c imag(f(1:3,1,1,nz/2+1)) = real part of mode nx/2,0,nz/2
c imag(f(1:3,1,ny/2+1,nz/2+1)) = real part of mode nx/2,ny/2,nz/2
c using jpl storage convention, as described in:
c E. Huang, P. C. Liewer, V. K. Decyk, and R. D. Ferraro, "Concurrent
c Three-Dimensional Fast Fourier Transform Algorithms for Coarse-Grained
c Distributed Memory Parallel Computers," Caltech CRPC Report 217-50,
c December 1993.
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, indz, nyi, nyp, nxhd, nyd, nzd
      integer nxhyzd, nxyzhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd,nzd), mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2
      integer nz, nzh, nz2, nxyz, nxhyz, nyt, nrz
      integer i, j, k, l, n, jj, j1, j2, k1, k2, l1, ns, ns2, km, kmr
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 30 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 30
      do 20 n = nyi, nyt
      do 10 i = 1, nxh
      t1 = f(1,i,n,l1)
      t2 = f(2,i,n,l1)
      t3 = f(3,i,n,l1)
      f(1,i,n,l1) = f(1,i,n,l)
      f(2,i,n,l1) = f(2,i,n,l)
      f(3,i,n,l1) = f(3,i,n,l)
      f(1,i,n,l) = t1
      f(2,i,n,l) = t2
      f(3,i,n,l) = t3
   10 continue
   20 continue
   30 continue
c finally transform in z
      nrz = nxyz/nz
      do 80 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 n = nyi, nyt
      do 40 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 110 n = 2, nzh
      if (nyi.eq.1) then
         do 90 jj = 1, 3
         t1 = f(jj,1,1,nz2-n)
         f(jj,1,1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,1,n) + t1),
     1                               real(f(jj,1,1,n) - t1))
         f(jj,1,1,n) = 0.5*cmplx(real(f(jj,1,1,n) + t1),
     1                           aimag(f(jj,1,1,n) - t1))
   90    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 100 jj = 1, 3
         t1 = f(jj,1,nyh+1,nz2-n)
         f(jj,1,nyh+1,nz2-n) = 0.5*cmplx(aimag(f(jj,1,nyh+1,n) + t1),
     1                                  real(f(jj,1,nyh+1,n) - t1))
         f(jj,1,nyh+1,n) = 0.5*cmplx(real(f(jj,1,nyh+1,n) + t1),
     1                              aimag(f(jj,1,nyh+1,n) - t1))
  100    continue
      endif
  110 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  120 do 150 n = 2, nzh
      if (nyi.eq.1) then
         do 130 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,1,nz2-n)),real(f(jj,1,1,nz2-n)))
         f(jj,1,1,nz2-n) = conjg(f(jj,1,1,n) - t1)
         f(jj,1,1,n) = f(jj,1,1,n) + t1
  130    continue
      endif
      if ((nyi.le.nyh+1).and.(nyt.ge.nyh+1)) then
         do 140 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,nyh+1,nz2-n)),
     1              real(f(jj,1,nyh+1,nz2-n)))
         f(jj,1,nyh+1,nz2-n) = conjg(f(jj,1,nyh+1,n) - t1)
         f(jj,1,nyh+1,n) = f(jj,1,nyh+1,n) + t1
  140    continue
      endif
  150 continue
c bit-reverse array elements in z
      nrz = nxhyz/nz
      do 180 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 180
      do 170 n = nyi, nyt
      do 160 i = 1, nxh
      t1 = f(1,i,n,l1)
      t2 = f(2,i,n,l1)
      t3 = f(3,i,n,l1)
      f(1,i,n,l1) = f(1,i,n,l)
      f(2,i,n,l1) = f(2,i,n,l)
      f(3,i,n,l1) = f(3,i,n,l)
      f(1,i,n,l) = t1
      f(2,i,n,l) = t2
      f(3,i,n,l) = t3
  160 continue
  170 continue
  180 continue
c first transform in z
      nrz = nxyz/nz
      do 230 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 200 n = nyi, nyt
      do 190 i = 1, nxh
      t2 = t1*f(1,i,n,j2)
      t3 = t1*f(2,i,n,j2)
      t4 = t1*f(3,i,n,j2)
      f(1,i,n,j2) = f(1,i,n,j1) - t2
      f(2,i,n,j2) = f(2,i,n,j1) - t3
      f(3,i,n,j2) = f(3,i,n,j1) - t4
      f(1,i,n,j1) = f(1,i,n,j1) + t2
      f(2,i,n,j1) = f(2,i,n,j1) + t3
      f(3,i,n,j1) = f(3,i,n,j1) + t4
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GSPOST3L(part,q,qm,nop,idimp,nxv,nyv,nxyzv)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c q(j,k,l) = charge density at grid point j,k,l
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 6
c nxv = first virtual dimension of charge array, must be >= nx+1
c nyv = second virtual dimension of charge array, must be >= ny+1
c nxyzv = actual dimension of charge array, must be >= nxv*nyv*(nz+1)
      implicit none
      integer nop, idimp, nxv, nyv, nxyzv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxyzv)
c local data
      integer j, nnn, mmm, lll, nxyv, nn, mm, ll, mp, lp
      real dxn, dyn, dzn, dxp, dyp, dzp, amx, amy, amz, dx1, dx2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmm = part(2,1)
      lll = part(3,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmm)
      dzn = part(3,1) - real(lll)
      nxyv = nxv*nyv
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      nnn = part(1,j)
      mmm = part(2,j)
      lll = part(3,j)
      dxp = qm*dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j) - real(nnn)
      dyn = part(2,j) - real(mmm)
      dzn = part(3,j) - real(lll)
      amx = qm - dxp
      amy = 1.0 - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzp
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c deposit charge
      dxp = q(mm) + amx*amz
      dx2 = q(mm+1) + amy*amz
      dx3 = q(mp) + dyp*amz
      amz = q(mp+1) + dx1*amz
      amx = q(ll) + amx*dzp
      amy = q(ll+1) + amy*dzp
      dyp = q(lp) + dyp*dzp
      dzp = q(lp+1) + dx1*dzp
      q(mm) = dxp
      q(mm+1) = dx2
      q(mp) = dx3
      q(mp+1) = amz
      q(ll) = amx
      q(ll+1) = amy
      q(lp) = dyp
      q(lp+1) = dzp
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      dxp = qm*dxn
      amx = qm - dxp
      amy = 1.0 - dyn
      mm = mm + nn
      dx1 = dxp*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzn
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
      q(mm) = q(mm) + amx*amz
      q(mm+1) = q(mm+1) + amy*amz
      q(mp) = q(mp) + dyp*amz
      q(mp+1) = q(mp+1) + dx1*amz
      q(ll) = q(ll) + amx*dzn
      q(ll+1) = q(ll+1) + amy*dzn
      q(lp) = q(lp) + dyp*dzn
      q(lp+1) = q(lp+1) + dx1*dzn
      return
      end
c-----------------------------------------------------------------------
      subroutine GSPUSH3L(part,fxyz,qbm,dt,ek,idimp,nop,nx,ny,nz,nxv,nyv
     1,nxyzv,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 94 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = position z of particle n
c part(4,n) = velocity vx of particle n
c part(5,n) = velocity vy of particle n
c part(6,n) = velocity vz of particle n
c fxyz(1,j,k,l) = x component of force/charge at grid (j,k,l)
c fxyz(2,j,k,l) = y component of force/charge at grid (j,k,l)
c fxyz(3,j,k,l) = z component of force/charge at grid (j,k,l))
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny/nz = system length in x/y/z direction
c nxv = first virtual dimension of field array, must be >= nx+1
c nyv = second virtual dimension of field array, must be >= ny+1
c nxyzv = actual dimension of field array, must be >= nxv*nyv*(nz+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      implicit none
      integer idimp, nop, nx, ny, nz, nxv, nyv, nxyzv, ipbc
      real qbm, dt, ek
      real part, fxyz
      dimension part(idimp,nop), fxyz(3,nxyzv)
c local data
      integer j, nnn, mmm, lll, nop1, nxyv, nn, mm, ll, mp, lp
      real qtm, edgelx, edgely, edgelz, edgerx, edgery, edgerz
      real dxn, dyn, dzn, dxp, dyp, dzp, amx, amy, amz, dx, dy, dz, dx1
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmm = part(2,1)
      lll = part(3,1)
      dxn = part(1,1) - real(nnn)
      dyn = part(2,1) - real(mmm)
      dzn = part(3,1) - real(lll)
      nop1 = nop - 1
      qtm = qbm*dt
      nxyv = nxv*nyv
c set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgelz = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      edgerz = real(nz)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgelz = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
         edgerz = real(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      endif
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      nnn = part(1,j+1)
      mmm = part(2,j+1)
      lll = part(3,j+1)
      dxp = dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j+1) - real(nnn)
      dyn = part(2,j+1) - real(mmm)
      dzn = part(3,j+1) - real(lll)
      amx = 1.0 - dxp
      amy = 1.0 - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzp
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
c find acceleration
      dx = amz*(amx*fxyz(1,mm) + amy*fxyz(1,mm+1) + dyp*fxyz(1,mp)
     1        + dx1*fxyz(1,mp+1))
     2   + dzp*(amx*fxyz(1,ll) + amy*fxyz(1,ll+1) + dyp*fxyz(1,lp)
     3        + dx1*fxyz(1,lp+1))
      dy = amz*(amx*fxyz(2,mm) + amy*fxyz(2,mm+1) + dyp*fxyz(2,mp)
     1        + dx1*fxyz(2,mp+1))
     2   + dzp*(amx*fxyz(2,ll) + amy*fxyz(2,ll+1) + dyp*fxyz(2,lp)
     3        + dx1*fxyz(2,lp+1))
      dz = amz*(amx*fxyz(3,mm) + amy*fxyz(3,mm+1) + dyp*fxyz(3,mp)
     1        + dx1*fxyz(3,mp+1))
     2   + dzp*(amx*fxyz(3,ll) + amy*fxyz(3,ll+1) + dyp*fxyz(3,lp)
     3        + dx1*fxyz(3,lp+1))
c new velocity
      dx = part(4,j) + qtm*dx
      dy = part(5,j) + qtm*dy
      dz = part(6,j) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j))**2 + (dy + part(5,j))**2
     1            + (dz + part(6,j))**2
      part(4,j) = dx
      part(5,j) = dy
      part(6,j) = dz
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
      dz = part(3,j) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j)
            part(6,j) = -part(6,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(4,j) = -part(4,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(5,j) = -part(5,j)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
      part(3,j) = dz
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
      amx = 1.0 - dxn
      amy = 1.0 - dyn
      mm = mm + nn
      dx1 = dxn*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1.0 - dzn
      ll = mm + nxyv
      amy = dxn*amy
      lp = mp + nxyv
c find acceleration
      dx = amz*(amx*fxyz(1,mm) + amy*fxyz(1,mm+1) + dyp*fxyz(1,mp)
     1        + dx1*fxyz(1,mp+1))
     2   + dzn*(amx*fxyz(1,ll) + amy*fxyz(1,ll+1) + dyp*fxyz(1,lp)
     3        + dx1*fxyz(1,lp+1))
      dy = amz*(amx*fxyz(2,mm) + amy*fxyz(2,mm+1) + dyp*fxyz(2,mp)
     1        + dx1*fxyz(2,mp+1))
     2   + dzn*(amx*fxyz(2,ll) + amy*fxyz(2,ll+1) + dyp*fxyz(2,lp)
     3        + dx1*fxyz(2,lp+1))
      dz = amz*(amx*fxyz(3,mm) + amy*fxyz(3,mm+1) + dyp*fxyz(3,mp)
     1        + dx1*fxyz(3,mp+1))
     2   + dzn*(amx*fxyz(3,ll) + amy*fxyz(3,ll+1) + dyp*fxyz(3,lp)
     3        + dx1*fxyz(3,lp+1))
c new velocity
      dx = part(4,nop) + qtm*dx
      dy = part(5,nop) + qtm*dy
      dz = part(6,nop) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,nop))**2 + (dy + part(5,nop))**2
     1            + (dz + part(6,nop))**2
      part(4,nop) = dx
      part(5,nop) = dy
      part(6,nop) = dz
c new position
      dx = part(1,nop) + dx*dt
      dy = part(2,nop) + dy*dt
      dz = part(3,nop) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,nop)
            part(6,nop) = -part(6,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(4,nop) = -part(4,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(5,nop) = -part(5,nop)
         endif
         if (dz.lt.edgelz) dz = dz + edgerz
         if (dz.ge.edgerz) dz = dz - edgerz
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      part(3,nop) = dz
c normalize kinetic energy
   20 ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
