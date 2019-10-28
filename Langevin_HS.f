!*********************************************************************************************************************
!* This program calculates the mean-square displacement for a monodisperse system of hard-spheres in 3D by using
!* Brownian Dynamics simulations
!* 18-Nov-2015
!* The diameter is used as unit length.
!*********************************************************************************************************************
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter(mp=1024) !* Número máximo de partículas
	  parameter(mt=100000) !* Tamaño de la ventana de tiempo
	  parameter(nm=2**8) !* Partición de r para funciones que dependen de r
      dimension x(mp),y(mp),z(mp) !* Posiciones de las partículas
      dimension fx(mp),fy(mp),fz(mp) !* Fuerzas
      dimension t(mt),wt(mt),ft(mt) !* tiempo, desplazamiento cuadrático medio, función de dispersión intermedia (auto)
      dimension r(nm),g(nm),q(nm),sq(nm)!*,h(nm)* !* Partición de r, función de distribución radial, vector de onda, factor de estructura, transformada de Fourier de g - 1
      dimension cfx(mt,mp),cfy(mt,mp),cfz(mt,mp) !* Coordenadas i-ésima partícula al tiempo t
      !* El common es para constantes. De preferencia, variables no físicas. Su scope es todo el código, incluyendo dentro de subrutinas
!* La dinamica molecular vive en el micro canonico y es um problema implementar el termostato. DB vive en el canonico y ese no es problema
!* La temperatura reducida siempre es 1.
      common/sys/rho,phi !* rho^* = rho*sigma^3, fracción de empaquetamiento V_p/V = (pi/6)*rho^*
      common/box/boxl,rc,np !* Longitud de la caja, radio de corte, número de partículas en la simulación
      common/time/deltat !* Paso en el tiempo. Sigue siendo parámetro no físico
	  common/poths/dlr,dla,a2,tem !* Parámetros del mapeo continuo de HS: Lambda repulsiva, lambda atractiva, A del paper, inverso de epsilon estrella (paper, también)
      pi=4.d0*datan(1.d0)
c number of particles
	  np=6**3
c packing fraction
      phi=0.20
      rho=6*phi/pi
	  d=(1./rho)**(1.d0/3.d0) !* distancia promedio entre partículas. Meramente geométrica.
c time step
	   dtt=0.00001 !* Delta t estrella
c number of configurations to thermalize
	   nct=500000 !* Tiempo de espera hasta que termalize
c box length
	   boxl=(dfloat(np)/rho)**(1.d0/3.d0) !* Cast a float. np debe tener una raíz cúbica entera
	   rc=boxl/2.d0 !* Radio de corte
	   dr=rc/nm
       dq=pi/rc
	   print*,'The length of the box is: ',boxl
	   print*,'The mean interparticle distance is: ',d
	   print*,'Cut radius: ',rc

c pair potentials
c effective hard-sphere potential
       dlr=50.d0 !* lambda repulsiva
       dla=49.d0 !* lambda atractiva
       a2=(dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
       tem=1.4737d0 !* 1/e estrella, del paper de Báez

c subroutine for the initial configuration
	   iseed=1234567890 !* Semilla para los aleatorios
c initial configuration for the centers of mass
	   call iniconf(x,y,z,d)
	   open(20,file='iniconfHS.dat',status='unknown')
	   do i=1,np
	      write(20,400)x(i),y(i),z(i) !* El 400 es el formato. Se especifica en la línea 166
	   enddo
	   close(20)
c initial force on the particles
       call force(x,y,z,fx,fy,fz,ener) !* ener es una de las salidas. Está declarada en el implicit

c Energy of the initial configuration
	   print*,'E/N = ',ener/np

       deltat=dtt

c The system will thermalize

c Periodic boundary conditions; pbc > 0
       pbc=1.d0 !* Periodic boundary conditions if this is greater than 0

	  open(51,file='energy.dat',status='unknown')
      do istep=1,nct
         call position(x,y,z,fx,fy,fz,pbc)
         call force(x,y,z,fx,fy,fz,enerpot)
         epotn=enerpot/dfloat(np)
	     if (mod(istep,10000) .eq. 0) print*,istep,epotn,'T'
	     if (mod(istep,100) .eq. 0) write(51,400)dfloat(istep),enerpot/np
      enddo

	  close(51)

      open(60,file='finalcon.dat',status='unknown')
      do i=1,np 
         write(60,400)x(i),y(i),z(i)
      enddo  
      close(60)

c Thermalization ends
	  do i=1,nm
         g(i)=0.d0
         !*h(i)=0.d0
         sq(i)=0.d0
      enddo

      do i=1,mt
         t(i)=0.d0
         wt(i)=0.d0
         ft(i)=0.d0
      enddo

      ncep=10
      ncp=1000000
      nprom=0
      nconf=ncp

      pbc=0.d0

      do i=1,nconf
         call position(x,y,z,fx,fy,fz,pbc)
         call force(x,y,z,fx,fy,fz,enerpot)
         if (mod (i,10000) .eq. 0) print*,i,enerpot/np,'Average'
         if (mod (i,ncep) .eq. 0) then
            nprom=nprom+1
            t(nprom)=deltat*ncep*(nprom-1)
	        do j=1,np
               cfx(nprom,j)=x(j)
               cfy(nprom,j)=y(j)
               cfz(nprom,j)=z(j)
	        enddo
	        call gr(x,y,z,g,dr)
         endif
      enddo

c      goto 101

      open(65,file='gr0p20_BD.dat',status='unknown')
      write(65,400)r(1),g(1)

c      print*,dr,nprom

      do i=2,nm
         r(i)=(i-1)*dr
         q(i)=(i-1)*dq
         dv=4.d0*pi*r(i)**2.*dr
         fnorm=boxl**3./(np**2*nprom*dv)
         graux=g(i)*fnorm
         hraux=graux-1.d0
         g(i)=graux
         !*h(i)=hraux
         write(65,400)r(i),graux,hraux
      enddo
      close(65)

c      call structuref(h,r,q)

c      open(60,file='sq0p005.dat',status='unknown')
c      do i=2,nm
c         sq(i)=1.d0+rho*h(i)
c         write(60,400)q(i),sq(i)
c      enddo
c      close(60)

c101      call difusion(nprom,cfx,cfy,cfz,t,wt,ft)

c      open(80,file='wt0p005.dat',status='unknown')

c      do i=1,(nconf/ncep)-1
c         write(80,400)t(i+1),wt(i),ft(i)
c      enddo
 
c      close(80)

100   format(3f15.7)
400	  format(3f16.8) !* Cada fila del archivo tiene 3 números de dieciséis cifras enteras y ocho decimales
      end
c
c This subroutine calculates the initial configuration in 3D.		 
      subroutine iniconf(xc,yc,zc,d)
      implicit integer*8(i-n),real*8(a-h,o-z)
	  parameter(mp=1024)
	  dimension xc(mp),yc(mp),zc(mp)
	  common/box/boxl,rc,np

	  xc(1)=-(boxl-d)/2.d0
	  yc(1)=-(boxl-d)/2.d0
	  zc(1)=-(boxl-d)/2.d0

      do i=2,np
         xc(i)=xc(i-1)+d
         yc(i)=yc(i-1)
         zc(i)=zc(i-1)
         if (xc(i) .gt. boxl/2.) then
            xc(i)=xc(1)
            yc(i)=yc(i-1)+d
            if (yc(i) .gt. boxl/2.) then
               xc(i)=xc(1)
               yc(i)=yc(1)
               zc(i)=zc(i-1)+d
            endif
         endif
      enddo

	  return
      end
c
      subroutine force(x,y,z,fx,fy,fz,ener)
      implicit integer*8(i-n),real*8 (a-h,o-z)
      parameter(mp=1024)
      dimension x(mp),y(mp),z(mp)
      dimension fx(mp),fy(mp),fz(mp)
	  common/poths/dlr,dla,a2,tem
      common/box/boxl,rc,np
      pi=4.d0*datan(1.d0)
      ener=0.d0
      do i=1,np 
         fx(i)=0.d0 
         fy(i)=0.d0
	     fz(i)=0.d0
      enddo
c pair contribution
      do i=1,np-1 
         do j=i+1,np
            uij=0.d0
            fij=0.d0
            fxij=0.d0
            fyij=0.d0
            fzij=0.d0
            xij=x(i)-x(j)
            yij=y(i)-y(j)
		    zij=z(i)-z(j)
            xij=xij-boxl*dnint(xij/boxl) 
            yij=yij-boxl*dnint(yij/boxl) 
	        zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2) 
            if (rij .lt. rc) then
               call hardsphere(rij,uij,xij,yij,zij,fxij,fyij,fzij)
               ener=ener+uij
               fx(i)=fx(i)+fxij 
               fy(i)=fy(i)+fyij
	           fz(i)=fz(i)+fzij
               fx(j)=fx(j)-fxij 
               fy(j)=fy(j)-fyij
	           fz(j)=fz(j)-fzij
            endif 
         enddo 
      enddo
      return
      end

      subroutine hardsphere(rij,uij,xij,yij,zij,fxij,fyij,fzij)
      implicit integer*8(i-n),real*8 (a-h,o-z)
      parameter(mp=1024)
      common/poths/dlr,dla,a2,tem
      common/box/boxl,rc,np
      if (rij .lt. (dlr/dla)**(1./(dlr-dla))) then
         uij=(a2/tem)*((1./rij)**dlr-(1./rij)**dla)+1./tem
         fij=dlr*(1./rij)**(dlr+1.d0)-dla*(1./rij)**(dla+1.d0)
         fij=(a2/tem)*fij
      else
         uij=0.d0
         fij=0.d0
      endif
      fxij=fij*xij/rij
      fyij=fij*yij/rij
      fzij=fij*zij/rij
      return
      end

c This subroutine calculates the new position of the particles
      subroutine position(x,y,z,fx,fy,fz,pbc)
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter(mp=1024)
      dimension x(mp),y(mp),z(mp)
      dimension fx(mp),fy(mp),fz(mp)
	  common/box/boxl,rc,np
      common/times/deltat
      sigma=dsqrt(2.d0*deltat)
      do i=1,np
         dx=sigma*gasdev() !* gasdev es la distribución gaussina
         dy=sigma*gasdev()
	     dz=sigma*gasdev()
         x(i)=x(i)+fx(i)*deltat+dx
         y(i)=y(i)+fy(i)*deltat+dy
	     z(i)=z(i)+fz(i)*deltat+dz
         if (pbc .gt. 0.d0) then
            x(i)=x(i)-boxl*dnint(x(i)/boxl)
            y(i)=y(i)-boxl*dnint(y(i)/boxl)
            z(i)=z(i)-boxl*dnint(z(i)/boxl)
         endif
      enddo
      return
      end

c This sobroutine computes the mean-square displacement
      subroutine difusion(nprom,cfx,cfy,cfz,t,wt,ft)
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter(mp=1024,mt=100000)
      dimension t(mt) 
      dimension wt(mt),ft(mt),f10t(mt)
      dimension cfx(mt,mp),cfy(mt,mp),cfz(mt,mp)
      common/box/boxl,rc,np
c Mean-square displacement and intermediate scattering function
      dk=6.6d0
      do i=1,nprom-1 
         dif2=0.d0
         dself=0.d0
         do j=1,nprom-i 
            do k=1,np 
               dx=cfx(j+i,k)-cfx(j,k) 
               dy=cfy(j+i,k)-cfy(j,k)
	           dz=cfz(j+i,k)-cfz(j,k)
               dif2=dif2+dx*dx+dy*dy+dz*dz
               aux=dsqrt(dx*dx+dy*dy+dz*dz)
               aux=dsin(dk*aux)/(dk*aux)
               dself=dself+aux
            enddo 
         enddo 
         aux2=(np*(nprom-i)) 
         dif2=dif2/aux2
         dself=dself/aux2
         wt(i)=wt(i)+dif2
         ft(i)=ft(i)+dself
      enddo 
      return 
      end
c
c This subroutine calculates the g(r)
       subroutine gr(x,y,z,g,dr)
       implicit integer*8(i-n),real*8(a-h,o-z)
       parameter(mp=1024)
       parameter(nm=2**8)
       dimension x(mp),y(mp),z(mp),g(nm)
       common/box/boxl,rc,np
       do i=1,np-1
          do j=i+1,np
             xij=x(j)-x(i)
             yij=y(j)-y(i)
             zij=z(j)-z(i)
             xij=xij-boxl*dnint(xij/boxl)
             yij=yij-boxl*dnint(yij/boxl)
             zij=zij-boxl*dnint(zij/boxl)
             rij2=xij*xij+yij*yij+zij*zij
             rij=dsqrt(rij2)
             if (rij .lt. rc) then
c jfix(x) obtiene el entero de x
                nbin=jfix(rij/dr)+1
                if (nbin .le. nm) then
                   g(nbin)=g(nbin)+2.
                endif
             endif
          enddo
       enddo
       
       return
       end

	!***********************************************************************

      function gasdev() 
      implicit integer*8(i-n),real*8(a-h,o-z)
      real*8 grnd
      save iset,gset 
      data iset/0/ 
      if (iset.eq.0) then 
10       v1=2.d0*grnd()-1.d0 
         v2=2.d0*grnd()-1.d0 
         rsq=v1*v1+v2*v2 
         if (rsq.ge.1.d0.or.rsq.eq.0.d0) goto 10 
         fac=dsqrt(-2.d0*dlog(rsq)/rsq) 
         gset=v1*fac 
         gasdev=v2*fac 
         iset=1 
      else 
         gasdev=gset 
         iset=0 
      endif 
      return 
      end 
          

!************************************************************************
      subroutine sgrnd(seed)
!*
      implicit integer(a-z)
!*
!* Period parameters
      parameter(N     =  624)
!*
      dimension mt(0:N-1)
!*                     the array for the state vector
      common /block/mti,mt
      save   /block/
!*
!*      setting initial seeds to mt[N] using
!*      the generator Line 25 of Table 1 in
!*      [KNUTH 1981, The Art of Computer Programming
!*         Vol. 2 (2nd Ed.), pp102]
!*
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
!*
      return
      end

!************************************************************************
      double precision function grnd()
!*
      implicit integer(a-z)
!*
!* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
!*                                    constant vector a
      parameter(UMASK = -2147483648)
!*                                    most significant w-r bits
      parameter(LMASK =  2147483647)
!*                                    least significant r bits
!* Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
!*
      dimension mt(0:N-1)
!*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
!*                     mti==N+1 means mt[N] is not initialized
!*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!*                        mag01(x) = x * MATA for x=0,1
!*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
!*
      if(mti.ge.N) then
!*                       generate N words at one time
        if(mti.eq.N+1) then
!*                            if sgrnd() has not been called,
          call sgrnd(4357)
!*                              a default initial seed is used
        endif
!*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
!*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!*
!ccc ancienne version (produit des 1.0000)
!c      if(y.lt.0) then
!c        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!c      else
!c        grnd=dble(y)/(2.0d0**32-1.0d0)
!c      endif
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32)
      else
        grnd=dble(y)/(2.0d0**32)
      endif
!*
      return
      end


	function ranf(idum)
	implicit integer*8(i-n),real*8(a-h,o-z)
c      integer idum,IA,IM,IQ,IR,MASK
c      real ran0,AM
      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *MASK=123459876)
c      integer k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ranf=AM*idum
      idum=ieor(idum,MASK)
      return
      end
c
      subroutine structuref(h,r,q)
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter(nm=2**8)
      dimension haux(nm,1),snn(nm)
      dimension r(nm),q(nm),h(nm)
      dimension a0(1),a1(1),a2(1)
      common/sys/rho
      common/box/boxl,rc,np
      pi=4.d0*datan(1.d0)
      do i=1,nm
         haux(i,1)=h(i)
      enddo
c      call counti(haux,r,nm,a0,a1,a2)
      call fftm(haux,rc,nm,1)
c      call discount(haux,q,nm,a0,a1,a2)
      do i=1,nm
         h(i)=haux(i,1)
      enddo
      return
      end
c
      subroutine counti(f,r,n,a0,a1,a2)
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter (nm=2**8)
      dimension f(nm,1),r(nm),a0(1),a1(1),a2(1),m(1)
	dimension diam(1)
      common/sys/rho

	  diam(1)=1.d0

      dr=r(2)-r(1)
      do k=1,1
         do i=1,n
            if (r(i).lt.diam(k)) then
               m(k)=i
            endif
         enddo
      enddo

      do k=1,1
      mk=m(k)
      fpd=-274.d0*f(mk+1,k)+600.d0*f(mk+2,k)-600.d0*f(mk+3,k)
      fpd=fpd+400.d0*f(mk+4,k)-150.d0*f(mk+5,k)+24.d0*f(mk+6,k)
      fpd=fpd/(120.d0*dr)
      fppd=225.d0*f(mk+1,k)-770.d0*f(mk+2,k)+1070.d0*f(mk+3,k)
      fppd=fppd-780.d0*f(mk+4,k)+305.d0*f(mk+5,k)-50.d0*f(mk+6,k)
      fppd=fppd/(60.d0*dr**2)
      fd=fpd*(diam(k)-r(mk+1))+f(mk+1,k)
      fpi=274.d0*f(mk,k)-600.d0*f(mk-1,k)+600.d0*f(mk-2,k)
      fpi=fpi-400.d0*f(mk-3,k)+150.d0*f(mk-4,k)-24.d0*f(mk-5,k)
      fpi=fpi/(120.d0*dr)
      fppi=225.d0*f(mk,k)-770.d0*f(mk-1,k)+1070.d0*f(mk-2,k)
      fppi=fppi-780.d0*f(mk-3,k)+305.d0*f(mk-4,k)-50.d0*f(mk-5,k)
      fppi=fppi/(60.d0*dr**2)
      fi=fpi*(diam(k)-r(mk))+f(mk,k)
      delf=fd-fi
      delfp=fpd-fpi
      delfpp=fppd-fppi
      a0(k)=delf-diam(k)*delfp+diam(k)**2*delfpp/2.d0
      a1(k)=delfp-diam(k)*delfpp
      a2(k)=delfpp/2.d0
      do i=1,mk
         f(i,k)=f(i,k)+a0(k)+a1(k)*r(i)+a2(k)*r(i)**2
      enddo
      enddo
      return
      end
c 
      subroutine discount(f,q,n,a0,a1,a2)
      implicit integer*8(i-n),real*8(a-h,o-z)
      parameter (nm=2**8)
      dimension f(nm,1),q(nm),fa(nm,1),a0(1),a1(1),a2(1)
	dimension diam(1)
      common/sys/rho

      pi=4.d0*datan(1.d0)

	diam(1)=1.d0

      do k=1,1
      do i=1,n
         if (q(i).gt.0.d0) then
            part1=a0(k)+2.d0*diam(k)*a1(k)+3.d0*diam(k)**2*a2(k)
            part1=part1-6.d0*a2(k)/q(i)**2
            part1=part1*dsin(q(i)*diam(k))/q(i)**3
            part2=-diam(k)*a0(k)+2.d0*a1(k)/q(i)**2-diam(k)**2*a1(k)
            part2=part2+6.d0*diam(k)*a2(k)/q(i)**2
            part2=(part2-diam(k)**3*a2(k))*dcos(q(i)*diam(k))/q(i)**2
            part3=-2.d0*a1(k)/q(i)**4
            fa(i,k)=4.d0*pi*(part1+part2+part3)
            f(i,k)=f(i,k)-fa(i,k)
         else
            f(i,k)=0.d0
         endif
      enddo
      enddo
      return
      end
c
C **************************************
c
      SUBROUTINE FFTM(FA,RMAX,N,II)
      implicit integer*4(I-N), REAL*8(A-H,O-Z)      
      PARAMETER (nm=2**8)
      DIMENSION FA(NM,1),F(NM)
      DO 10 I=1,1
         DO 20 K=1,N
20          F(K)=FA(K,I)
         CALL FFT(F,RMAX,N,II)
         DO 30 K=1,N
30          FA(K,I)=F(K)
10    CONTINUE
      RETURN
      END
c
C **************************************
c
      SUBROUTINE FFT(F,RMAX,N,II)
      PARAMETER (nm=2**8)
C SI II=1 TF, SI II=2 TFI.
      REAL*8 F,A,RMAX
c	integer*8 I,N,II
      DIMENSION F(NM)
      DO 10 I=1,N
         F(I)=(I-1)*F(I)
10    CONTINUE
      CALL SINFT(F,N,NM)
      IF (II.EQ.1) GOTO 20
      A=N*(1.D0/2.D0/RMAX**3)
      GOTO 30
20    A=4*RMAX**3/(1.D0*N**2)
      F(1)=0.D0
30    DO 40 I=2,N
      F(I)=A*F(I)/(1.D0*(I-1))
40    CONTINUE
      RETURN
      END
c
C *****************************
c
      SUBROUTINE SINFT(Y,N,NMAX)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,Y
c	integer*8 I,J,N,M,NMAX
      DIMENSION Y(NMAX)
      PI=4*DATAN(1.D0)
      THETA=PI/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.D0
      M=N/2
      DO 11 J=1,M
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
      Y1=WI*(Y(J+1)+Y(N-J+1))
      Y2=0.5*(Y(J+1)-Y(N-J+1))
      Y(J+1)=Y1+Y2
      Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL REALFT(Y,M,+1,NMAX)
      SUM=0.0
      Y(1)=0.5*Y(1)
      Y(2)=0.D0
      DO 12 J=1,N-1,2
      SUM=SUM+Y(J)
      Y(J)=Y(J+1)
      Y(J+1)=SUM
12    CONTINUE
      RETURN
      END
c
C ****************************************
c
      SUBROUTINE REALFT(DATA,N,ISIGN,NMAX)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,DATA
c	integer*8 I,N,NMAX,ISIGN
      DIMENSION DATA(2*NMAX)
      PI=4*DATAN(1.D0)
      THETA=PI/DBLE(N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
      C2=-0.5
      CALL FOUR1(DATA,N,+1,NMAX)
      ELSE
      C2=0.5
      THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3
c     DO 11 I=2,N/2+1
      do 11 i=2,n/2
      I1=2*I-1
      I2=I1+1
      I3=N2P3-I2
      I4=I3+1
      WRS=SNGL(WR)
      WIS=SNGL(WI)
      H1R=C1*(DATA(I1)+DATA(I3))
      H1I=C1*(DATA(I2)-DATA(I4))
      H2R=-C2*(DATA(I2)+DATA(I4))
      H2I=C2*(DATA(I1)-DATA(I3))
      DATA(I1)=H1R+WRS*H2R-WIS*H2I
      DATA(I2)=H1I+WRS*H2I+WIS*H2R
      DATA(I3)=H1R-WRS*H2R+WIS*H2I
      DATA(I4)=-H1I+WRS*H2I+WIS*H2R
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
      H1R=DATA(1)
      DATA(1)=H1R+DATA(2)
      DATA(2)=H1R-DATA(2)
      ELSE
      H1R=DATA(1)
      DATA(1)=C1*(H1R+DATA(2))
      DATA(2)=C1*(H1R-DATA(2))
      CALL FOUR1(DATA,N,-1,NMAX)
      ENDIF
      RETURN
      END
c
C ***************************************
c
      SUBROUTINE FOUR1(DATA,NN,ISIGN,NMAX)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,DATA
c	integer*8 I,J,M,N,NN,NMAX,ISIGN
      DIMENSION DATA(2*NMAX)
      PI=4*DATAN(1.D0)
      N=2*NN
      J=1
      DO 11 I=1,N,2
      IF (J.GT.I) THEN
      TEMPR=DATA(J)
      TEMPI=DATA(J+1)
      DATA(J)=DATA(I)
      DATA(J+1)=DATA(I+1)
      DATA(I)=TEMPR
      DATA(I+1)=TEMPI
      ENDIF
      M=N/2
1     IF ((M.GE.2).AND.(J.GT.M)) THEN
      J=J-M
      M=M/2
      GOTO 1
      ENDIF
      J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
      ISTEP=2*MMAX
      THETA=2*PI/(ISIGN*MMAX)
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      DO 13 M=1,MMAX,2
      DO 12 I=M,N,ISTEP
      J=I+MMAX
      TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
      TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
      DATA(J)=DATA(I)-TEMPR
      DATA(J+1)=DATA(I+1)-TEMPI
      DATA(I)=DATA(I)+TEMPR
      DATA(I+1)=DATA(I+1)+TEMPI
12    CONTINUE
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
13    CONTINUE
      MMAX=ISTEP
      GOTO 2
      ENDIF
      RETURN
      END
c
