c Plot contours corresponding to orbits in a specified potential.
c This version plots electric potential and allows for a sign of charge.
c And plots ion phase space contours, and allows delphi!=0.
      integer nv,nx,ncl
      parameter (nv=101,nx=101,ncl=100)
      real phi(nx),f(nx,nv),en(nx,nv),cworka(nv,nx),x(nx),v(nv)
      real vinfin(nx,nv),vmodinfin(nx,nv)
      real vrange,xrange,zclv(ncl)
      real xp(2),yp(2)
      real mass

      cs=0.018
      xrange=4
      phi0=.575
      phi0=.8
      delphi=.2
      vh=1.5
c Electrons
      charge=-1.
      mass=1/1846.
c Ions
!      charge=1.
!      mass=1836
      vrange=(4.+charge)/sqrt(mass)
!      vhole=.065
      phizero=0.-abs(delphi/2.)
      if(vrange.lt.0.1)phizero=.0001

      do i=1,nx
         x(i)=xrange*(-1. + 2.*(i-1.)/(nx-1))
         phi(i)=phifunc(phi0,x(i),delphi)
         do j=1,nv
            v(j)=vrange*(-1. + 2.*(j-1.)/(nv-1))
            en(i,j)=0.5*mass*(v(j))**2+charge*phi(i)
c            vinfin(i,j)=sign(sqrt(2*en(i,j)/mass),v(j))
            vmodinfin(i,j)=(sqrt(2*en(i,j)/mass))
            vinfin(i,j)=vmodinfin(i,j)
            if(v(j).lt.0)vinfin(i,j)=0.
            f(i,j)=exp(-en(i,j)-delphi/2.)
         enddo
      enddo
      
c      write(*,*)x,v
c      write(*,*)f
      call pfset(3)
      call multiframe(3,1,1)
      call pltinit(-xrange,xrange
     $     ,min(phizero,phi0*1.05),max(0.,phi0*1.05))
c      call autoplot(x,phi,nx)
      call polyline(x,phi,nx)
      call charsize(cs,cs)
      call ticset(0.015,0.015,-0.04,-0.02,0,0,0,0)
      call axis()
      call axis2()
      call axlabels('','Potential   !Af!@')
      call polyline([-xrange,xrange],[0.,0.],2)
!      call legendline(.08,.5,258,'Electric Potential')
c Panel 2 electrons
      call pltinit(-xrange,xrange,-vrange,vrange)
      call charsize(cs,cs)
      icl=2
      zclv(1)=1.001
      zclv(2)=.999
!     icsw=1+16+32+64    ! Takes for ever. May be broken.
      icsw=1+16+32
      call accisgradinit(65000,50000,50000,65535,65535,65535)
      call contourl(f,cworka,nx,nx,nv,zclv,icl,x,v,icsw) 
      call charsize(cs,cs)
      icl=-999 
      zclv(1)=16.
      icsw=1
      call contourl(f,cworka,nx,nx,nv,zclv,icl,x,v,icsw) 
      call ticset(0.015,0.015,-0.04,-0.02,0,0,0,0)
      call axis()
      call axis2()
      call charsize(cs,cs)
      call axlabels('','velocity !Bv!de!d!@')
      call legendline(.08,.07,258,
     $     '!Bf!@(!Bv,x!@) contours = energy contours = orbits')
c Annotations
      xp(1)=-3.5
      xp(2)=-3.5
      yp(1)=-vrange
      yp(2)=3.*vrange
      call jdrwstr(wx2nx(xp(2))+.05,wy2ny(.45*vrange),'>',0)
      call jdrwstr(wx2nx(xp(2))+.05,wy2ny(-.45*vrange),'<',0)
      call jdrwstr(wx2nx(xp(2))+.05,wy2ny(vrange*.8),'Passing',1.)
      if(.false.)then
      call dashset(1)
      call polyline(xp,yp,2)
      call jdrwstr(wx2nx(xp(2)),wy2ny(yp(2))-.02,'!Bx!@!d1!d',1.2)
      xp(1)=-xp(1)
      xp(2)=-xp(2)
      call jdrwstr(wx2nx(xp(2)),wy2ny(yp(2))-.02,'!Bx!@!d2!d',-1.2)
      call polyline(xp,yp,2)
      endif
      call jdrwstr(wx2nx(0.),wy2ny(0.),'Trapped',1.)
c Panel 3 ions. First recalculate.
      charge=1.
      mass=1
      vrange=(4.+charge)/sqrt(mass)*.8
      phizero=0.
      if(vrange.lt.0.1)phizero=.0001

      do i=1,nx
         x(i)=xrange*(-1. + 2.*(i-1.)/(nx-1))
         phi(i)=phifunc(phi0,x(i),delphi)
         isigmax=nint(sign(1.,x(i)))
!         write(*,*)'x(i),isigmax',x(i),isigmax
        do j=1,nv
            v(j)=vrange*(-1. + 2.*(j-1.)/(nv-1))
            en(i,j)=0.5*mass*(v(j))**2+charge*phi(i)
            isigmainf=isigmax   ! Inward going or reflected ions
            if(en(i,j).ge.phi0.and.sign(1.,v(j)).eq.isigmax)then
               isigmainf= -isigmainf ! Outward going passing ions.
            endif
            vmodinfin(i,j)=sqrt((2*max(0.,en(i,j)-isigmainf*delphi/2.))
     $           /mass)
            vinfin(i,j)=-isigmainf*vmodinfin(i,j)
            f(i,j)=exp(-(vinfin(i,j)-vh)**2*0.5*mass)
!            if(i.eq.1)write(*,*)v(j),vinfin(i,j),f(i,j),isigmax
!     $           ,isigmainf
         enddo
      enddo
      call pltinit(-xrange,xrange,-vrange,vrange)
      call charsize(cs,cs)
      call accisgradinit(-15000,25000,-200000,65535,75535,65535)
      icl=2
      zclv(1)=0.
      zclv(2)=1.
      icsw=1+16+32
      call contourl(f,cworka,nx,nx,nv,zclv,icl,x,v,icsw) 
      call jdrwstr(wx2nx(-xrange*.9),wy2ny(vrange/2),'>',0)
      call jdrwstr(wx2nx(xrange*.9),wy2ny(vrange/2),'>',0)
      call jdrwstr(wx2nx(-xrange*.9),wy2ny(-vrange/2),'<',0)
      call jdrwstr(wx2nx(xrange*.9),wy2ny(-vrange/2),'<',0)
      call color(15)
      call gradlegend(0.,1.,1.07,0.,1.07,1.,.02,.true.)
      call legendline(1.01,0.5,258,'!Bf!di!d!@')
      icl=-999 
      zclv(1)=17.
      icsw=1
      call contourl(vmodinfin,cworka,nx,nx,nv,zclv,icl,x,v,icsw) 
      
      call ticset(0.015,0.015,-0.04,-0.02,0,0,0,0)
      call axis()
      call axis2()
      call charsize(cs,cs)
      call jdrwstr(wx2nx(xrange/2.),wy2ny(0.),'Reflected',1.)         
      call jdrwstr(wx2nx(-xrange/2.),wy2ny(0.),'Reflected',-1.)         
      call axlabels('position !Bx!@','velocity !Bv!di!d!@')
!      call legendline(.08,.05,258,
!     $     '!Bf!@(!Bv,x!@) contours = energy contours = orbits')



      call pltend()

      end
c************************************************************
      real function phifunc(phi0,x,delphi)
      phifunc=phi0*exp(-x**2) +tanh(x/1.)*delphi/2.
      end
c************************************************************
