module AsymHill
  integer, parameter :: npts=200,nofv=400,nvh=100,ncmax=4,nsmax=20
  real, dimension(npts) :: x,phiofx,denofx,deninteg,dFdelx
  real, dimension(nofv) :: vofv,fofv,Pfofv,Forceofv
  real, dimension(nvh) :: vha,Forcevh,delFdxvh,dphivh,fofvh
  real, dimension(ncmax) :: vs=[1.5,-1.5,0.,0.],vt=1.,dc=[1.,.5,0.,0.],vss
  real, dimension(nvh,nsmax) :: Fcvhns,dFvhns,fvvhns,dpvhns
  real, dimension(nsmax) :: fv0ns,dF0ns,vh0ns
  real :: phimax=.2,phi=.1,xmax=9.8
  character*30 string,argument
  real :: vh=0.,Te=1.,denave,vh0,Fvh0,dFdx0,vhmin=-4.9,vhmax=4.9
  integer :: index,nc=2,ns=1,pfint=0
  logical :: local=.false.,lcd=.true.,ltestnofx=.false.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine parseAsymargs
    ig=1
    do ia=1,iargc()
       call getarg(ia,argument)
       if(argument(1:2).eq.'-p')read(argument(3:),*)phimax
       if(argument(1:2).eq.'-e')read(argument(3:),*)phi
       if(argument(1:2).eq.'-g'.and.ig.le.ncmax)then
          nc=ig
          read(argument(3:),*,end=101)vs(ig),vt(ig),dc(ig)
101       ig=ig+1
       endif
       if(argument(1:2).eq.'-T')read(argument(3:),*)Te
       if(argument(1:2).eq.'-s')read(argument(3:),*)ns
       if(argument(1:2).eq.'-f')ltestnofx=.true.
       if(argument(1:2).eq.'-h')goto 120
       if(argument(1:2).eq.'-v')read(argument(3:),*,end=102)vhmin,vhmax
       if(argument(1:2).eq.'-w')then
          read(argument(3:),*,end=102)pfint
          call pfset(pfint)
       endif
102    continue
    enddo
    return
120 write(*,*)'Usage: AsymHill [-p<phimax> -g<vshift>,<vthermal>,<density> ...'
    write(*,*)' -g...  Enter Maxwellian parameters, starting at first. Currently:'
    write(*,'(''   '',$)')
    write(*,11)(i,' [',vs(i),vt(i),dc(i),']',i=1,nc)
    write(*,'(a,f6.3)')'  -p...  Set peak potential             [',phimax
    write(*,'(a,f6.3)')'  -e...  Set test potential             [',phi
    write(*,'(a,f6.3)')'  -T...  Set electron temperature       [',Te
    write(*,'(a,i4)'  )'  -s...  Set number of shifts in scan   [',ns
    write(*,'(a,l4)'  )'  -f...  Display f(v) etc               [',ltestnofx
    write(*,'(a,l4)'  )'  -w...  Postscript write (0,3,-3)      [',pfint
    write(*,'(a,2f6.2)')'  -v...  Set vh range                   [',vhmin,vhmax
    call exit
11  format(i2,a,3f6.3,a,i2,a,3f6.3,a,i2,a,3f6.3,a)
  end subroutine parseAsymargs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testfvhill(delphi)
  call multiframe(4,1,1)
  call dcharsize(.02,.02)
  do isigma=-1,1,2
     call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phi,isigma,local,&
          nofv,vofv,fofv,Pfofv)
     call autoplot(vofv,fofv,nofv)
     call axis2
     call axlabels('','f(v)')
     call iwrite(isigma,iwidth,string)
     call legendline(.1,.8,258,'!As!@!dx!d='//string(1:iwidth))
     call fwrite(phi,iwidth,2,string)
     call legendline(.1,.6,258,'!Af!@='//string(1:iwidth))
     call fwrite(phimax,iwidth,2,string)
     call legendline(.1,.4,258,'!Ay!@='//string(1:iwidth))
     call autoplot(vofv,Pfofv,nofv)
     call axis2
     if(isigma.eq.-1)call axlabels('','P(v)=!AJ!@f(v)dv')
     if(isigma.eq.+1)call axlabels('v-v!dh!d','P(v)=!AJ!@f(v)dv')
  enddo
  call multiframe(0,0,0)
  call pltend
end subroutine testfvhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddenofx(delphi)
! Integrate over the potential hill to find the density, force, and
! delF/delx (arrays of x), for the current vh and other parameters. 
  deninteg=0.
  dFdelx=0.
  do i=1,npts
     x(i)=-xmax+2.*xmax*(i-1.)/(npts-1.)
     isigma=int(sign(1.,x(i)))
     dph=isigma*delphi/2.
     xcor=x(i)/4.
!     xcor=x(i)/4./sqrt(1.-dph/phimax)  !Symmetrizes curvature at origin.
     phiofx(i)=(phimax-dph)/cosh(xcor)**4+dph
     if(phiofx(i).gt.phimax)then
        write(*,*)'phi>phimax',phiofx(i),phimax,delphi,i,dph,xcor
        stop
     endif
     call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiofx(i),isigma,local, &
          nofv,vofv,fofv,Pfofv)
     denofx(i)=Pfofv(nofv)
     if(i.gt.1)deninteg(i)=deninteg(i-1)+&      !\int n dphi
          (denofx(i)+denofx(i-1))*(phiofx(i)-phiofx(i-1))*.5
     if(i.gt.1)dFdelx(i)=dFdelx(i-1)-&          !-\int (dn/dx) dphi
          (denofx(i)-denofx(i-1))*(phiofx(i)-phiofx(i-1)) &
          /(x(i)-x(i-1))
  enddo
  if(lcd)write(*,'(a,f8.5,a,f9.5,a,f9.5,a,f9.5)')'vh=',vh,' delphi=',delphi,&
       ' Ion force=',-deninteg(npts),' delF/delx=',dFdelx(npts)
end subroutine finddenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotdenofx
  call multiframe(3,1,0)
  call autoplot(x,phiofx,npts)
  call axis2
  call polyline([x(1),x(npts)],[0.,0.],2)
  call fwrite(vh,iwidth,2,string)
  call legendline(.1,1.1,258,'v!dh!d='//string(1:iwidth))
  call axlabels('','!Af!@(x)')
  call autoplot(x,denofx,npts)
  call axis2
  call polyline([x(1),x(npts)],[denave,denave],2)
  call axlabels('','n!di!d(x)')
  call winset(.true.)
  call color(4)
!  call polyline(x,denave*exp(phiofx/Te),npts)
  call polymark([x(1),x(npts)], &
       denave*[exp(phiofx(1)/Te),exp(phiofx(npts)/Te)],2,1)
  call color(15)
  call autoplot(x,deninteg,npts)
  call axis2
  call axlabels('x','!AJ!@n!di!dd!Af!@')
  call polyline([x(1),x(npts)],[0.,0.],2)
  call pltend
  call multiframe(0,0,0)
end subroutine plotdenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotex
  call multiframe(4,1,1)
  call minmax(phiofx,npts,pmin,pmax)

  call pltinit(-10.,10.,pmin*2.,pmax*1.4)
  call legendline(.95,.1,258,'(a)')
  call axis
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(1)
  call axlabels('','Ion Energy');call axptset(0.,0.);call ticrev()
  call polyline([-9.,-9.],[pmin,pmax],2)
  call polyline([-9.2,-8.8],[pmax,pmax],2)
  call jdrwstr(wx2nx(-9.),wy2ny(pmin)+.05,'Reflected!A_!@',1.05)
  call polyline([9.,9.],[phiofx(npts),pmax],2)
  call polyline([9.2,8.8],[pmax,pmax],2)
  call jdrwstr(wx2nx(9.),wy2ny(phiofx(npts))+.04,'!A^!@Reflected',-1.05)
  call polyline([0.,0.],[pmax,pmax*1.4],2)
  call jdrwstr(wx2nx(.2),wy2ny(pmax)+.02,'Passing!A^_!@',1.)
  call color(15)
  call polyline(x,phiofx,npts)
  call polyline([x(1),x(npts)],[0.,0.],2)
  call axlabels('','!Af!@(x)')
  
  call pltinit(-10.,10.,-pmax*1.1,pmax)
  call legendline(.95,.1,258,'(b)')
  call axis
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(4)
  call axlabels('','Electron Energy');call axptset(0.,0.);call ticrev()
  call polyline([-9.,-9.],[-pmin,pmax],2)
  call jdrwstr(wx2nx(-9.),wy2ny(pmin)+.05,'Passing!A^_!@',1.05)
  call polyline([0.,0.],[-pmax,-phiofx(npts)],2)
  call jdrwstr(wx2nx(0.),wy2ny(-pmax+.04),'Trapped!A^!@ ',-1.0)
  call drcstr(' !A_!@')
  call polyline([-.2,.2],[-phiofx(npts),-phiofx(npts)],2)
  call polyline([9.,9.],[-phiofx(npts),-pmin],2)
  call jdrwstr(wx2nx(9.),wy2ny(-phiofx(npts))+.01,'!A^!@Reflected',-1.05)
  call polyline([9.2,8.8],[-pmin,-pmin],2)
  call color(15)
  call polyline(x,-phiofx,npts)
  call polyline([x(1),x(npts)],[0.,0.],2)
  call axlabels('','-!Af!@(x)')

  call autoinit(x,denofx,npts)
  call axis
  call legendline(.95,.1,258,'(c)')
  call color(1)
  call polyline(x,denofx,npts)
  call axlabels('','n!di!d(x)')
  call color(15)
  call polyline([x(1),x(npts)],[denave,denave],2)
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(4)
  call axlabels('','n!de!d(!A+;!@) !A1!@');call axptset(0.,0.);call ticrev()
  call polymark([x(1),x(npts)], &
       denave*[exp(phiofx(1)/Te),exp(phiofx(npts)/Te)],2,1)
  call color(15)
!  call polyline(x,denave*exp(phiofx/Te),npts)

  call autoplot(x,deninteg,npts)
  call legendline(.95,.1,258,'(d)')    
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(1)
  call axptset(0.,0.);call ticrev()
  call axlabels('x','!AJ!@n!di!dd!Af!@')
  call polyline([x(1),x(npts)],[0.,0.],2)
  call pltend
  call multiframe(0,0,0)
end subroutine plotex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddelphi(delphi)
! Iterate to find delphi consistent with specified ion distributions.
! But prevent |delphi| from exceeding 2*phimax because that's improper
  real, dimension(2) :: deninf
  if(lcd)write(*,*)'j   delphi      dpp   delphi-dpp   n(-)      n(+)'
  dpp=delphi
  do j=1,10
     do i=1,2
        isigma=-3+2*i
        phiinf=isigma*delphi/2.
        call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiinf,isigma,local,&
             nofv,vofv,fofv,Pfofv)
        deninf(i)=Pfofv(nofv)
     enddo
!     delphi=delphi+(1-.5/j)*(alog(deninf(2)/deninf(1))-delphi/Te)
     delphi=delphi+(1-.5/j)*(alog(deninf(2)/deninf(1))*min(1.,Te)-delphi/max(1.,Te))
     if(lcd)write(*,'(i2,5f10.5)')j,delphi,dpp,delphi-dpp,deninf
     if(abs(delphi).gt.2.*phimax)then  !Prevent improper delphi
        write(*,*)'WARNING: Improper delphi suppressed',delphi,phimax
        delphi=sign(2.*phimax,delphi)
        exit
     endif
     if(abs(delphi-dpp).lt.1.e-4)exit
     dpp=delphi
  enddo
  denave=deninf(1)*exp(0.5*delphi/Te)
!delphi=0.
end subroutine finddelphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scanvh(index,delphi)
! Scan vhmin to vhmax to construct arrays of F and delF/delx
! Return the highest index of vh at which F crosses from + to -.
  Fprior=0.
  index=0
  fofvh=0.
  do i=1,nvh
     vh=vhmin+(i-1.)*(vhmax-vhmin)/(nvh-1.)
     do k=1,nc             ! Store the finfty of vh.
        vt2x2=vt(k)**2*2.
        gcoef=1./sqrt(vt2x2*3.1415926)
        fofvh(i)=fofvh(i)+dc(k)*gcoef*exp(-(vh-vs(k))**2/vt2x2)
     enddo
     vha(i)=vh
     lcd=.false.
     call finddelphi(delphi)
     dphivh(i)=delphi
     call finddenofx(delphi)
     Force=-deninteg(npts)+denave*Te*(exp(delphi/(2.*Te))-exp(-delphi/(2.*Te)))
     if(Fprior.gt.0.and.Force.le.0.)index=i
     Forcevh(i)=Force
     delFdxvh(i)=dFdelx(npts)
     Fprior=Force
  enddo
end subroutine scanvh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initvofv
! Initialize velocity range using current vh as maximal.
  vmax=4.4+max(abs(maxval(vs)-vh),abs(minval(vs)-vh))
  do i=1,nofv
     vofv(i)=-vmax+2.*vmax*(i-1.)/(nofv-1.)
  enddo
end subroutine initvofv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timenofx
! This routine takes 8 seconds for 1M densities, when npts=400.
  if(nc.eq.2)dc(2)=.5
  if(nc.eq.2)vs(nc)=-vs(1)
  ncycle=2500
  call initvofv
  write(*,*)'Timing routine starting'
  do j=1,ncycle
  do i=1,npts
     x(i)=-xmax+2.*xmax*(i-1.)/(npts-1.)
     phiofx(i)=phimax/cosh(x(i)/4.)**4
     isigma=int(sign(1.,x(i)))
     call fvhill(nc,dc,vs,vt-vh,phimax,0.,phiofx(i),isigma,local, &
          nofv,vofv,fofv,Pfofv)
     denofx(i)=Pfofv(nofv)
     if(.not.abs(denofx(i)).ge.0)write(*,*)i,Pfofv(nofv)
  enddo
  enddo
  write(*,'(a,i5,a,i5,a,i8,a)')'Found denofx for',npts,' points a total of '&
       ,ncycle,' times: ',npts*ncycle,' densities'
end subroutine timenofx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forcedivplot
  call multiframe(2,1,0)
  call dcharsize(.02,.02)
  call autoplot(vha,fofvh,nvh)
  call axlabels('','f!d!A;!@!d(v!d!A;!@!d)')
  call axis2
  if(dFdx0.lt.0)call polymark(vh0,Fvh0,1,10)
   ! Electron force
  call autoplot(vha,denave*Te*(exp(dphivh/(2.*Te))-exp(-dphivh/(2.*Te))),nvh)
  call axlabels('v!d!A;!@!d or v!dh!d','Forces(v!dh!d)')
  call axis2
  call legendline(.05,.8,0,'F!de!d')
  call winset(.true.)
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
  if(index.ne.0)call polyline([vh0,vh0],[-10.,10.],2)
  call color(1)
  call polyline(vha,Forcevh,nvh)      ! Total force
  call legendline(.05,.7,0,'F!di!d+F!de!d')
  call color(2)
  call polyline(vha,-denave*dphivh/Te+Forcevh,nvh)   ! Ion force
  call legendline(.05,.9,0,'F!di!d')
  call pltend
  call multiframe(0,0,0)
end subroutine forcedivplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forceplotns
  call multiframe(3,1,1)
  call dcharsize(.02,.02)
  call minmax(fvvhns(1,1),nvh,symin,symax)
  call pltinit(vhmin,vhmax,symin,symax*1.1)
  call axis
  call axis2
  call axlabels('','f!d!A;!@!d(v!d!A;!@!d)')
  do i=1,ns
     call color(i)
     call polyline(vha,fvvhns(1,i),nvh)
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),fv0ns(i),1,10)
  enddo
  call color(15)
  call minmax(Fcvhns(1,1),nvh,symin,symax)
  call pltinit(vhmin,vhmax,symin*1.1,symax*1.1)
  call axis
  call axis2
  call axlabels('','F(v!dh!d)')  
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
!  call legendline(.05,.7,0,'F!di!d+F!de!d')
  do i=1,ns
     call color(i)
     call polyline(vha,Fcvhns(1,i),nvh)      ! Total force
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),0,1,10)
  enddo
  call color(15)
  call minmax(dFvhns(1,1),nvh,symin,symax)
  call pltinit(vhmin,vhmax,symin*1.1,symax*1.1)
  call axis
  call axis2
  call axlabels('v!dh!d  or  v!d!A;!@!d','!Ad!@F/!Ad!@x(v!dh!d)')  
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
  do i=1,ns
     call color(i)
     call polyline(vha,dFvhns(1,i),nvh)      ! Total force
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),dF0ns(i),1,10)
  enddo
  call pltend
  call multiframe(0,0,0)  
end subroutine forceplotns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module AsymHill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate velocity distributions and densities for shifted Maxwellian
! distant distributions on a single-humped repelling potential hill of
! maximum potential phimax, at potential phi, on side isigma of the hill.
! Also return \int fdv in Pfofv. 
! On entry:
! nc is the number of Maxwellian components 
! dc(nc), vs(nc), vt(nc) is each component's density, vshift, sqrt(T/M)
! phimax is hill's peak potential relative to mean distant potential
! delphi is potential difference phi(+inf)-phi(-inf) across the hill. 
! phi is the potential at which to find the distribution.
! isigma is the side of the hill on which this potential lies.
! vofv(nofv) is the monotonic array of velocities to use.
! All velocities are relative to the hill's rest frame.
! On exit:
! fofv(nofv) is the total velocity distribution at vofv
! Pfofv(nofv) is the integral of fofv dv. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fvhill(nc,dc,vs,vt,phimax,delphi,phi,isigma,local,&
     nofv,vofv,fofv,Pfofv)
  real :: dc(nc),vs(nc),vt(nc),phimax,phi,delphi
  integer :: isigma, nofv
  real, dimension(nofv) :: vofv, fofv, Pfofv
  logical :: local
  
  if(isigma.ne.1.and.isigma.ne.-1)stop 'isigma must be +-1'
! Decide the meaning of f(v) when delphi !=0.  
  if(local)then ! f(v) is relative to different phi_infty
     phix2=2.*(phi    -isigma*delphi/2.)  ! 2*(phi   -phiinf)
     phimx2=2.*(phimax-isigma*delphi/2.)  ! 2*(phimax-phiinf)
     delphix2=2.*isigma*delphi  
     if(phix2.lt.0)then
        write(*,*)'Impossible phi (<phiinf)',phi,isigma,delphi
        stop
     endif
  else           ! f(v)= f( sigma*sqrt(v^2+2*phi) )
     phix2=2.*phi
     phimx2=2.*phimax
     delphix2=0.
  endif
  if(phimax.lt.phi)then
     write(*,*)'Normally phi must be < phimax',phi,phimax
     vthresh=0.
  elseif(phimx2.lt.0)then
     write(*,*)'No Potential peak this side',isigma,phimx2
     vthresh=0.
  else
     vthresh=isigma*sqrt(phimx2-phix2)
  endif
  Pfofv=0.
  fofv=0.
  do j=1,nc
     vt2x2=vt(j)**2*2.
     gcoef=1./sqrt(vt2x2*3.1415926)
     vdiffm=vofv(1)
     fim=0.
     pam=0.
     pa=0.
     do i=1,nofv
        vsign=sign(1.,vofv(i))
        vdiff=vofv(i)-vthresh
        if(int(vsign).ne.isigma)then              ! Moving inward
           vinf=vsign*sqrt(max(0.,phix2+vofv(i)**2))
        elseif((phix2+vofv(i)**2).gt.phimx2)then  ! Passing
           vinf=vsign*sqrt(max(0.,phix2+delphix2+vofv(i)**2))
        else                                      ! Reflected
           vinf=-vsign*sqrt(max(0.,phix2+vofv(i)**2))
        endif
        fi=gcoef*exp(-(vinf-vs(j))**2/vt2x2)   ! f(v)=finf(vinf)
        fofv(i)=fofv(i)+dc(j)*fi
        if(vdiffm*vdiff.le.0.)then
! Crossing threshold at which the discontinuity in f occurs, which is
! where v=isigma*sqrt(phimx2-phix2), i.e. vdiff=0; to correctly integrate:
           pa=pam+ dc(j)*(abs(vdiffm)*fim+abs(vdiff)*fi) &
                /(abs(vdiffm)+abs(vdiff))*(vofv(i)-vofv(i-1))
        else
           if(i.gt.1)pa=pam+dc(j)*0.5*(fi+fim)*(vofv(i)-vofv(i-1))
        endif
        Pfofv(i)=Pfofv(i)+pa
        pam=pa
        vdiffm=vdiff
        fim=fi
     enddo
  enddo
end subroutine fvhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Asym
use AsymHill

call parseAsymargs
call initvofv

if(ltestnofx)then
   call finddelphi(delphi)
   call testfvhill(delphi)
   write(*,'(a,9f6.3)')'testnofx: vshift,vtherm,dens',(vs(i),vt(i),dc(i),i=1,nc)
   call finddenofx(delphi)
!   call plotdenofx
   call plotex
endif

vss=vs
vsfac=1.
if(ns.ge.3)vsfac=2.
do i=1,ns
   vs=vss*vsfac*(i-min(1,ns-1))/(ns-min(1,ns-1))
   index=2
   call scanvh(index,delphi)
   Fcvhns(:,i)=Forcevh
   dFvhns(:,i)=delFdxvh
   fvvhns(:,i)=fofvh
   dpvhns(:,i)=dphivh
   dF0ns(i)=1. 
   if(i.eq.1)write(*,*)&
        'At   vh(i-1)  vh(i)  F(i-1)  F(i)   delF/dx,   F=0  and dF/dv_h<0'
   if(index.ne.0)then
      vh0=(vha(index-1)*abs(Forcevh(index))+vha(index)*abs(Forcevh(index-1)))/ &
           (abs(Forcevh(index-1))+abs(Forcevh(index)))
      fvh0=(fofvh(index-1)*abs(Forcevh(index)) &
           +fofvh(index)*abs(Forcevh(index-1)))/ &
           (abs(Forcevh(index-1))+abs(Forcevh(index)))
      dFdx0=(delFdxvh(index-1)*abs(Forcevh(index)) &
           +delFdxvh(index)*abs(Forcevh(index-1)))/ &
           (abs(Forcevh(index-1))+abs(Forcevh(index)))
      write(*,'(i4,5f8.4)')index,vha(index-1),vha(index+1),&
           Forcevh(index-1),Forcevh(index),dFdx0
      vh=vh0
      dF0ns(i)=dFdx0
      fv0ns(i)=fvh0
      vh0ns(i)=vh0
      vh=vh0
      if(ns.lt.3)then
!         write(*,*)'Calling finddenofx, plotdenofx',i,ns
         call finddelphi(delphi)
         call finddenofx(delphi)
         call plotdenofx
      endif
   endif
   if(ns.le.3)call forcedivplot
enddo
call forceplotns
!call timenofx ! For timing comment the above.
end program Asym
