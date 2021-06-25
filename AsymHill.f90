module AsymHill
  integer, parameter :: npts=200,nofv=200,nvh=100,ncmax=4,nsmax=20
  real, dimension(npts) :: x,phiofx,denofx,deninteg,dFdelx
  real, dimension(nofv) :: vofv,fofv,Pfofv,Forceofv
  real, dimension(nvh) :: vha,Forcevh,delFdxvh,dphivh,fofvh
  real, dimension(ncmax) :: vs=[1.5,-1.5,0.,0.],vt=1.,dc=[1.,.5,0.,0.],vss
  real, dimension(nvh,nsmax) :: Fcvhns,dFvhns,fvvhns,dpvhns
  real, dimension(nsmax) :: fv0ns,dF0ns,vh0ns,delphins
  real :: delphi,phimax=.2,phi=.1,xmax=12.,delphi0
  character*30 string,argument
  real :: vh=0.,Te=1.,denave,vh0,fvh0,dFdx0,vhmin=-4.,vhmax=4.,vhn,vhx,vh0r
  integer :: index,nc=2,ns=1,pfint=0,isigma
  logical :: local=.false.,lcd=.true.,ltestnofx=.false.,ldenion=.false.
  logical :: lrefinethreshold=.true.,lfdconv,lexplain=.false.
  logical :: ldenalt=.true.
! BKGint arrays etc. Reminder u is v/sqrt(2). 
  integer, parameter :: nphi=200
  real, dimension(-nphi:nphi) :: phiofi,xofphi,edenofphi,Vminus,x2ofphi
  real, dimension(0:nphi) :: edentrap,etilden,fofphi,u0ofphi,usofphi,edenuntrap

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
       if(argument(1:2).eq.'-x')lexplain=.not.lexplain
       if(argument(1:2).eq.'-h')goto 120
       if(argument(1:2).eq.'-v')read(argument(3:),*,end=102)vhmin,vhmax
       if(argument(1:2).eq.'-w')then
          read(argument(3:),*,end=102)pfint
          call pfset(pfint)
       endif
102    continue
    enddo
    vhx=vhmax
    vhn=vhmin
    if(ns.lt.3)lrefinethreshold=.false.
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
    write(*,'(a,l4)'  )'  -x...  Plot f(v) explanation          [',lexplain
    write(*,'(a,l4)'  )'  -w...  Postscript write (0,3,-3)      [',pfint
    write(*,'(a,2f6.2)')'  -v...  Set vh range                   [',vhmin,vhmax
    call exit
11  format(i2,a,3f6.3,a,i2,a,3f6.3,a,i2,a,3f6.3,a)
  end subroutine parseAsymargs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddenofx
! Integrate over the potential hill to find the density, force, and
! delF/delx (arrays of x), for the current vh and other parameters. 
  deninteg=0.
  dFdelx=0.
  dx=2.*xmax/(npts-1.)
  do i=1,npts
!     x(i)=-xmax+2.*xmax*(i-1.)/(npts-1.)
     x(i)=-xmax+(i-1.)*dx
     isigma=int(sign(1.,x(i)))
     dph=isigma*delphi/2.
     xcor=x(i)/4.
!     xcor=x(i)/4./sqrt(1.-dph/phimax)  !Symmetrizes curvature at origin.
     phiofx(i)=(phimax-dph)/cosh(xcor)**4+dph
     if(ldenalt)then   ! This density must be consistent with delphi.
        call denhill(nc,dc,vs-vh,vt,vh,phimax,phiofx(i),isigma,nofv,denofx(i))
     else
        call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiofx(i),isigma,local, &
             nofv,vofv,fofv,Pfofv)
        denofx(i)=Pfofv(nofv)
     endif
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
subroutine finddelphi
! Iterate to find delphi consistent with specified ion distributions.
  real, dimension(2) :: deninf
  integer, parameter :: niter=10
  real, dimension(niter) :: residit
  if(lcd)write(*,*)'j   dpp       delphi      delphi-dpp    n(-)        n(+)',&
       '   resid'
  lfdconv=.false.
  dpp=delphi
  ddp=0.
  residp=0.
  do j=1,niter
     do i=1,2
        isigma=-3+2*i
        phiinf=isigma*delphi/2.
        if(ldenalt)then ! Alternative ni calculation.
           call denhill(nc,dc,vs-vh,vt,vh,phimax,phiinf,isigma,nofv,deninf(i))
        else
           call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiinf,isigma,local,&
                nofv,vofv,fofv,Pfofv)
!           write(*,'(a,i4,5f10.6)')'Findelphi',j,phiinf,deninf(i),Pfofv(nofv)
           deninf(i)=Pfofv(nofv)
        endif
     enddo
     resid=(alog(deninf(2)/deninf(1))-delphi/Te)
     residit(j)=resid
     if(j.eq.1)then
        delphi=delphi+.5*resid*min(1.,Te)
     else
        if((resid-residp).eq.0)then
           write(*,'(a,i3,5e12.4)')'resid unchanged',j,resid,delphi,dpp,ddp
           delphi=delphi+.5*resid*min(1.,Te)
        else
           ddp=resid*(delphi-dpp)/(resid-residp)
        endif
        dpp=delphi
!        delphi=delphi-ddp
        delphi=delphi-sign(min(abs(ddp),phimax/2.),ddp)
     endif
     if(lcd)write(*,'(i2,6f11.7)')j,dpp,delphi,delphi-dpp,deninf,resid
     if(abs(resid).lt.0.5e-6)then
        lfdconv=.true.
        goto 2
     endif
     residp=resid
  enddo
  write(*,'(a,7f11.7)')'finddelphi unconverged',resid,residp,delphi,dpp,ddp
!  write(*,'(10f8.4)')residit
2 continue
  denave=deninf(1)*exp(0.5*delphi/Te)
!delphi=0.
end subroutine finddelphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scanvh(index)
! Scan vhmin to vhmax to construct arrays of F and delF/delx
! Return the highest index of vh at which F crosses from + to -.
! Then refine the search in the vicinity of F=0 to find vh0.
  index=0
  fofvh=0.
  vssize=maxval(vt)
  vhn=vhmin
  vhx=vhmax
  vhm=(vhmax+vhmin)/2.
  vhd=(vhmax-vhmin)/2.
  vht=minval(vt)
  vht=vhd/3.
  vh1=0
  vh2=0
  do j=1,3            ! Iterate for refinement of vh0
     Force=0.
     Fmax=0.
     vh=vhn
     indexref=0
     do i=1,nvh          ! Scan over the current range
        vhprior=vh
        delphiprior=delphi
        Fprior=Force
        vh=vhn+(i-1.)*(vhx-vhn)/(nvh-1.)
        lcd=.false.
        call finddelphi
        call finddenofx
        Force=-deninteg(npts) &
             +denave*Te*(exp(delphi/(2.*Te))-exp(-delphi/(2.*Te)))
!        if(lfdconv &
        if(vhd-abs(vh-vhm).gt.vht  &    ! Not too near the ends
             .and.Fprior.gt.0.and.Force.lt.0..and.indexref.eq.0)then
           indexref=i
           vh0=(vhprior*abs(Force)+vh*abs(Fprior))/ &
                (abs(Fprior)+abs(Force))
           delphi0=(delphiprior*abs(Force)+delphi*abs(Fprior))/ &
                (abs(Fprior)+abs(Force))
           delphi0=delphi
           vh1=vhprior
           vh2=vh
        endif
        if(j.eq.1)then    ! Store coarse scan.
           vha(i)=vh
           dphivh(i)=delphi
           Forcevh(i)=Force
           delFdxvh(i)=dFdelx(npts)
           if(dFdelx(npts).le.0.)then
              index=indexref
              fvh0=(fofvh(index-1)*abs(Forcevh(index)) &
                   +fofvh(index)*abs(Forcevh(index-1)))/ &
                   (abs(Forcevh(index-1))+abs(Forcevh(index)))
              dFdx0=(delFdxvh(index-1)*abs(Forcevh(index)) &
                   +delFdxvh(index)*abs(Forcevh(index-1)))/ &
                   (abs(Forcevh(index-1))+abs(Forcevh(index)))
           endif
           do k=1,nc             ! Store the finfty of vh.
              vt2x2=vt(k)**2*2.
              gcoef=1./sqrt(vt2x2*3.1415926)
              fofvh(i)=fofvh(i)+dc(k)*gcoef*exp(-(vh-vs(k))**2/vt2x2)
           enddo
        endif
     enddo
!     write(*,*)'vh0 etc',vh0,Force,delphi0,deninteg(npts)
     if(indexref.gt.1)then  ! Refine the equilibrium vh0.
        vhn=vh1
        vhx=vh2
!        write(*,*)'vhn,vhx',vhn,vhx
     else
        exit ! from refinement j-loop.
     endif
  enddo
!  write(*,*)vh0,index,vha(index-1),vha(index),Forcevh(index-1),Forcevh(index)
end subroutine scanvh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scanspacing
! Iterate scanvh calls over different spacings  
! If ns .eq. 1 or 2 report vh and do plots
  vss=vs
  vsfac=1.
  if(ns.ge.3)vsfac=2.
  do i=1,ns
     vsfaci=vsfac*(i-min(1.,ns-1.))/(ns-min(1.,ns-1.))
     vs=vss*vsfaci
     call scanvh(index)
     Fcvhns(:,i)=Forcevh
     dFvhns(:,i)=delFdxvh
     fvvhns(:,i)=fofvh
     dpvhns(:,i)=dphivh
     dF0ns(i)=1. 
     if(i.eq.1)write(*,*)&
          ' At   vh(i-1)  vh(i)  F(i-1)  F(i)   delF/dx,  vh0',&
          '   vshift  Delphi/psi^2'
     if(index.ne.0)then
! vh0 is instead refined in scanvh
!     vh0=(vha(index-1)*abs(Forcevh(index))+vha(index)*abs(Forcevh(index-1)))/ &
!           (abs(Forcevh(index-1))+abs(Forcevh(index)))
!        dphivh0=(dphivh(index-1)*abs(Forcevh(index)) &
!             +dphivh(index)*abs(Forcevh(index-1)))/ &
!             (abs(Forcevh(index-1))+abs(Forcevh(index)))
        dF0ns(i)=dFdx0
        fv0ns(i)=fvh0
        vh0ns(i)=vh0
        vh=vh0
        call finddelphi
        write(*,'(i4,7f8.4,f9.5)')index,vha(index-1),vha(index),&
             Forcevh(index-1)/phimax**2,Forcevh(index)/phimax**2, &
             dFdx0/phimax**2,vh0,vs(1),delphi0/phimax**2
        delphins(i)=delphi
        if(ns.lt.3)then
           lcd=.true.
           call finddenofx
           call plotdenofx
           call findtrappedf
           call otherxofphi
        endif
     else
        write(*,'(a,a,f8.4)')'     No stable equilibrium found for vs(1)=    ',&
             '     ',vs(1)
     endif
     if(ns.le.3)call forcedivplot
  enddo
  if(lrefinethreshold)then  ! Refinement of threshold: 
! Find vs1fac and vs2fac bracketting threshold
  vhminst=vhmin
  vhmaxst=vhmax
  if(ns.gt.2)then
     do i=1,ns
        if(dF0ns(i).ge.0)then
           vs1fac=vsfac*(i-min(1,ns-1))/(ns-min(1,ns-1))
        elseif(vs1fac.ne.0.)then
           vs2fac=vsfac*(i-min(1,ns-1))/(ns-min(1,ns-1))
           goto 2
        else
           write(*,*)'ERROR: dF0ns negative before positive',dF0ns(i)
        endif
     enddo
     write(*,*)'Failed to bracket vsfacs:',vs1fac,vs2fac
     return
2    continue
     write(*,*)'Refining vshift threshold'
     vh0=vh0ns(i)
     vs1fac=0.
     dvh=(vhmax-vhmin)/5.
     vhmin=vh0-dvh
     vhmax=vh0+dvh
     do i=1,ns
        call bisectspacing(vs1fac,vs2fac)
     enddo
     vh0r=vh0
     dF0ns(ns+1)=dFdx0
     vh0ns(ns+1)=vh0r
     vhmin=vhminst
     vhmax=vhmaxst
     call scanvh(index) ! Needed to reset various arrays to vhmin/max
     Fcvhns(:,ns+1)=Forcevh     ! Set the profiles etc. for refined soln.
     dFvhns(:,ns+1)=delFdxvh
     fvvhns(:,ns+1)=fofvh
     dpvhns(:,ns+1)=dphivh
     vh0ns(ns+1)=vh0r
     fvh0=(fofvh(index-1)*abs(Forcevh(index)) &
          +fofvh(index)*abs(Forcevh(index-1)))/ &
          (abs(Forcevh(index-1))+abs(Forcevh(index)))
     fv0ns(ns+1)=fvh0
  endif
  endif ! Refinement
end subroutine scanspacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bisectspacing(vs1fac,vs2fac)
! Given vs1 and vs2 spacing factors that are below and above the 
! stability threshold, bisect to refine the threshold estimate.
  vs3fac=(vs1fac+vs2fac)/2.
  vs=vss*vs3fac
  call scanvh(index)   ! Sets the vh0
  dFdx0=(delFdxvh(index-1)*abs(Forcevh(index)) &
       +delFdxvh(index)*abs(Forcevh(index-1)))/ &
       (abs(Forcevh(index-1))+abs(Forcevh(index)))
  if(index.ne.0.and.dFdx0.le.0.)then
     vs2fac=vs3fac
     dvh=(vhmax-vhmin)/1.
     vhmin=vh0-dvh/2.
     vhmax=vh0+dvh/2.
     vh=vh0 
     call finddelphi  ! For the print out:
     write(*,'(i4,7f8.4,f9.5)')index,vha(index-1),vha(index),&
          Forcevh(index-1)/phimax**2,Forcevh(index)/phimax**2, &
             dFdx0/phimax**2,vh0,vs(1),delphi/phimax**2
  else
     vs1fac=vs3fac
  endif
end subroutine bisectspacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initvofv
! Initialize velocity range using current vh as maximal.
  vmax=4.4+max(abs(maxval(vs)-vh),abs(minval(vs)-vh))
  do i=1,nofv
     vofv(i)=-vmax+2.*vmax*(i-1.)/(nofv-1.)
  enddo
end subroutine initvofv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findtrappedf   
! Find the trapped distrib on prescribed phi side. And plot it. 2 plots.
  external denionhill
  isigma=nint(sign(1.,delphi))   ! Side of hill with highest phi_infty.
  coshlen=4.              ! Set the shape adjustment parameters.
  toplen=-10.             ! Negligible elongation when large and negative.
  psi=phimax-isigma*delphi/2.
  um=0.
  call BGKint(nphi,psi,um,xmax,coshlen,toplen,denion,phiofi(0),usofphi &
  ,xofphi(0),edenofphi(0),edenuntrap,edentrap,etilden,fofphi,u0ofphi,denionhill)
!  write(*,*)'Return from BGKint',u0ofphi(0),u0ofphi(nphi)
!  write(*,'(i4,2f8.4)')(i,xofphi(i),etilden(i),i=0,nphi)
  call multiframe(3,1,0)
  call autoplot(xofphi(0),phiofi(0),nphi+1)
  call axlabels('','!Af!@-!Af!d;!d!@')
  call autoplot(xofphi(0),edenofphi(0),nphi+1)
  call axlabels('','n!de!d')
  call autoplot(xofphi(0),etilden,nphi+1)
  call axlabels('x','          !pn!q!o~!o!q!de!d')
  call pltend
  call multiframe(0,0,0)
  call autoplot(u0ofphi*sqrt(2.),fofphi/sqrt(2.),nphi+1)
  call axlabels('v!de!d  [/(T!de!d/m!de!d)!u1/2!u]','f!det!d(v!de!d)')
  call pltend
end subroutine findtrappedf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine otherxofphi  ! Solve Poisson's equation on the other side.
! and plot the comparison.
  isigma=-nint(sign(1.,delphi))
  phistep=(phimax+isigma*delphi/2.)/nphi      ! The spacing of the first side.
  nextra=max(2,min(nint(abs(delphi)/phistep),nphi))! The number of extra values.
  phistep=abs(delphi)/nextra
!  write(*,*)'isigma,nextra,delphi',isigma,nextra,delphi
  do i=-nextra,-1
     phiofi(i)=i*phistep
     edenofphi(i)=exp((phiofi(i)+abs(delphi)/2.)/Te) ! Maxwellian e-density
  enddo  ! phiofi is now the full potential relative to first side's phiinfty.
  do i=-nextra,nphi
     phiofi(i)=phiofi(i)+abs(delphi) ! Refer to this side for denionhill.
  enddo
  ! Integrate from phi=phimax:phiinfty, i= nphi:-nextra, x= 0:isigma*infty
  rhoprior=denionhill(phiofi(nphi))-edenofphi(nphi)
  ! Initialize nonzero Vminus at phi=phimax to correct first step.
  Vminus(nphi)=(phiofi(nphi)-phiofi(nphi-1))*rhoprior/9.
  x2ofphi(nphi)=0.
  ninfty=-nextra
  do i=nphi-1,-nextra,-1
     rhothis=denionhill(phiofi(i))-edenofphi(i)     ! Charge density
     Vminus(i)=Vminus(i+1)-(phiofi(i)-phiofi(i+1))*(rhoprior+rhothis)/2.
     if(Vminus(i).le.0)then
        write(*,*)'Reached zero Vminus',i,Vminus(i)
!        x2ofphi(i)=isigma*xmax
        x2ofphi(i)=xmax
        ninfty=i
        exit
     endif
!     x2ofphi(i)=x2ofphi(i+1)-isigma*(phiofi(i)-phiofi(i+1))* &
     x2ofphi(i)=x2ofphi(i+1)-(phiofi(i)-phiofi(i+1))* &
          (sqrt(1./(2.*Vminus(i)))+sqrt(1./(2.*Vminus(i+1))))/2.
     rhoprior=rhothis
     if(lcd.and.i.le.2)write(*,'(i4,4f10.6)')i,phiofi(i)-abs(delphi/2)&
          &,rhothis,Vminus(i),x2ofphi(i)
  enddo
  if(lcd)write(*,*)'  i     phi      rho       Vminus       x'
! For illustrative reasons make the penultimate point an exponential decay
! This is a fudge to cover up the slight inaccuracy is the V integration.

  x2ofphi(ninfty+1)=x2ofphi(ninfty+2)+1.*alog(2.)
  x2ofphi(ninfty)=xmax
  phiofi=phiofi-abs(delphi/2.)  ! Refer phi to zero.
  
!  write(*,*)'Final phiinfs',delphi/2.,phiofi(0),phiofi(ninfty)
  call autoplot(x2ofphi(ninfty),phiofi(ninfty),nphi-ninfty+1)
  call axis2
  call axlabels('|x|','!Af!@(x)')
  call legendline(.6,.9,0,' -!As!@!dm!d')
  call dashset(1)
  call color(1)
  call legendline(.6,.8,0,'  !As!@!dm!d')
  call polyline(xofphi(0),phiofi(0),nphi+1)
  call dashset(0)
  call pltend
end subroutine otherxofphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotfvden
call multiframe(2,4,0)
call dcharsize(.025,.025)
call pltinit(vha(1),vha(nvh),0.,0.5)
call ticlabtog
call axis
call ticlabtog
call axis2
call polyline(vha,fofvh,nvh)
call legendline(.05,.9,258,'!Bf!di!A;!@!d(!Bv!@!d!A;!@!d)')
call legendline(.85,.08,258,'!Bv!@!d!A;!@!d')
call legendline(.05,.8,258,'0.4')
if(dFdx0.lt.0)call polymark(vh0,fvh0,1,10)

call setframe(2)
call minmax(denofx,npts,dmin,dmax)
call fitrange(dmin,dmax,3,ipow,fac10,delta,first,xlast)
call pltinit(x(1),x(npts),first,xlast)
call ticlabtog
call axis
call ticlabtog
call axis2
do i=1,2
   yl=first+i*delta
   call fwrite(yl,iwidth,3,string)
   call jdrwstr(wx2nx(x(npts)),wy2ny(yl),string(1:iwidth),-1.2)
enddo
call polyline(x,denofx,npts)
call jdrwstr(wx2nx(x(npts)),wy2ny(.1*first+.9*xlast),'!Bn!di!@!d(!Bx!@)!@',-1.2)
call legendline(.9,.08,258,'!Bx!@')
call multiframe(0,0,0)
end subroutine plotfvden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotdenofx
! Plot phi,ni,V of x.
  call multiframe(3,1,0)
  call autoplot(x,phiofx,npts)
  call legendline(.9,.9,258,'(a)')
  call axis2
  call polyline([x(1),x(npts)],[0.,0.],2)
  call fwrite(vh,iwidth,2,string)
  call legendline(.1,1.1,258,'v!dh!d='//string(1:iwidth))
  call axlabels('','!Af!@(x)')
  call autoplot(x,denofx,npts)
  call legendline(.9,.9,258,'(b)')
  call axis2
  call polyline([x(1),x(npts)],[denave,denave],2)
  call axlabels('','n!di!d(x)')
  if(ldenion)then    ! Testing of denionhill.
     ixf=int(npts*.35)
     isigma=-1
     den1=denionhill(phiofx(ixf)-isigma*delphi/2.)
     isigma=1
     den2=denionhill(phiofx(npts-ixf)-isigma*delphi/2.)
     call polymark([x(ixf),x(npts-ixf)],[den1,den2],2,2)
  endif
  call winset(.true.)
  call color(4)
!  call polyline(x,denave*exp(phiofx/Te),npts)
  call polymark([x(1),x(npts)], &
       denave*[exp(phiofx(1)/Te),exp(phiofx(npts)/Te)],2,1)
  call color(15)
  call autoplot(x,deninteg,npts)
  call legendline(.9,.9,258,'(c)')
  call axis2
  call axlabels('x','!AJ!@n!di!dd!Af!@')
  call polyline([x(1),x(npts)],[0.,0.],2)
  call pltend
  call multiframe(0,0,0)
end subroutine plotdenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotexplain
! Explanatory plot with reflected/passing labels etc.
  call multiframe(4,1,1)
  call minmax(phiofx,npts,pmin,pmax)

  call pltinit(-xmax,xmax,pmin*2.,pmax*1.4)
  call legendline(.95,.1,258,'(a)')
  call axis
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(1)
  call axlabels('','Ion Energy');call axptset(0.,0.);call ticrev()
  call polyline([-.91,-.91]*xmax,[pmin,pmax],2)
  call polyline([-.93,-.89]*xmax,[pmax,pmax],2)
  call jdrwstr(wx2nx(-.9*xmax),wy2ny(pmin)+.05,'Reflected!A_!@',1.05)
  call polyline([.9,.9]*xmax,[phiofx(npts),pmax],2)
  call polyline([.92,.88]*xmax,[pmax,pmax],2)
  call jdrwstr(wx2nx(.9*xmax),wy2ny(phiofx(npts))+.03,'!A^!@Reflected',-1.05)
  call polyline([0.,0.],[pmax,pmax*1.4],2)
  call jdrwstr(wx2nx(.02*xmax),wy2ny(pmax)+.01,'Passing!A^_!@',1.)
  call color(15)
  call polyline(x,phiofx,npts)
  call polyline([x(1),x(npts)],[0.,0.],2)
  call axlabels('','!Af!@(x)')
  
  call pltinit(-xmax,xmax,-pmax*1.1,pmax)
  call legendline(.95,.1,258,'(b)')
  call axis
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(4)
  call axlabels('','Electron Energy');call axptset(0.,0.);call ticrev()
  call polyline([-.91,-.91]*xmax,[-pmin,pmax],2)
  call jdrwstr(wx2nx(-.9*xmax),wy2ny(-pmin)+.02,'Passing!A^_!@',1.05)
  call polyline([0.05,0.05],[-pmax,-phiofx(npts)],2)
  call jdrwstr(wx2nx(0.),wy2ny(-pmax+.04),'Trapped!A^!@ ',-1.0)
  call drcstr(' !A_!@')
  call polyline([-.15,.25],[-phiofx(npts),-phiofx(npts)],2)
  call polyline([.9,.9]*xmax,[-phiofx(npts),-pmin],2)
  call jdrwstr(wx2nx(.9*xmax),wy2ny(-phiofx(npts))+.015,'!A^!@Reflected',-1.05)
  call polyline([.92,.88]*xmax,[-pmin,-pmin],2)
  call color(15)
  call polyline(x,-phiofx,npts)
  call polyline([x(1),x(npts)],[0.,0.],2)
  call axlabels('','-!Af!@(x)')

  call minmax(denofx,npts,dmin,dmax)
  call pltinit(-xmax,xmax,dmin-(dmax-dmin)*.07,dmax+(dmax-dmin)*.07)
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

  call minmax(deninteg,npts,dmin,dmax)
  call pltinit(-xmax,xmax,dmin,dmax+(dmax-dmin)*.07)
  call axis
  call polyline(x,deninteg,npts)
  call legendline(.95,.1,258,'(d)')    
  call axptset(1.,0.);call ticrev();call altyaxis(1.,1.0);call color(1)
  call axptset(0.,0.);call ticrev()
  call axlabels('x','!AJ!@n!di!dd!Af!@')
  call polyline([x(1),x(npts)],[0.,0.],2)
  call pltend
  call multiframe(0,0,0)
end subroutine plotexplain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine plotfvhill
! Plot f(v) for sigma=+-, and P(v)
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
end subroutine plotfvhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forcedivplot
! Plot the different contributions to F from fofvh Forcevh as fn of vha.
! And if there's an equilibrium mark it. 
  call multiframe(2,1,3)
  call dcharsize(.02,.02)
  call autoplot(vha,fofvh,nvh)
  call axlabels('','f!di!A;!@!d(v!d!A;!@!d)')
  call legendline(.6,-.1,258,'v!d!A;!@!d')
  call axis2
  call legendline(.9,.9,258,'(a)')
  if(dFdx0.lt.0)call polymark(vh0,fvh0,1,10)
   ! Electron force
  call autoplot(vha,denave*Te*(exp(dphivh/(2.*Te))-exp(-dphivh/(2.*Te))),nvh)
  call axlabels('','Forces(v!dh!d)')
  call legendline(.6,-.1,258,'v!dh!d')
  call axis2
  call legendline(.05,.8,0,'F!de!d')
  call legendline(.9,.9,258,'(b)')
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
  call legendline(.9,.9,258,'(a)')
  call axlabels('','f!d!A;!@!d(v!d!A;!@!d)')
  do i=1,ns
     call color(i)
     call polyline(vha,fvvhns(1,i),nvh)
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),fv0ns(i),1,10)
  enddo
  call color(15)
  if(lrefinethreshold)then
     call polymark(vha,fvvhns(1:nvh,ns+1),nvh,4)
     dsy=(symax-symin)/10.
     call polyline([vh0ns(ns+1),vh0ns(ns+1)],&
          [fv0ns(ns+1)+dsy,fv0ns(ns+1)-dsy],2)
  endif
  
  call minmax(Fcvhns(1:nvh,1)/phimax**2,nvh,symin,symax)
  call pltinit(vhmin,vhmax,symin*1.1,symax*1.1)
  call legendline(.9,.9,258,'(b)')
  call axis
  call axis2
  call axlabels('','F(v!dh!d)/!Ay!@!u2!u')  
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
!  call legendline(.05,.7,0,'F!di!d+F!de!d')
  do i=1,ns
     call color(i)
     call polyline(vha,Fcvhns(1:nvh,i)/phimax**2,nvh)      ! Total force
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),0,1,10)
  enddo
  call color(15)
  if(lrefinethreshold)then
     call polymark(vha,Fcvhns(1:nvh,ns+1)/phimax**2,nvh,4)
     dsy=(symax-symin)/10.
     call polyline([vh0ns(ns+1),vh0ns(ns+1)],&
          [dsy,-dsy],2)
  endif
  call minmax(dFvhns(1:nvh,1)/phimax**2,nvh,symin,symax)
  call pltinit(vhmin,vhmax,symin*1.1,symax*1.1)
  call legendline(.9,.9,258,'(c)')
  call axis
  call axis2
  call axlabels('v!dh!d  or  v!d!A;!@!d','!Ad!@F/!Ad!@x(v!dh!d)/!Ay!@!u2!u')  
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
  do i=1,ns
     call color(i)
     call polyline(vha,dFvhns(1:nvh,i)/phimax**2,nvh)      ! Total force
     if(dF0ns(i).lt.0)call polymark(vh0ns(i),dF0ns(i)/phimax**2,1,10)
  enddo
  call color(15)
  if(lrefinethreshold)then
  call polymark(vha,dFvhns(1:nvh,ns+1)/phimax**2,nvh,4)
  dsy=(symax-symin)/10.
  call polyline([vh0ns(ns+1),vh0ns(ns+1)],&
       [dF0ns(ns+1)/phimax**2+dsy,dF0ns(ns+1)/phimax**2-dsy],2)
  endif
  call pltend
  call multiframe(0,0,0)  
end subroutine forceplotns
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
!     vthresh=isigma*sqrt(max(0.,phimx2-phix2)) ! prevent rounding NANs.
  if(phimx2.lt.phix2)then    ! Move the threshold far away.
     vthresh=40.*(vofv(nofv)-vofv(1))
  else
     vthresh=isigma*sqrt(phimx2-phix2)
  endif
  Pfofv=0.
  fofv=0.
  vdiffm=vofv(1)-vthresh
  fim=0.
  pam=0.
  pa=0.
  sqtwopi=sqrt(2.*3.1415926)
  do i=1,nofv
     vsign=sign(1.,vofv(i))
     vdiff=vofv(i)-vthresh
     enx2=phix2+vofv(i)**2
     if(int(vsign).ne.isigma)then              ! Moving inward
        vinf=vsign*sqrt(max(0.,enx2))
     elseif((phix2+vofv(i)**2).gt.phimx2)then  ! Passing
        vinf=vsign*sqrt(max(0.,delphix2+enx2))
     else                                      ! Reflected
        vinf=-vsign*sqrt(max(0.,enx2))
     endif
     fi=0.
     do j=1,nc   ! f(v)=finf(vinf)
        fi=fi+dc(j)*exp(-0.5*(vinf-vs(j))**2/vt(j)**2)/(vt(j)*sqtwopi)
     enddo
!     if(enx2.lt.0)fi=0.    ! Zero f at negative energy??
     fofv(i)=fi
     if(i.gt.1)pa=pam+ (abs(vdiffm)*fim+abs(vdiff)*fi) &
             /(abs(vdiffm)+abs(vdiff))*(vofv(i)-vofv(i-1))
! Actually using this interpolation always smooths ftrapped glitches.
! But I am not 100% certain it is always justified.
     Pfofv(i)=Pfofv(i)+pa
     pam=pa
     vdiffm=vdiff
     fim=fi
  enddo
end subroutine fvhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate density for shifted Maxwellian
! distant distributions on a single-humped repelling potential hill of
! maximum potential phimax, at potential phi, on side isigma of the hill.
! This version uses a velocity mesh for integration that always has a 
! point at the discontinuity of fv, and does not return fv, only density.
! On entry:
! nc is the number of Maxwellian components 
! dc(nc), vs(nc), vt(nc) is each component's density, vshift, sqrt(T/M)
! phimax is hill's peak potential relative to mean distant potential
! phi is the potential at which to find the distribution.
! isigma is the side of the hill on which this potential lies.
! All velocities are relative to the hill's rest frame.
! On exit:
! theden is the integral of fv dv. I.e. the density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine denhill(nc,dc,vs,vt,vh,phimax,phi,isigma,nofv,theden)
  real :: dc(nc),vs(nc),vt(nc),phimax,phi,theden
  integer :: isigma,nofv

  
  if(isigma.ne.1.and.isigma.ne.-1)stop 'isigma must be +-1'
  phix2=2.*phi
  phimx2=2.*phimax
  vthresh=isigma*sqrt(max(0.,phimx2-phix2)) ! prevent rounding NANs.
  vinf=vthresh ! Silence warning
  vmax=4.4+max(abs(maxval(vs)-vh),abs(minval(vs)-vh))
  sqtwopi=sqrt(2.*3.1415926)
  nvhere=nofv+1-mod(nofv,2) ! Odd number of points so thresh is at (nvhere+1)/2
  ithresh=(nvhere+1)/2
  theden=0.
  fim=0.
  vm=(ithresh-1)*(-vmax-vthresh)/(ithresh-1.)+vthresh
  do i=1,nvhere
     if(i.le.ithresh)then  ! Lower and upper velocity ranges
        vi=(ithresh-i)*(-vmax-vthresh)/(ithresh-1.)+vthresh
     else                  ! vthresh to +-vmax
        vi=(i-ithresh)*(vmax-vthresh)/(nvhere-ithresh)+vthresh
     endif
     vsign=sign(1.,vi)
     enx2=phix2+vi**2
     if(i.eq.ithresh)then  ! Crossing threshold. Use prior sign this step
        vinf=sign(sqrt(max(0.,enx2)),vinf)
     elseif(int(vsign).ne.isigma)then        ! Moving inward
        vinf=vsign*sqrt(max(0.,enx2))        
     elseif(enx2.gt.phimx2)then              ! Passing
        vinf=vsign*sqrt(max(0.,enx2))
     else                                     ! Reflected
        vinf=-vsign*sqrt(max(0.,enx2))
     endif
     fi=0.
     do j=1,nc   ! Sum over maxwellian components.
        fi=fi+dc(j)*exp(-0.5*(vinf-vs(j))**2/vt(j)**2)/(vt(j)*sqtwopi)
     enddo
     theden=theden+ (fim+fi)*0.5*(vi-vm)
     vinf=-vinf
     if(i.eq.ithresh)then   ! Set the new prior fi for next step
        fi=0.
        do j=1,nc   ! Sum over maxwellian components.
           fi=fi+dc(j)*exp(-0.5*(vinf-vs(j))**2/vt(j)**2)/(vt(j)*sqtwopi)
        enddo
     endif
!     write(*,'(i5,5f10.6)')i,vthresh,vi,fi,theden
     fim=fi
     vm=vi
  enddo
end subroutine denhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the ion density at a potential phiminf relative to phiinf on
! the side isigma. To be called by BGKint.
real function denionhill(phiminf)
  use AsymHill
  phimzero=phiminf+isigma*delphi/2.
  if(ldenalt)then   ! This density must be consistent with delphi.
     call denhill(nc,dc,vs-vh,vt,vh,phimax,phimzero,isigma,nofv,denionhill)
  else
     call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phimzero,isigma,local, &
          nofv,vofv,fofv,Pfofv)
     denionhill=Pfofv(nofv)
  endif
end function denionhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Asym
use AsymHill

call parseAsymargs
call initvofv

if(.false.)then
! Testing force calculation
!vh=5.e-6
call finddenofx
do i=1,npts/2
write(*,*)i,denofx(i)-denofx(npts+1-i),deninteg(i),deninteg(npts+1-i)
enddo
call finddelphi
! This shows that we are in rounding noise at moderately low psi.
endif

if(ltestnofx)then
   call finddelphi
   call plotfvhill
!   write(*,'(a,9f6.3)')'plotfvhill: vshift,vthm,dens',(vs(i),vt(i),dc(i),i=1,nc)
   call finddenofx
   if(lexplain)call plotexplain
   call plotdenofx
endif

if(ns.gt.0)then
call scanspacing
call forceplotns
endif


call scanvh(index)
call plotfvden
call pltend

end program Asym
