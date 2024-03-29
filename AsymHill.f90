module AsymHill
  integer, parameter :: npts=200,nofv=200,nvh=100,ncmax=4,nsmax=20,npsi=10
  real, dimension(npts) :: x,phiofx,denofx,deninteg,dFdelx
  real, dimension(nofv) :: vofv,fofv,Pfofv,Forceofv
  real, dimension(nvh) :: vha,Forcevh,delFdxvh,dphivh,fofvh
  real, dimension(nvh) :: analFrci,analdendiff,analdendift
  real, dimension(nvh) ::  directddiff,directddift
  real, dimension(ncmax) :: vs=[1.6,-1.6,0.,0.],vt=1.,dc=[.7,.3,0.,0.],vss
  real, dimension(nvh,nsmax) :: Fcvhns,dFvhns,fvvhns,dpvhns
  real, dimension(nsmax) :: fv0ns,dF0ns,vh0ns,delphins
  real, dimension(npsi) :: psia,delpha,analdp
  real :: delphi,phimax=.2,phi=.1,xmax=12.,delphi0,Forcediff,dendiff,diffint
  integer :: jit
  character*30 string,argument
  real :: vh=0.,Te=1.,denave,vh0,fvh0,dFdx0,vhmin=-4.,vhmax=4.,vhn,vhx,vh0r
  integer :: index,nc=2,ns=1,pfint=0,isigma
  logical :: local=.false.,lcd=.true.,ltestnofx=.false.,ldenion=.false.
  logical :: lrefinethreshold=.true.,lfdconv,lexplain=.false.
  logical :: ldenalt=.true.,lfvden=.false.,lfover=.false.,lfindforce=.false.
! lfindforce breaks some things at the moment.        
  logical :: ldscale=.false.,ltesting=.false.
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
       if(argument(1:2).eq.'-m')lfvden=.true.
       if(argument(1:2).eq.'-x')lexplain=.not.lexplain
       if(argument(1:2).eq.'-d')ldenalt=.not.ldenalt
       if(argument(1:2).eq.'-o')lfover=.not.lfover
       if(argument(1:2).eq.'-a')ldscale=.not.ldscale
       if(argument(1:2).eq.'-b')lfindforce=.not.lfindforce
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
    write(*,'(a,l4)'  )'  -m...  Multi-T,n scan                 [',lfvden
    write(*,'(a,l4)'  )'  -x...  Plot f(v) explanation          [',lexplain
    write(*,'(a,l4)'  )'  -d...  Other density algorithm toggle [',ldenalt
    write(*,'(a,l4)'  )'  -o...  f(v) overlay plot toggle       [',lfover
    write(*,'(a,l4)'  )'  -a...  delphi scale analysis          [',ldscale
    write(*,'(a,l4)'  )'  -b...  Use findforce                  [',lfindforce
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
  dpdr=0.
  dx=2.*xmax/(npts-1.)
  ximid=npts/2+.5*mod(npts+1,2)
  do i=1,npts
!     x(i)=-xmax+(i-1.)*dx
     x(i)=(i-ximid)*dx     ! Symmetrized x version
     isigma=int(sign(1.,x(i)))
     dph=isigma*delphi/2.
     xcor=abs(x(i)/4.)
!     xcor=x(i)/4./sqrt(1.-dph/phimax)  !Symmetrizes curvature at origin.
     phiofx(i)=(phimax-dph)/cosh(xcor)**4+dph
     if(ldenalt)then   ! This density must be consistent with delphi.
        call denhill(nc,dc,vs-vh,vt,vh,phimax,phiofx(i),isigma,nofv&
             &,denofx(i),dendiff)
     else
        call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiofx(i),isigma,local, &
             nofv,vofv,fofv,Pfofv)
        denofx(i)=Pfofv(nofv)
     endif
     if(i.eq.1)then
        deninteg(i)=denofx(i)*(phiofx(i)-dph)
     else
        deninteg(i)=deninteg(i-1)+&      !\int n dphi
             (denofx(i)+denofx(i-1))*(phiofx(i)-phiofx(i-1))*.5
     endif
     if(i.eq.1)then
        dFdelx(i)=-denofx(i)*(phiofx(i)-phiofx(i-1))/(x(i)-x(i-1))
     else
        dFdelx(i)=dFdelx(i-1)-&          !-\int (dn/dx) dphi
          (denofx(i)-denofx(i-1))*(phiofx(i)-phiofx(i-1)) &
          /(x(i)-x(i-1))
     endif
!     if(x(i).gt.0)write(*,*)i,x(i),x(1+npts-i)   ! Show symmetry.
!     if(x(i).gt.0)write(*,*)i,x(i),x(1+npts-i),phiofx(i),phiofx(1+npts-i)
  enddo
  deninteg(npts)=deninteg(npts)-denofx(npts)*(phiofx(npts)-dph)
  if(lcd)write(*,'(a,f8.5,a,f9.5,a,f9.5,a,f9.5)')'vh=',vh,' delphi=',delphi,&
       ' Ion force=',-deninteg(npts),' delF/delx=',dFdelx(npts)
end subroutine finddenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddelphi
! Iterate to find delphi consistent with specified ion distributions.
  real, dimension(2) :: deninf
  integer, parameter :: niter=15
  real, dimension(niter) :: residit,delphiit
  if(lcd)write(*,*)'j   dpp       delphi      delphi-dpp    n(-)        n(+)',&
       '   resid'
  lfdconv=.false.
  dpp=delphi
  ddp=0.
  resid=0.
  do jit=1,niter
     residp=resid
     do i=1,2
        isigma=-3+2*i
        phiinf=isigma*delphi/2.
        if(ldenalt)then ! Alternative ni calculation.
           call denhill(nc,dc,vs-vh,vt,vh,phimax,phiinf,isigma,nofv&
                &,deninf(i),dendiff)
        else
           call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiinf,isigma,local,&
                nofv,vofv,fofv,Pfofv)
!           write(*,'(a,i4,5f10.6)')'Findelphi',jit,phiinf,deninf(i),Pfofv(nofv)
           deninf(i)=Pfofv(nofv)
           if(i.eq.2)dendiff=deninf(2)-deninf(1)
        endif
     enddo
     resid=alog(1.+dendiff/deninf(1))-delphi/Te
!        resid=(alog(deninf(2)/deninf(1))-delphi/Te)
     residit(jit)=resid
     if(jit.eq.1)then
        delphi=delphi+.5*resid*min(1.,Te)
     else
        if((resid-residp).eq.0)then
           write(*,'(a,i3,5e12.4)')'resid unchanged',jit,resid,delphi,dpp,ddp
           dpp=delphi
           delphi=delphi+.5*resid*min(1.,Te)
        else
           dpdr=(delphi-dpp)/(resid-residp)
           ddp=resid*dpdr
           dpp=delphi
           delphi=delphi-sign(min(abs(ddp),phimax/2.),ddp)
        endif
     endif
     delphiit(jit)=delphi
     if(lcd)write(*,'(i2,6f11.7)')jit,dpp,delphi,delphi-dpp,deninf,resid
     if(abs(resid).lt.0.1e-6)then
        lfdconv=.true.
        goto 2
     endif
  enddo
  if(.not.abs((delphi-dpp)/delphi).lt.1.e-10)then
     write(*,*)'finddelphi unconverged. residit,delphiit:'
     write(*,'(5f11.7)')(residit(j),j=1,niter)
     write(*,'(5f11.7)')(delphiit(j),j=1,niter)
     write(*,'(a,7f11.7)')' delphi,dpp,ddp',delphi,dpp,ddp
     call automark(residit(1),delphiit(1),niter-1,1)
     call pltend
  endif
2 continue
  denave=deninf(1)*exp(0.5*delphi/Te)
!  denave=(deninf(1)+deninf(2))/2.   ! Makes very little difference.
end subroutine finddelphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findforce
! Integrate the intrinsic ion force symmetrically over phi>abs(delphi/2)
! giving \int dendiff d\phi, from phimax to abs(delphi/2.)
  npt2=npts/2
  diffint=0.
  dendiff=0.
  phip=0
! Instrinsic ion force
  do i=1,npt2
     dendiffp=dendiff
     phi=(phimax-abs(delphi/2.))/(npt2-1.)*(npt2-i)+abs(delphi/2.)
     do j=1,2         ! dendiff between positive and negative isigma
        isigma=-3+2*j
        call denhill(nc,dc,vs-vh,vt,vh,phimax,phi,isigma,nofv,denhere,dendiff)
     enddo
     diffint=diffint-(dendiff+dendiffp)*0.5*(phi-phip)  !Ion force integral.
     phip=phi
  enddo
! Extrinsic ion force. Maybe?
if(.false.)then
  do i=1,npt2
     phi=abs(delphi/2)*(1.-2.*(i-1.)/(npt2-1))
     isigma=-nint(sign(1.,delphi))
     call denhill(nc,dc,vs-vh,vt,vh,phimax,phi,isigma,nofv,denhere,dendiff)
     diffint=-(dendiff+dendiffp)*0.5*(phi-phip)  !Ion force integral.
     phip=phi
     dendiffp=dendiff
  enddo
endif   
! Extrinsic i-e force difference. 
! Naively the ion density times \delphi using neutrality
! balanced against the electron pressure difference. 
  Fext=denave*(-delphi+Te*2.*sinh(delphi/(2.*Te)))
!  write(*,*)'vh,denave,Fext',vh,denave,Fext
  Forcediff=diffint + Fext
!  write(*,*)'findforce',Forcediff,delphi,diffint,Fext
end subroutine findforce
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
  nref=1
  do j=1,nref   ! Maybe iterate for additional refinement of vh0. Not usually
     Force=0.
     Fmax=0.
     vh=vhn
     indexref=0
     do i=1,nvh          ! Scan over the current vh range
        vhprior=vh
        fvh0prior=fvh0
        delphiprior=delphi
        Fprior=Force
        vh=vhn+(i-1.)*(vhx-vhn)/(nvh-1.)
        lcd=.false.
        call finddelphi
        call analcomp(i)
        if(lfindforce)then
           call findforce
           Force=Forcediff
        else
           call finddenofx
           Force=-deninteg(npts) &
                +denave*Te*2.*sinh(delphi/(2.*Te))
!             +denave*Te*(exp(delphi/(2.*Te))-exp(-delphi/(2.*Te)))
! Fill in diffint for intrinsic force diff. For testing
           Fext=denave*(-delphi+Te*2.*sinh(delphi/(2.*Te)))
           diffint=Force - Fext
        endif
        if(vhd-abs(vh-vhm).gt.vht  &    ! Not too near the ends
             .and.Fprior.gt.0.and.Force.lt.0..and.index.eq.0)then
           indexref=i
           vh0=(vhprior*abs(Force)+vh*abs(Fprior))/ &
                (abs(Fprior)+abs(Force))
           delphi0=(delphiprior*abs(Force)+delphi*abs(Fprior))/ &
                (abs(Fprior)+abs(Force))
           vh1=vhprior
           vh2=vh
        endif
        if(j.eq.1)then    ! Store coarse scan.
           if(lfindforce)call finddenofx ! Needed for dFdelx at the moment.
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
!     write(*,*)'vh0, F, dp',vh0,Force,delphi0
     if(indexref.gt.1)then  ! Refine the equilibrium vh0.
        vhn=vh1
        vhx=vh2
!        write(*,*)'vhn,vhx',vhn,vhx
     else
        exit ! from refinement j-loop.
     endif
  enddo
  if(indexref.gt.1.and.nref.eq.1)then  ! Refine the equilibrium vh0.
     vhn=vh1
     vhx=vh2
     do i=1,20
!        write(*,*)'vhn,vhx,Frc',vhn,vhx,Forcediff
        call bisectvh(vhn,vhx)
     enddo
     delphi0=delphi
  endif
end subroutine scanvh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bisectvh(vhn,vhx)
  ! Given vhn and vhx at which Force is +ve and -ve, bisect to find
  ! a better pair.
  vh=(vhn+vhx)/2.
  vh0=vh
  call finddelphi
  if(lfindforce)then
     call findforce
     Force=Forcediff
  else
     call finddenofx
     Force=-deninteg(npts) &
          +denave*Te*2.*sinh(delphi/(2.*Te))
!             +denave*Te*(exp(delphi/(2.*Te))-exp(-delphi/(2.*Te)))
  endif
  if(Force.lt.0.)then
     vhx=vh
  else
     vhn=vh
  endif
end subroutine bisectvh
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
!     vs1fac=0.  ! This seems strange. Why do it?
     dvh=(vhmax-vhmin)/5.
     vhmin=vh0-dvh
     vhmax=vh0+dvh
     do i=1,ns*2
!        write(*,*)vs1fac,vs2fac
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
     edenofphi(i)=denave*exp((phiofi(i)+abs(delphi)/2.)/Te) ! Maxwellian ne.
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
real function f1prime(v)
  f1prime=0.
  do i=1,nc
     vdiff=(v-vs(i))
     f1prime=f1prime-dc(i)*vdiff/vt(i)**3*exp(-(vdiff/(vt(i)))**2/2.)
  enddo
  f1prime=f1prime/sqrt(2.*3.1415926)
end function f1prime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function f3prime(v)
  f3prime=0.
  do i=1,nc
     vdiff=(v-vs(i))
     f3prime=f3prime+dc(i)*(3.-(vdiff/vt(i))**2)*vdiff/vt(i)**5&
          *exp(-(vdiff/(vt(i)))**2/2.)
  enddo
  f3prime=f3prime/sqrt(2.*3.1415926)
end function f3prime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotfvden
! For 2 components with symmetric velocities +-(vsin(1)+vt(2))
! Scan the second component's density n2 and temperature T2 (=vt(2)^2),
! plotting in a 2D multiframe arrangement the f(v) and ni(x)
! for psi values psimax*10^{-(i-1)} for i=1,npscan
nn2=3
nT2=3
npscan=3
den2max=0.5
T2max=1.
phimxin=phimax
vs1in=vs(1)
vs2in=vs(2)
vt1in=vt(1)
ncin=nc
nc=2

call multiframe(nn2,2*nT2,0)
call dcharsize(.025,.025)
do k=1,nT2
   T2=T2max*(0.2+0.8/max(nT2-1,1)*(k-1))
   vt(1)=1.
   vt(2)=sqrt(T2)
   vs(1)=vs1in+vt(2)
   vs(2)=-vs(1)
   if(k.gt.1)write(*,*)
   write(*,'(a,2f10.4,a,$)')'vs(1),vt(2)',vs(1),vt(2),',   dc(2)'
   do j=1,nn2
      den2=den2max*(0.2+0.8/max(nn2-1,1)*(nn2-j))
      dc(2)=den2
      dc(1)=1.-den2
      write(*,'(f10.4,$)')dc(2)
      do i=1,npscan
         phimax=phimxin/10.**(i-1)
         call fwrite(T2,iwidth,1,argument)
         string='T!d2!d='
         string(lentrim(string)+1:)=argument(1:iwidth)
         call scanvh(index)
! Left subframe
         iframe=(mod(j-1,nn2)+2*nn2*(k-1))
         call setframe(iframe)
         if(i.eq.1)then
            call pltinit(vha(1),vha(nvh),0.,0.5)
         else
            call scalewn(vha(1),vha(nvh),0.,0.5,.false.,.false.)
         endif
         call ticlabtog
         call axis
         call ticlabtog
         call legendline(.15,.9,258,'!Bf!di!A;!@!d(!Bv!@!d!A;!@!d)')
         call legendline(.6,.08,258,'!Bv!@!d!A;!@!d')
         call legendline(.05,.8,258,'0.4')
         if(j.eq.nn2)call legendline(.8,-.1,258,string)
         if(i.eq.1)call polyline(vha,fofvh,nvh)
         call color(max(1,2*(i-1))) !call color(2*i-1)
         if(dFdx0.lt.0)call polymark(vh0,fvh0,1,10)
         call color(15)
         call fwrite(den2,iwidth,1,argument)
         string='n!d2!d='
         string(lentrim(string)+1:)=argument(1:iwidth)
         if(k.eq.1)call legendline(-.5,0.5,258,string)

! Right subframe
         call setframe(iframe+nn2)
         if(i.eq.1)then
            denmax=0.6
            denmin=-0.2
            dendelta=-2*denmin
            call fitrange(denmin,denmax,3,ipow,fac10,delta,first,xlast)
         endif
!         call pltinit(x(1),x(npts),first,xlast)
         call pltinit(x(1),x(npts),denmin-.05,denmax)
         call ticlabtog
         call xaxis(0.,0.)
         call ticlabtog
         call axptset(1.,0.)
         call ticrev
         call yaxis(first,delta)
         call ticrev
         call axptset(0.,0.)
         call axbox
!         call legendline(-.2,.65,258,string)
         do ii=1,2   ! y-axis scale
            yl=denmin+(i-1)*dendelta
            if(yl.lt..8*denmax)then
               call fwrite(yl,iwidth,1,string)
               call jdrwstr(wx2nx(x(npts)),wy2ny(yl),string(1:iwidth),-1.4)
            endif
         enddo
         call color(max(1,2*(i-1)));call dashset(2*(i-1))
         if(dFdx0.lt.0)call polyline(x,(denofx-1)/phimax,npts)
         call fwrite(phimax,iwidth,3,string)
         if(j.eq.1.and.k.eq.i)call legendline(-.2,1.1,0,' !Ay!@='&
              &//string(1:iwidth))
         call color(15);call dashset(0)
         call jdrwstr(wx2nx(x(npts)),wy2ny(.85*denmax)&
              &,'[!Bn!di!@!d(!Bx!@)-1]/!Ay!@!@',-1.2)
         call legendline(.4,.08,258,'!Bx!@')
      enddo
   enddo
enddo
call pltend
call multiframe(0,0,0)
write(*,*)
phimax=phimxin
vs(1)=vs1in
vs(2)=vs2in
vt(1)=vt1in

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
  if(lfover)then ! Overplot on frame 1 the f of frame 3.
     call setframe(2) ! Ensure the autoinit does not clear frames.
     call autoinit(vofv,fofv,nofv) ! Set scaling
     call setframe(0) ! Set to plot on frame 1.
     call dashset(4)
     call polyline(vofv,fofv,nofv)
     call dashset(0)
  endif
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
!  call autoplot(vha,Te*2.*sinh(dphivh/(2.*Te)),nvh)
  call autoinit(vha,Te*2.*sinh(dphivh/(2.*Te)),nvh)
  call axis
  call axis2
  call axlabels('','Forces(v!dh!d)')
  call legendline(.6,-.1,258,'v!dh!d')
  call legendline(.9,.9,258,'(b)')
  call winset(.true.)
  call polyline([vha(1),vha(nvh)],[0.,0.],2)
  if(index.ne.0)call polyline([vh0,vh0],[-10.,10.],2)
  call polyline(vha,Forcevh,nvh)      ! Total force
  call legendline(.05,.7,0,'F!di!d+F!de!d')
  call color(1);call dashset(1)
  call polyline(vha,Te*2.*sinh(dphivh/(2.*Te)),nvh)
  call legendline(.05,.8,0,'F!de!d')
  call color(2);call dashset(2)
  call polyline(vha,-denave*dphivh/Te+Forcevh,nvh)   ! Ion force
  call legendline(.05,.9,0,'F!di!d')
  call dashset(0)
  if(ltesting)then
     call winset(.false.)
     call color(4)
     call polyline(vha,analFrci,nvh)  ! Analytic total force.
     call legendline(.95,.8,0,'analFrci')
     call color(5)
     call polyline(vha,Te*analdendift,nvh)  ! Analytic Density difference.
     call legendline(.95,.6,0,'analdendift')
     call color(6)
     call polyline(vha,Te*directddift,nvh)  ! Direct density difference
     call legendline(.95,.4,0,'directdendift')
  endif
  call pltend
  call color(15)
  call multiframe(0,0,0)
  if(ltesting)then
     call autoplot(vha,analdendiff,nvh)
     call axlabels('vh','delta-n intrinsic')
     call color(1)
     call polyline(vha,directddiff,nvh)
     call pltend
  endif
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
  call color(4)
  if(ltesting.and.ns.eq.1.)call polyline(vha,analFrci/phimax**2,nvh)
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
  call axlabels('v!dh!d  or  v!d!A;!@!d','[!Ad!@F/!Ad!@x(v!dh!d)]/!Ay!@!u2!u')  
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
subroutine analcomp(i)
! Compare the densities and forces between numerical and analytical forms.
!  cfp1=2.*(1.+abs(delphi)*alog(2.*sqrt(2.*phimax)/sqrt(abs(delphi))))
!  anddenf=-cfp1*f1prime(vh)*2.*phimax-(1./6)*f3prime(vh)*(2.*phimax)**2
! The above gloss makes no significant difference.  
  anddenf=-2.*f1prime(vh)*2.*phimax-(1./6)*f3prime(vh)*(2.*phimax)**2
  andFrcf=-f1prime(vh)*(2.*phimax)**2-f3prime(vh)/9.*(2.*phimax)**3
  isigma=nint(sign(1.,delphi))
!  phi=abs(delphi)/2.
!  call denhill(nc,dc,vs-vh,vt,vh,phimax,phi,isigma,nofv,denacross,dendiff)
  denacross=denionhill(0.)
  isigma=-nint(sign(1.,delphi))
!  call denhill(nc,dc,vs-vh,vt,vh,phimax,phi,isigma,nofv,denup,dendiff)
  denup=denionhill(abs(delphi))
  directddiff(i)=isigma*(denup-denacross)
!  directddiff(i)=isigma*dendiff            ! Direct ion dendiff intrinsic
!  phi=-abs(delphi)/2.
!  call denhill(nc,dc,vs-vh,vt,vh,phimax,phi,isigma,nofv,dendown,dendiff)
  dendown=denionhill(0)
  directddift(i)=sign(1.,delphi)*(denacross-dendown)
  ddenext=sign(1.,delphi)*(denup-dendown)
  analdendift(i)=anddenf+ddenext           ! Analytic Total ni difference.
  analdendiff(i)=anddenf           ! Analytic intrinsic ni difference.
  analFrci(i)=andFrcf!+denave*(-delphi+Te*2*sinh(delphi/(2.*Te)))!Total frc.
!  write(*,'(7f10.5)')vh,anddenf,delphi,andFrcf,delphi,denacross-dendown&
!       ,(denup-dendown)/delphi
!  write(*,'(6f10.5)')vh,anddenf,delphi-ddenext,andFrcf,diffint
end subroutine analcomp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testdendiff
  ! Do a direct comparison of the anal dendiff and numerical.
  delphi=0.01   ! adjust by hand
  do i=1,nvh
     vh=vhmin+(vhmax-vhmin)*(i-1.)/(nvh-1.)
     vha(i)=vh
     analdendiff(i)=-2.*f1prime(vh)*2.*phimax-(1./6)*f3prime(vh)*(2.*phimax)**2
     isigma=nint(sign(1.,delphi))
     denplus=denionhill(0.)
     isigma=-isigma
     denminus=denionhill(abs(delphi))
     directddiff(i)=-isigma*(denplus-denminus)
  enddo
  call autoplot(vha,analdendiff,nvh)
  call color(4)
  call polyline(vha,directddiff,nvh)
  call pltend
  call exit
end subroutine testdendiff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotanalcomp
   call multiframe(2,1,3)
   call dcharsize(.02,.02)
   phimaxin=phimax
   do i=1,npsi
      phimax=phimaxin*(float(i)/npsi)**2
      call scanvh(index)
      if(i.eq.1.)then
         call pltinit(vha(1),vha(nvh),0.,0.4)
         call axis
         call axis2
         call legendline(.9,.9,258,'(b)')
         call legendline(-.12,.65,258,'!Bf!@!di!A;!@!d(!Bv!@)')
         call legendline(.6,-.1,258,'!Bv!@')
         call polyline(vha,fofvh,nvh)
         call polymark(vha(index),fofvh(index),1,10)
!         call autoplot(x,phiofx/phimax,npts)
         call pltinit(x(1),x(npts),-.1,1.05)
         call legendline(.9,.9,258,'(c)')
         call axis; call axis2
         call dashset(1)
         call polyline(x,phiofx/phimax,npts)
         call fwrite(phimax,iwidth,3,string)
         call legendline(0.05,.9,0,' !Ay!@='//string(1:iwidth))
         call dashset(0)
         call legendline(-.1,.65,258,'!Af!@/!Ay!@')
         call legendline(.6,-.1,258,'!Bx!@')
      elseif(i.eq.npsi)then
         call polyline(x,phiofx/phimax,npts)
         call fwrite(phimax,iwidth,3,string)
         call legendline(0.05,.75,0,' !Ay!@='//string(1:iwidth))
         call dashset(4)
         call polyline([x(1),x(npts)],[0.,0.],2)
         call dashset(0)
         call multiframe(0,0,0)
         call pltend
      endif
! Find dn/dphi      
      isigma=-nint(sign(1.,delphi0))
      dphi=1.*delphi0
      call denhill(nc,dc,vs,vt,vh,phimax,-dphi/2.,isigma,nofv,theden,dendiff)
      call denhill(nc,dc,vs,vt,vh,phimax, dphi/2.,isigma,nofv,othden,dendiff)
      dnidphi=dendiff/dphi
      dnidphi=.3
!      write(*,*)'dnidphi,Te=',dnidphi,Te
      psia(i)=phimax
      delpha(i)=delphi0
      f3p=f3prime(vh0)
      sigmaD=sign(1.,f3p)
      Ts=1/(1/Te-dnidphi)
! Third version
!      Ts=1.3*Te   ! Additive dnidphi is better.
      g=sigmaD*Ts*(2./9.)*f3p*phimax
      analdp(i)=sigmaD*phimax*(1+g)*(1.-sqrt(1-2.*g/(1+g)))
!      analdp(i)=sigmaD*phimax*g ! This is practically just as good.
      write(*,'(a,6f10.6)')'psi,vh0,dphi,anal',phimax,vh0,delphi0&
           ,analdp(i)
   enddo
!   call autoplot(psia,delpha,npsi)
   call charsize(.02,.02)
   call lautoplot(psia,abs(delpha),npsi,.true.,.true.)
   call legendline(.1,.9,258,'(a)')
   call axlabels('!Ay!@','!ADf!@')
   call axis2
   call color(4)
   call dashset(1)
   call polyline(psia,abs(analdp),npsi)
   call dashset(0)
   call pltend
   call exit
end subroutine plotanalcomp
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
! phi is the potential at which to find the distribution/density.
! isigma is the side of the hill on which this potential lies.
! All velocities are relative to the hill's rest frame.
! On exit:
! theden is the integral of fv dv. I.e. the density
! dendiff is the difference between this theden and the last.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine denhill(nc,dc,vs,vt,vh,phimax,phi,isigma,nofv,theden,dendiff)
  implicit double precision (a-h,o-z)
  real :: dc(nc),vs(nc),vt(nc),phimax,phi,theden,vh,dendiff
  integer :: isigma,nofv
  double precision :: dtheden=0.
  one=1.
  if(isigma.ne.1.and.isigma.ne.-1)stop 'isigma must be +-1'
  phix2=2.*phi
  phimx2=2.*phimax
  vthresh=isigma*sqrt(max(0.,phimx2-phix2)) ! prevent rounding NANs.
  vinf=vthresh ! Silence warning
  vmax=4.4+max(abs(maxval(vs)-vh),abs(minval(vs)-vh))
  sqtwopi=sqrt(2.*3.1415926)
  nvhere=nofv+1-mod(nofv,2) ! Odd number of points so thresh is at (nvhere+1)/2
  ithresh=(nvhere+1)/2
  ddenlast=dtheden
  dtheden=0.
  fim=0.
  vm=(ithresh-1)*(-vmax-vthresh)/(ithresh-1.)+vthresh
  do i=1,nvhere
     if(i.le.ithresh)then  ! Lower and upper velocity ranges
        vi=(ithresh-i)*(-vmax-vthresh)/(ithresh-1.)+vthresh
     else                  ! vthresh to +-vmax
        vi=(i-ithresh)*(vmax-vthresh)/(nvhere-ithresh)+vthresh
     endif
     vsign=sign(one,vi)
     enx2=phix2+vi**2
     if(i.eq.ithresh)then  ! Crossing threshold. Use prior sign this step
        vinf=sign(sqrt(max(0.,enx2)),vinf)
     elseif(nint(vsign).ne.isigma)then        ! Moving inward
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
     dtheden=dtheden+ (fim+fi)*0.5*(vi-vm)
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
  dendiff=real(dtheden-ddenlast)
  theden=real(dtheden)
end subroutine denhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the ion density at a potential phiminf relative to phiinf on
! the side isigma. To be called by BGKint.
real function denionhill(phiminf)
  use AsymHill
  phimzero=phiminf+isigma*delphi/2.
  if(ldenalt)then   ! This density must be consistent with delphi.
     call denhill(nc,dc,vs-vh,vt,vh,phimax,phimzero,isigma,nofv&
          &,denionhill,dendiff)
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

!call testdendiff

! delphi scale analysis rarely used:
if(ldscale.and..not.ltestnofx)call plotanalcomp

! Plots of non-equilibrium case -f
if(ltestnofx)then
   call finddelphi
   call plotfvhill
!   write(*,'(a,9f6.3)')'plotfvhill: vshift,vthm,dens',(vs(i),vt(i),dc(i),i=1,nc)
   call finddenofx
   if(lexplain)call plotexplain
   call plotdenofx
endif

! The workhorse scanning of ns vshifts (-s)
if(ns.gt.0)then
call scanspacing
call forceplotns
endif

if(lfvden)then  ! Multi-T,n scan.
call plotfvden
endif

end program Asym
