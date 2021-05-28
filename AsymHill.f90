module AsymHill
  integer, parameter :: npts=200,nofv=200,nvh=20,nc=2
  real, dimension(npts) :: x,phiofx,denofx,deninteg,dFdelx
  real, dimension(nofv) :: vofv,fofv,Pfofv,Forceofv
  real, dimension(nvh) :: vha,Forcevh,delFdxvh
  real, dimension(nc) :: vs=[1.5,-1.5],vt=1.,dc=[1.,.0]
  character*10 string
  real :: vh=0.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testfvhill(phimax,delphi,phi)
  call multiframe(4,1,1)
  call dcharsize(.02,.02)
  do isigma=-1,1,2
     call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phi,isigma,&
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
     call axlabels('v-v!dh!d','P(v)=!AJ!@f(v)dv')
  enddo
  call multiframe(0,0,0)
  call pltend
end subroutine testfvhill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddenofx(phimax,delphi)
! Integrate over the potential hill to find the density, force, and
! delF/delx (arrays of x), for the current vh and other parameters. 
  xmax=10.
  deninteg=0.
  dFdelx=0.
  do i=1,npts
     x(i)=-xmax+2.*xmax*(i-1.)/(npts-1.)
     isigma=int(sign(1.,x(i)))
     dph=isigma*delphi/2.
     xcor=x(i)/4.
!     xcor=x(i)/4./sqrt(1.-dph/phimax)  !Symmetrizes curvature at origin.
     phiofx(i)=(phimax-dph)/cosh(xcor)**4+dph
     call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiofx(i),isigma, &
          nofv,vofv,fofv,Pfofv)
     denofx(i)=Pfofv(nofv)
     if(i.gt.1)deninteg(i)=deninteg(i-1)+&      !\int n dphi
          (denofx(i)+denofx(i-1))*(phiofx(i)-phiofx(i-1))*.5
     if(i.gt.1)dFdelx(i)=dFdelx(i-1)+&          !\int n (dphi/dx)dphi
          (denofx(i)+denofx(i-1))*(phiofx(i)-phiofx(i-1))**2*.5 &
          /(x(i)-x(i-1))
  enddo
  write(*,'(a,f8.5,a,f9.5,a,f9.5)')'vh=',vh,' Ion force=',-deninteg(npts),' delF/delx=',dFdelx(npts)
end subroutine finddenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotdenofx
  denave=(denofx(1)+denofx(npts))/2
  call multiframe(3,1,0)
  call autoplot(x,phiofx,npts)
  call axis2
  call polyline([x(1),x(npts)],[0.,0.],2)
  call axlabels('','!Af!@(x)')
  call autoplot(x,denofx,npts)
  call axis2
  call polyline([x(1),x(npts)],[denave,denave],2)
  call axlabels('','n!di!d(x)')
  call winset(.true.)
  call color(4)
  call polyline(x,denave*exp(phiofx),npts)
  call color(15)
  call autoplot(x,deninteg,npts)
  call axis2
  call axlabels('x','!AJ!@n!di!dd!Af!@')
  call polyline([x(1),x(npts)],[0.,0.],2)
  call pltend
  call multiframe(0,0,0)
end subroutine plotdenofx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finddelphi(phimax,delphi)
  real, dimension(2) :: deninf
! Iterate to find delphi consistent with specified ion distributions.
!  write(*,*)'j   delphi      n(-)      n(+)'
  do j=1,4
     do i=1,2
        isigma=-3+2*i
        phiinf=isigma*delphi/2.
!     write(*,*)phimax,delphi,phi
        call fvhill(nc,dc,vs-vh,vt,phimax,delphi,phiinf,isigma,&
             nofv,vofv,fofv,Pfofv)
        deninf(i)=Pfofv(nofv)
     enddo
     write(*,'(i2,5f10.5)')j,delphi,deninf
     delphi=alog(deninf(2)/deninf(1))    ! Times Te
!     delphi=delphi/3.
  enddo
!delphi=0.
end subroutine finddelphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scanvh(vhmin,vhmax,index,phimax,delphi)
! Scan vhmin to vhmax to construct arrays of F and delF/delx
! Return the highest index of vh at which F crosses from + to -.
  Fprior=0.
  do i=1,nvh
     vh=vhmin+(i-1.)*(vhmax-vhmin)/(nvh-1.)
     call finddelphi(phimax,delphi)
     call finddenofx(phimax,delphi)
     Force=-deninteg(npts)
     if(Fprior.gt.0.and.Force.le.0.)index=i
     vha(i)=vh
     Forcevh(i)=Force
     delFdxvh(i)=dFdelx(npts)
!     write(*,*)i,vha(i),Forcevh(i),delFdxvh(i)
     Fprior=Force
  enddo
end subroutine scanvh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initvofv
! Initialize velocity range using current vh as maximal.
  vmax=4.5+max(abs(maxval(vs)-vh),abs(minval(vs)-vh))
  do i=1,nofv
     vofv(i)=-vmax+2.*vmax*(i-1.)/(nofv-1.)
  enddo
end subroutine initvofv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timenofx
! This routine takes 8 seconds for 1M densities, when npts=400.
  if(nc.eq.2)dc(2)=.5
  if(nc.eq.2)vs(nc)=-vs(1)
  phimax=.2
  xmax=10.
  ncycle=2500
  call initvofv
  write(*,*)'Timing routine starting'
  do j=1,ncycle
  do i=1,npts
     x(i)=-xmax+2.*xmax*(i-1.)/(npts-1.)
     phiofx(i)=phimax/cosh(x(i)/4.)**4
     isigma=int(sign(1.,x(i)))
     call fvhill(nc,dc,vs,vt-vh,phimax,0.,phiofx(i),isigma, &
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
subroutine fvhill(nc,dc,vs,vt,phimax,delphi,     &
     phi,isigma,nofv,vofv,fofv,Pfofv)
  real :: dc(nc),vs(nc),vt(nc),phimax,phi,delphi
  integer :: isigma, nofv
  real, dimension(nofv) :: vofv, fofv, Pfofv

  if(isigma.ne.1.and.isigma.ne.-1)stop 'isigma must be +-1'
  phix2=2.*(phi    -isigma*delphi/2.)  ! 2*(phi   -phiinf)
  phimx2=2.*(phimax-isigma*delphi/2.)  ! 2*(phimax-phiinf)
  if(phimax.lt.phi)then
     write(*,*)'Normally phi must be < phimax',phi,phimax
     vthresh=0.
  elseif(phimx2.lt.0)then
     write(*,*)'No Potential peak this side',isigma,phimx2
     vthresh=0.
  elseif(phix2.lt.0)then
     write(*,*)'Impossible phi (<phiinf)',phi,isigma,delphi
     stop
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
           vinf=vsign*sqrt(phix2+vofv(i)**2)
           fi=gcoef*exp(-(vinf-vs(j))**2/vt2x2)   ! f(v)=finf(vinf)
        elseif((phix2+vofv(i)**2).gt.phimx2)then ! Passing
           vinf=vsign*sqrt(phix2+vofv(i)**2+2.*isigma*delphi)
           fi=gcoef*exp(-(vinf-vs(j))**2/vt2x2)   ! f(v)=finf(vinf)
        else   ! Reflected
           vinf=vsign*sqrt(phix2+vofv(i)**2)
           fi=gcoef*exp(-(-vinf-vs(j))**2/vt2x2)   ! f(v)=finf(vinf)
        endif
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
! Optional main for testing. 
use AsymHill
phimax=0.3
phi=.1
delphi=0.0
vh=1.
call initvofv
call finddelphi(phimax,delphi)
call testfvhill(phimax,delphi,phi)
write(*,'(a,8f8.3)')'testnofx: vshift,vtherm,dens',vs,vt,dc
call finddenofx(phimax,delphi)
call plotdenofx
vhmin=-5.;vhmax=+5
index=2
call scanvh(vhmin,vhmax,index,phimax,delphi)
write(*,*)index,vha(index-1),vha(index+1),Forcevh(index-1),Forcevh(index)

call autoplot(vha,Forcevh,nvh)
call axis2
call pltend
!call timenofx ! For timing comment the above.
end program
