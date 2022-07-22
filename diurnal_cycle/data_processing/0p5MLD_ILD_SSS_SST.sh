#!/bin/sh
set -xe

# Taken from /lustre/f2/dev/ncep/JieShun.Zhu/CFSm501dy/CFSm501dy1980010200/0DIA/FILE/0p5MLD_ILD_SSS_SST.sh 
#---------------------------------------------
#source /opt/cray/pe/modules/3.2.10.5/init/sh
#README before running it on gaea, run the following commands for loading proper inter compiler and netcdf.
####module purge
#module swap PrgEnv-gnu/6.0.3 PrgEnv-intel/6.0.3
#module load cray-hdf5/1.10.2.0
#module swap cray-netcdf/4.4.0 cray-netcdf/4.6.1.3
#---------------------------------------------

iyrb=2004
iyre=2004
disot=0.5

cat >MLD_ILDp.f<<fEOF
      Program main
      INCLUDE 'netcdf.inc'

      integer im,jm,km
      parameter (im=720,jm=410,km=50)
      parameter (ndtc=7)
      real t(im,jm,km),s(im,jm,km)
      real mld(im,jm),ild(im,jm),sss(im,jm)
      real zt(km)
      data zt/0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5, 
     &       12.5,17.5,25.,35.,45.,55.,65.,75.,85.,95.,105., 
     &       115.,125.,135.,145.,155.,165.,175.,185.,195.,205., 
     &       215.,225.,238.4779,262.2945,303.0287, 
     &       366.7978, 459.091, 584.6193, 747.187, 
     &       949.5881, 1193.53, 1479.588, 1807.187, 
     &       2174.619, 2579.091, 3016.798, 3483.029, 
     &       3972.294, 4478.478/
      real dtc(ndtc)
      data dtc/2.5,5.,10.,15.,20.,25.,28./
      integer days(12)
      data days/31,28,31,30,31,30,31,31,30,31,30,31/
      character*4 cyr
      character*2 cnm(0:31)
      data cnm/'00',
     &         '01','02','03','04','05','06','07','08','09','10',
     &         '11','12','13','14','15','16','17','18','19','20',
     &         '21','22','23','24','25','26','27','28','29','30',
     &         '31'/
      data undef/-1.0000000E+20/
      character*200 indir
C============
      do iyr=$iyrb,$iyre
        write(cyr,'(i4.4)')iyr
      do imo=1,12
        ic=0
        open(11,file='/lustre/f2/dev/ncep/Jack.Reeveseyre/'//
     &      'MLDILD/MLDILDSSSSST_hr.'//cyr 
     &      //cnm(imo), 
     &      form='unformatted',access='direct',recl=im*jm)

        nday=days(imo)
        if(imo.eq.2)then
          if(mod(iyr,4).eq.0)nday=29 
        endif 
      do id=1,nday
      do ih=0,23
        indir='/lustre/f2/dev/ncep/JieShun.Zhu/CFSm501hr/'//
     &        'CFSm501hr1980010200/DATA/ocn_'//cyr//'_'//cnm(imo)// 
     &        '_'//cnm(id)//'_'//cnm(ih)//'.nc'
        print*,trim(indir)
        status=nf_open(trim(indir),nf_nowrite,ncid)
        iret=NF_INQ_VARID(NCID,'temp',IDT)
        iret=NF_GET_VAR_REAL(NCID,IDT,t)
        iret=NF_INQ_VARID(NCID,'salt',IDS)
        iret=NF_GET_VAR_REAL(NCID,IDS,s)
        iret=NF_CLOSE(ncid)

!        print*,t(100,205,1),t(100,205,50),t(1,1,1) 
!        print*,s(100,205,1),s(100,205,50),s(1,1,1)
! mld (MLD)
          call mixed_layer(im,jm,km,t(:,:,:)+273.15,s(:,:,:),
     &                                           zt,mld,sss,undef)
! sfc isothm layer depth (ILD)
          call sfc_isothm_layer(im,jm,km,t(:,:,:)+273.15,
     &                                           zt,ild,undef)
          ic=ic+1
          write(11,rec=ic)mld
          ic=ic+1
          write(11,rec=ic)ild
          ic=ic+1
          write(11,rec=ic)sss
          ic=ic+1
          write(11,rec=ic)t(:,:,1)

      enddo !end ih
      enddo !end id 
        print*,ic
        print*,' '
      enddo !end im
        close(11)
      enddo !end iyr  
      print*,ild(:,100)-mld(:,100)
!      print*,ild(:,100)
      END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine isothm_layer(im,jm,km,dtc,temp,zlev,zisothm,undef)

      real, parameter :: c2k=273.15
      integer inumc,im,jm,km
      integer i,j,k
      real  dtc
      real, dimension(km) :: tz,zlev
      real, dimension(im,jm) :: zisothm
      real, dimension(im,jm,km) :: temp
      real  a,b,tc,undef

      tc=c2k+dtc

      do j=1,jm
      do i=1,im

         zisothm(i,j)=undef
         if (temp(i,j,1) .GE. tc) then
            do k=1,km
               tz(k)=temp(i,j,k)
            enddo
            do k=2,km
               if (tz(k) .LT. -3.0) go to 111
               if (tz(k) .LT. tc) then
                  a = (tz(k)-tc) / (tz(k)-tz(k-1))
                  b = (tc-tz(k-1)) / (tz(k)-tz(k-1))
                  zisothm(i,j)=a*zlev(k-1)+b*zlev(k)
                  go to 111
               endif
            enddo
         endif
 111     continue
      enddo
      enddo

      return
      end subroutine isothm_layer

      subroutine mixed_layer(im,jm,km,temp,salt,zlev,mld,sss,undef)

      real, parameter :: disot=$disot, c2k=273.15
      integer im,jm,km
      integer i,j,k,kmsk,kbm,kbp,krf
      real, dimension(km) :: zlev,plev,sa,ta,th,rho
      real, dimension(im,jm) :: mld,sss
      real, dimension(im,jm,km) :: temp,salt
      real  a,b,deltarho,dr,rb,undef

      do k=1,km
         plev(k) = press(zlev(k),980.0)
      enddo

      sss=undef
      do j=1,jm
      do i=1,im
         if(salt(i,j,1).gt.0.0)sss(i,j)=salt(i,j,1)

         kmsk = 0
         do k=1,km
            if (temp(i,j,k).GT.0.0 .AND. salt(i,j,k).GT.0.0) then
               ta(k) = temp(i,j,k)-c2k
               sa(k) = salt(i,j,k)
               kmsk = k
            endif
         enddo

         if (kmsk.EQ.0 .OR. ta(1).LT.-3.0) then
            mld(i,j)=undef
         else
            deltarho = (density(0.0,ta(1)-disot,sa(1))
     &                - density(0.0,ta(1),sa(1)))
            do k=1,kmsk
               th(k) = theta(plev(k),ta(k),sa(k),0.0)
               rho(k) = density(0.0,th(k),sa(k)) - 1000.0
            enddo
            krf = 1
            kbm = 0
            kbp = 0
            do k = krf,kmsk
               if ((rho(k)-rho(krf)) .GE. deltarho) then
                   kbp = k
                    exit
               endif
            enddo
            if (kbp .LE. 1) then
               mld(i,j) = undef
            else
               kbm = kbp - 1
               rb = rho(krf) + deltarho
               dr = rho(kbp) - rho(kbm)
               a = (rho(kbp) - rb) / dr
               b = (rb - rho(kbm)) / dr
               mld(i,j) = zlev(kbm)*a + zlev(kbp)*b
            endif
         endif
      enddo
      enddo

      end subroutine mixed_layer

      subroutine sfc_isothm_layer(im,jm,km,temp,zlev,sitd,undef)

      real, parameter :: disot=$disot
      integer im,jm,km
      integer i,j,k
      real, dimension(im,jm) :: sitd
      real, dimension(im,jm,km) :: temp
      real, dimension(km) :: zlev,tz
      real  a,b,tc,undef

      do j=1,jm
      do i=1,im

         sitd(i,j)=undef

         if (temp(i,j,1).GE.0.0) then
            tc=temp(i,j,1)-disot
            do k=1,km
               tz(k)=temp(i,j,k)
            enddo
            do k=2,km
               if (tz(k).LT.0.0) go to 112
               if (tz(k).LT.tc) then
                  a = (tz(k)-tc) / (tz(k)-tz(k-1))
                  b = (tc-tz(k-1)) / (tz(k)-tz(k-1))
                  sitd(i,j) = a*zlev(k-1) + b*zlev(k)
                  go to 112
               endif
            enddo
         endif

 112     continue

      enddo
      enddo

      return
      end subroutine sfc_isothm_layer

      subroutine ocean_heat(im,jm,km,temp,salt,zblev,zlev,ocnhc,undef)

      integer, parameter :: kmh=26
      real, parameter :: c2k=273.15
      integer im,jm,km
      integer i,j,k
      real, dimension(km) :: zblev,zlev,plev
      real, dimension(im,jm) :: ocnhc
      real, dimension(im,jm,km) :: salt,temp
      real  dptlyr,rk,sk,tk,undef
      real  pk,rhm,rhp,tempk,zk

      k=kmh
      zk=0.5*(300.0+zblev(k-1))
      pk=press(zk,980.0)
      rhm = (zk-zlev(k-1))/(zlev(k)-zlev(k-1))
      rhp = (zlev(k)-zk)/(zlev(k)-zlev(k-1))

      do k=1,km
         plev(k) = press(zlev(k),980.0)
      enddo

      do j=1,jm
      do i=1,im

         ocnhc(i,j)=undef

         if (temp(i,j,kmh).GT.0.0 .AND. salt(i,j,kmh).GT.0.0) then
            ocnhc(i,j)=0.
            do k=1,kmh-1
               tk=temp(i,j,k)-c2k
               sk=salt(i,j,k)
               rk=density(plev(k),tk,sk)
               if (k .eq. 1) then
                  dptlyr=zblev(k)
               else
                  dptlyr=zblev(k)-zblev(k-1)
               endif
               ocnhc(i,j)=ocnhc(i,j) + rk*temp(i,j,k)*dptlyr
            enddo
            k=kmh
            tempk=rhp*temp(i,j,k-1) + rhm*temp(i,j,k)
            tk=tempk - c2k
            sk=rhp*salt(i,j,k-1) + rhm*salt(i,j,k)
            rk=density(pk,tk,sk)
            dptlyr=300.0-zblev(k-1)
            ocnhc(i,j)=ocnhc(i,j)+rk*tempk*dptlyr
            ocnhc(i,j)=ocnhc(i,j)*3996.
         endif
      enddo
      enddo

      return
      end subroutine ocean_heat

      subroutine tchp26(im,jm,km,temp,salt,zblev,zlev,ocnhcp,undef)

      real, parameter :: c2k=273.15, t26=26.0
      integer im,jm,km
      integer i,j,k,k26
      real, dimension(km) :: zblev,zlev,plev,tz
      real, dimension(im,jm) :: ocnhcp,z26isothm
      real, dimension(im,jm,km) :: salt,temp
      real  dptlyr,rk,sk,tk,undef
      real  rhm,rhp,pk,zk
      real  a,b,skk,skm,tc,z26
      logical*1 lbms(im,jm)

      tc=c2k+t26

      do k=1,km
         plev(k) = press(zlev(k),980.0)
      enddo

      do j=1,jm
      do i=1,im

         z26isothm(i,j)=undef
         lbms(i,j)=.false.
         if (temp(i,j,1) .GE. tc) then
            do k=1,km
               tz(k)=temp(i,j,k)
            enddo
            k = 1
            do while (tz(k).GE.tc)
               k26 = k
               if (k.EQ.km) exit
               k = k + 1
            enddo
            k = k26
            if (tz(k) .GT. tc) then
               if (k.LT.km .AND. tz(k+1).GT.0.0) then
                  k = k + 1
                  a = (tz(k)-tc) / (tz(k)-tz(k-1))
                  b = (tc-tz(k-1)) / (tz(k)-tz(k-1))
                  z26isothm(i,j) = a*zlev(k-1) + b*zlev(k)
                  lbms(i,j)=.true.
               else if (k.GE.2 .AND. tz(k).LT.tz(k-1)) then
                  a = (tz(k)-tc) / (tz(k)-tz(k-1))
                  b = (tc-tz(k-1)) / (tz(k)-tz(k-1))
                  z26 = a*zlev(k-1) + b*zlev(k)
                  if (z26.LE.zblev(k)) then
                     z26isothm(i,j) = z26
                     lbms(i,j)=.true.
                  endif
              endif
            else if (tz(k).EQ.tc) then
               z26isothm(i,j) = zlev(k)
               lbms(i,j)=.true.
            endif
         endif

      enddo
      enddo

!
!---------- get ocean heat potential relative to 26C (TCHP) ------------
!
      do j=1,jm
      do i=1,im

         if (temp(i,j,1) .GT. 0.0) then
            ocnhcp(i,j)=0.0
         else
            ocnhcp(i,j)=undef
            cycle
         endif

         if (lbms(i,j)) then   ! we have water above 26c

            z26 = z26isothm(i,j)
!
!  case where Z26 is within the topmost layer
!
            if (z26 .LE. zblev(1)) then
               tk=temp(i,j,1)-c2k
               if (salt(i,j,1) .GT. 0.0) then
                   sk=salt(i,j,1)
               else
                   sk=35.   ! fake salinity
               endif
               rk=density(plev(1),tk,sk)
               dptlyr=z26
               ocnhcp(i,j) = rk*(tk-t26)*dptlyr*3996.
!
!  case where z26 is below the top layer and above the bottom
!
            else
               k26 = 1
               do k=2,km
                  if (z26.GT.zblev(k-1) .AND.  z26.LE.zblev(k)) k26=k
               enddo

               ocnhcp(i,j)=0.0
               do k=1,k26-1
                  tk=temp(i,j,k)-c2k
                  if (salt(i,j,K) .GT. 0.0) then
                     sk=salt(i,j,k)
                  else
                     sk=35.   ! fake salinity
                  endif
                  rk=density(plev(k),tk,sk)
                  if (k .EQ. 1) then
                     dptlyr=zblev(1)
                  else
                     dptlyr=zblev(k)-zblev(k-1)
                  endif
                  ocnhcp(i,j)=ocnhcp(i,j)+rk*(tk-26.0)*dptlyr
               enddo
               k=k26
               zk=0.5*(z26+zblev(k-1))
               pk=press(zk,980.0)
               rhm = (zk-zlev(k-1))/(zlev(k)-zlev(k-1))
               rhp = (zlev(k)-zk)/(zlev(k)-zlev(k-1))
               tk=rhp*temp(i,j,k-1) + rhm*temp(i,j,k) - c2k
               if (salt(i,j,k-1) .GT. 0.0) then
                  skm=salt(i,j,k-1)
               else
                  skm=35.   ! fake salinity
               endif
               if (salt(i,j,k) .GT. 0.0) then
                  skk=salt(i,j,k)
               else
                  skk=35.   ! fake salinity
               endif
               sk=(rhp*skm + rhm*skk)
               rk=density(pk,tk,sk)
               dptlyr=z26-zblev(k-1)
               ocnhcp(i,j)=ocnhcp(i,j)+rk*(tk-26.0)*dptlyr
               ocnhcp(i,j)=ocnhcp(i,j)*3996.
            endif
!
!  case where temperature is above 26C down to the bottom
!
         else if ((temp(i,j,1)-c2k) .GT. t26) then
            ocnhcp(i,j)=0.0
            do k=1,km
               if (temp(i,j,k) .GT. undef) then
                  tk=temp(i,j,k)-c2k
                  if (salt(i,j,k) .GT. 0.0) then
                     sk=salt(i,j,k)
                  else
                     sk=35.   ! fake salinity
                  endif
                  rk=density(plev(k),tk,sk)
                  if (k .EQ. 1) then
                     dptlyr=zblev(1)
                  else
                     dptlyr=zblev(k)-zblev(k-1)
                  endif
                  ocnhcp(i,j)=ocnhcp(i,j)+rk*(tk-26.0)*dptlyr
               endif
            enddo
            ocnhcp(i,j)=ocnhcp(i,j)*3996.
         endif

      enddo
      enddo

      return
      end subroutine tchp26
 
      function press(z, g)

!   copy from cfs_ocean_time.f and modified
!   depth (z) in meters and grav acc'l (g) in cm/sec**2

      integer, parameter :: itr=20
      integer i
      real p, a0, z, g, press
      real(kind=8) :: e, ae, es
!
      p = z*(1.0076+z*(2.3487e-6-z*1.2887e-11))
      e = zeta(p,g)-z
      ae = abs(e)
      es = ae*2.
      do i = 1,itr
        a0 = 0.972643+p*(1.32696e-5-p*(6.228e-12+p*1.885e-16))
        a0 = a0/(1.0+1.83e-5*p)
        p = p-((g+1.113e-4*p)/a0)*e*0.001
        es = ae
        e = zeta(p,g)-z
        ae = abs(e)
        if (ae .le. 0.01) exit
      enddo
!
      press = p
!
      end function press

      function zeta(p, glat)
!
!   copy from cfs_ocean_time.f and modified

      real p, glat, z, zeta

      z = ((-3.434e-12*p+1.113e-7)*p+0.712953)*p+
     &    14190.7*log(1.0+1.83e-5*p)
      z = (z/(glat+1.113e-4*p))*1000.

      zeta = z
!
      end function zeta

      function density(prs, tmp, sal)
!
!   copy from cfs_ocean_time.f and modified
!     Density is in units of kg/m**3  (1 g/cm**3 = 1000 kg/m**3)

      real density, prs, tmp, sal
      real p, t, s, kstp, k0, kw, d0, dw
!
      s = sal
      t = tmp
      p = prs/10.00
!
      kw = 19652.21+(148.4206-(2.327105-(1.360477e-2-
     &              5.155288e-5*t)*t)*t)*t
!
      k0 = kw+s*(54.6746-(0.603459-(1.09987e-2-6.1670e-5*t)*t)*t)  
     &        +sqrt(s*s*s)*(7.944e-2+(1.6483e-2-5.3009e-4*t)*t)
!
      kstp = k0+p*((3.239908+(1.43713e-3+(1.16092e-4-5.77905e-7*t)*t)*t)
     &        +s*(2.2838e-3-(1.0981e-5+1.6078e-6*t)*t)                 
     &        +sqrt(s*s*s)*1.91075e-4                                 
     &        +p*((8.50935e-5-(6.12293e-6-5.2787e-8*t)*t)            
     &        -s*(9.9348e-7-(2.0816e-8+9.1697e-10*t)*t))) 
!
      dw = 999.842594+(6.793952e-2-(9.095290e-3-(1.001685e-4 
     &        -(1.120083e-6-6.536332e-9*t)*t)*t)*t)*t
!
      d0 = dw+s*(0.824493-(4.0899e-3-(7.6438e-5-(8.2467e-7       
     &        -5.3875e-9*t)*t)*t)*t)                            
     &        -sqrt(s*s*s)*(5.72466e-3-(1.0227e-4-1.6546e-6*t)*t)
     &        +s*s*4.8314e-4
!
      density = d0/(1.0-p/kstp)

      end function density

      function theta(p, t, s, pref)

      real(kind=8), parameter :: sqrt2 = 0.7071067811865475
      real theta, p,t, s, pref
      real del_p, del_t1, del_t2, del_t3, del_t4, tp, th

      del_p = pref-p
      del_t1 = del_p*atg(p,t,s)
      tp = t+0.5*del_t1
      del_t2 = del_p*atg((p+0.5*del_p),tp,s)
      tp = t+(-0.5+sqrt2)*del_t1+(1.0-sqrt2)*del_t2
      del_t3 = del_p*atg((p+0.5*del_p),tp,s)
      tp = t-sqrt2*del_t2+(1.0+sqrt2)*del_t3
      del_t4 = del_p*atg(pref,tp,s)
      th = (t+(del_t1+(1.0-sqrt2)*del_t2*2.0 
     &   + (1.0+sqrt2)*del_t3*2.0+del_t4)/6.0)
      theta = th

      end function theta

      function atg(p, t, s)

      real atg, p, t, s, ds, a

      ds = s-35.0
      a = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p       
     &        +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t  
     &        +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p      
     &        +(-4.2393e-8*t+1.8932e-6)*ds                  
     &        +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

      atg = a

      end function atg
fEOF
ftn -c -I/opt/cray/pe/netcdf/4.4.0/INTEL/15.0/include MLD_ILDp.f
ftn -o MLD_ILDp MLD_ILDp.o -L/opt/cray/pe/netcdf/4.4.0/INTEL/15.0/lib -lnetcdf
./MLD_ILDp >& out&
rm MLD_ILDp.f MLD_ILDp.o MLD_ILDp

cat >./MLDILDSSSSST_hr.ctl<<fEOF
dset ^MLDILDSSSSST_hr.%y4%m2
undef -1.0000000E+20
options template
title MLDILDSSS1x1
xdef 720 linear -279.75 0.5
ydef 410 levels
    -80.75 -80.25 -79.75 -79.25 -78.75 -78.25 -77.75 -77.25 
    -76.75 -76.25 -75.75 -75.25 -74.75 -74.25 -73.75 -73.25 -72.75 
    -72.25 -71.75 -71.25 -70.75 -70.25 -69.75 -69.25 -68.75 -68.25 
    -67.75 -67.25 -66.75 -66.25 -65.75 -65.25 -64.75 -64.25 -63.75 
    -63.25 -62.75 -62.25 -61.75 -61.25 -60.75 -60.25 -59.75 -59.25 
    -58.75 -58.25 -57.75 -57.25 -56.75 -56.25 -55.75 -55.25 -54.75 
    -54.25 -53.75 -53.25 -52.75 -52.25 -51.75 -51.25 -50.75 -50.25 
    -49.75 -49.25 -48.75 -48.25 -47.75 -47.25 -46.75 -46.25 -45.75 
    -45.25 -44.75 -44.25 -43.75 -43.25 -42.75 -42.25 -41.75 -41.25 
    -40.75 -40.25 -39.75 -39.25 -38.75 -38.25 -37.75 -37.25 -36.75 
    -36.25 -35.75 -35.25 -34.75 -34.25 -33.75 -33.25 -32.75 -32.25 
    -31.75 -31.25 -30.75 -30.2500991821289 -29.7508850097656 
    -29.2525539398193 -28.755687713623 -28.26047706604 -27.7674942016602 
    -27.2769145965576 -26.7892990112305 -26.3048038482666 
    -25.8239707946777 -25.3469352722168 -24.8742179870605 
    -24.4059276580811 -23.9425525665283 -23.4841766357422 
    -23.0312595367432 -22.5838489532471 -22.1423664093018 
    -21.7068290710449 -21.2776222229004 -20.8547191619873 
    -20.4384670257568 -20.0287990570068 -19.6260223388672 
    -19.2300262451172 -18.8410739898682 -18.4590148925781 
    -18.0840625762939 -17.7160263061523 -17.3550758361816 
    -17.0009727478027 -16.6538486480713 -16.3134174346924 
    -15.9797677993774 -15.6525726318359 -15.3318786621094 
    -15.0173177719116 -14.7088975906372 -14.406210899353 -14.109227180481 
    -13.8175039291382 -13.5309762954712 -13.2491683959961 
    -12.9719848632812 -12.6989212036133 -12.4298524856567 
    -12.164249420166 -11.9019641876221 -11.64244556427 -11.3855266571045 
    -11.1306390762329 -10.8776025772095 -10.6258354187012 
    -10.375147819519 -10.1249504089355 -9.87504959106445 
    -9.62495040893555 -9.37504959106445 -9.12495040893555 
    -8.87504959106445 -8.62495040893555 -8.37504959106445 
    -8.12495040893555 -7.87504911422729 -7.62495088577271 
    -7.37504911422729 -7.12495088577271 -6.87504911422729 
    -6.62495088577271 -6.37504911422729 -6.12495088577271 
    -5.87504911422729 -5.62495088577271 -5.37504911422729 
    -5.12495088577271 -4.87504911422729 -4.62495088577271 
    -4.37504911422729 -4.12495088577271 -3.87504911422729 
    -3.62495088577271 -3.37504911422729 -3.12495088577271 
    -2.87504911422729 -2.62495088577271 -2.37504911422729 
    -2.12495088577271 -1.87504923343658 -1.62495076656342 
    -1.37504923343658 -1.12495076656342 -0.87504917383194 
    -0.62495082616806 -0.375049203634262 -0.124950811266899 
    0.124950811266899 0.375049203634262 0.62495082616806 0.87504917383194 
    1.12495076656342 1.37504923343658 1.62495076656342 1.87504923343658 
    2.12495088577271 2.37504911422729 2.62495088577271 2.87504911422729 
    3.12495088577271 3.37504911422729 3.62495088577271 3.87504911422729 
    4.12495088577271 4.37504911422729 4.62495088577271 4.87504911422729 
    5.12495088577271 5.37504911422729 5.62495088577271 5.87504911422729 
    6.12495088577271 6.37504911422729 6.62495088577271 6.87504911422729 
    7.12495088577271 7.37504911422729 7.62495088577271 7.87504911422729 
    8.12495040893555 8.37504959106445 8.62495040893555 8.87504959106445 
    9.12495040893555 9.37504959106445 9.62495040893555 9.87504959106445 
    10.1249504089355 10.375147819519 10.6258354187012 10.8776025772095 
    11.1306390762329 11.3855266571045 11.64244556427 11.9019641876221 
    12.164249420166 12.4298524856567 12.6989212036133 12.9719848632812 
    13.2491683959961 13.5309762954712 13.8175039291382 14.109227180481 
    14.406210899353 14.7088975906372 15.0173177719116 15.3318786621094 
    15.6525726318359 15.9797677993774 16.3134174346924 16.6538486480713 
    17.0009727478027 17.3550758361816 17.7160263061523 18.0840625762939 
    18.4590148925781 18.8410739898682 19.2300262451172 19.6260223388672 
    20.0287990570068 20.4384670257568 20.8547191619873 21.2776222229004 
    21.7068290710449 22.1423664093018 22.5838489532471 23.0312595367432 
    23.4841766357422 23.9425525665283 24.4059276580811 24.8742179870605 
    25.3469352722168 25.8239707946777 26.3048038482666 26.7892990112305 
    27.2769145965576 27.7674942016602 28.26047706604 28.755687713623 
    29.2525539398193 29.7508850097656 30.2500991821289 30.75 31.25 
    31.75 32.25 32.75 33.25 33.75 34.25 34.75 35.25 35.75 36.25 
    36.75 37.25 37.75 38.25 38.75 39.25 39.75 40.25 40.75 41.25 
    41.75 42.25 42.75 43.25 43.75 44.25 44.75 45.25 45.75 46.25 
    46.75 47.25 47.75 48.25 48.75 49.25 49.75 50.25 50.75 51.25 
    51.75 52.25 52.75 53.25 53.75 54.25 54.75 55.25 55.75 56.25 
    56.75 57.25 57.75 58.25 58.75 59.25 59.75 60.25 60.75 61.25 
    61.75 62.25 62.75 63.25 63.75 64.25 64.75 65.25 65.75 66.25 
    66.75 67.25 67.75 68.25 68.75 69.25 69.75 70.25 70.75 71.25 
    71.75 72.25 72.75 73.25 73.75 74.25 74.75 75.25 75.75 76.25 
    76.75 77.25 77.75 78.25 78.75 79.25 79.75 80.25 80.75 81.25 
    81.75 82.25 82.75 83.25 83.75 84.25 84.75 85.25 85.75 86.25 
    86.75 87.25 87.75 88.25 88.75 89.25 89.75 
tdef 17520 linear 00z01jan2002 1hr
zdef 1 levels 0
vars 4
mld 0 11,1, 0  **Mixed Layer depth
ild 0 11,1, 0  **Isothermal Layer depth
sss 0 11,1, 0  **Sea Surface Salinity
sst 0 11,1, 0  **Sea Surface Temperature
ENDVARS
fEOF
