      subroutine write_stat(nh,eners,enert,enerf,
     &   reynolds,reynolds_taylor,
     &   integral_scale,taylor_scale,kol_scale,
     &   vrms,energy,inj_energy,dissipation,
     &   turnover_time,time,nstat,
     &   nsid,iddims,idtims,idtime,idreyn,idreyn_tay,idinteg_sc,
     &   idtayl_sc,idkolm_sc,idrms_vel,iddissip,idinj,idturn,idener,
     &   idener_sp,idener_tr,idener_fl)

      implicit none
 
      include 'netcdf.inc'

      integer nh,nstat

      real   eners(0:nh)
      real   enert(0:nh)
      real   enerf(0:nh)

      real reynolds,reynolds_taylor
      real integral_scale,taylor_scale,kol_scale
      real vrms,energy,enstrophy,dissipation
      real turnover_time,time,inj_energy

      real   wavenumbers(0:nh)

      integer status,k
      integer spectr(2)
      integer start_sc,count_sc
      integer start_sp(2),count_sp(2)

      integer nsid,iddims,idtims,idtime,idreyn,idreyn_tay,idinteg_sc
      integer idtayl_sc,idkolm_sc,idrms_vel,iddissip,idinj,idturn,idener
      integer idener_sp,idener_tr,idener_fl,idwn


      if (nstat .eq. 0) then 

C====================================================================== 
C                         CREATE NetCDF FILE
C====================================================================== 

      status = nf_create('STAT.nc',nf_write,nsid)
      if (status.ne.nf_noerr) call handle_err(status)

C ---------------------------- Dimensions ----------------------------

      status = nf_def_dim(nsid,'dim_spectra',nh,iddims)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_dim(nsid,'dim_time',nf_unlimited,idtims)
      if (status.ne.nf_noerr) call handle_err(status)

C --------------------------- Scalar Statistics -----------------------

      status = nf_def_var(nsid,'time',
     &                    nf_float,1,idtims,idtime)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'reynolds',nf_float,1,idtims,idreyn)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'reynolds_taylor',nf_float,1,idtims,
     &                    idreyn_tay)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'integral_scale',nf_float,1,idtims,
     &                    idinteg_sc)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'taylor_scale',nf_float,1,idtims,
     &                    idtayl_sc)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'kolmogorov_scale',nf_float,1,idtims,
     &                    idkolm_sc)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'rms_velocity',nf_float,1,idtims,
     &                   idrms_vel)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'energy',nf_float,1,idtims,
     &                    idener)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'energy_dissipation',nf_float,1,idtims,
     &                   iddissip)
      if (status.ne.nf_noerr) call handle_err(status)

       status = nf_def_var(nsid,'energy_injection',nf_float,1,idtims,
     &                     idinj)
       if (status.ne.nf_noerr) call handle_err(status)

       status = nf_def_var(nsid,'turnover_time',
     &                    nf_float,1,idtims,idturn)
       if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'wavenumbers',nf_float,1,
     &                    iddims,idwn)
      if (status.ne.nf_noerr) call handle_err(status)


C --------------------------- Vector Statistics -----------------------

      spectr(1)=iddims
      spectr(2)=idtims


      status = nf_def_var(nsid,'energy_spectra',nf_float,2,
     &                    spectr,idener_sp)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'energy_transfer',nf_float,2,
     &                    spectr,idener_tr)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nsid,'energy_flux',nf_float,2,
     &                   spectr,idener_fl)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_enddef(nsid)

      do k=0,nh 
        wavenumbers(k) = float(k)
      enddo

      start_sc = 1
      count_sc = nh

      status = nf_put_vara_double(nsid,idwn,start_sc,
     &                          count_sc,wavenumbers)
      if (status.ne.nf_noerr) call handle_err(status)


      endif

C====================================================================== 
C                       WRITE in  NetCDF FILE
C====================================================================== 

      nstat = nstat+1

      start_sc = nstat
      count_sc = 1

      status = nf_put_vara_double(nsid,idtime,start_sc,
     &                          count_sc,time)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idreyn,start_sc,
     &                          count_sc,reynolds)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idreyn_tay,start_sc,
     &                          count_sc,reynolds_taylor)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idinteg_sc,start_sc,
     &                          count_sc,integral_scale)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idtayl_sc,start_sc,
     &                          count_sc,taylor_scale)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idkolm_sc,start_sc,
     &                          count_sc,kol_scale)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idrms_vel,start_sc,
     &                          count_sc,vrms)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idener,start_sc,
     &                          count_sc,energy)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,iddissip,start_sc,
     &                          count_sc,dissipation)
      if (status.ne.nf_noerr) call handle_err(status)

        status = nf_put_vara_double(nsid,idinj,start_sc,
     &                            count_sc,inj_energy)
       if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idturn,start_sc,
     &                          count_sc,turnover_time)
      if (status.ne.nf_noerr) call handle_err(status)

      start_sp(1) = 1
      start_sp(2) = nstat

      count_sp(1) = nh
      count_sp(2) = 1

      status = nf_put_vara_double(nsid,idener_sp,
     &                          start_sp,count_sp,eners)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idener_tr,
     &                          start_sp,count_sp,enert)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_put_vara_double(nsid,idener_fl,
     &                          start_sp,count_sp,enerf)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_sync(nsid)
      if (status.ne.nf_noerr) call handle_err(status)

      return
      end

C====================================================================== 

      subroutine handle_err(status)
      integer status
       print *,'NetCDF error: ',status
       stop
      return
      end
