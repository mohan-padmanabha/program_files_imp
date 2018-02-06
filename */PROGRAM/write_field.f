
      subroutine write_field(n1,n2,n3,nfield,time,vx,vy,vz)

      implicit none

      include 'netcdf.inc'

      integer n1,n2,n3,nfield
      real time
      real vx(n1+2,n2+1,n3+1),vy(n1+2,n2+1,n3+1),vz(n1+2,n2+1,n3+1)
     

      integer nid,field(3),status
      integer start(3),count(3),imap(3),stride(3)
      integer iddimx,iddimy,iddimz,idvx,idvy,idvz
      character(len=3) numf


C====================================================================== 
C                         CREATE NetCDF FILE
C====================================================================== 

 11   format(i3.3)

      write(numf,11)nfield

      status = nf_create('FIELD-'/ /numf/ /'.nc',
     &                   nf_write,nid)

      if (status.ne.nf_noerr) call handle_err(status)

      write(*,*)'===> CREATE FILE : ',
     &           'FIELD-'/ /numf/ /'.nc'
      write(*,*)' '

C =============================== Dimension ===========================

      status = nf_def_dim(nid,'resolution_x',n1,iddimx)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_dim(nid,'resolution_y',n2,iddimy)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_dim(nid,'resolution_z',nf_unlimited,iddimz)
      if (status.ne.nf_noerr) call handle_err(status)

C =============================== Variable ===========================

      field(1)=iddimx
      field(2)=iddimy
      field(3)=iddimz

      status = nf_def_var(nid,'velocity_x',nf_float,3,field,idvx)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nid,'velocity_y',nf_float,3,field,idvy)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_def_var(nid,'velocity_z',nf_float,3,field,idvz)
      if (status.ne.nf_noerr) call handle_err(status)

C =============================== Attribut ===========================

      status = nf_put_att_double(nid,nf_global,'Time',
     &                        nf_real,1,time)
      if (status.ne.nf_noerr) call handle_err(status)

      status = nf_enddef(nid)

C====================================================================== 
C                         SAVE Velocity field
C====================================================================== 

      stride(1)=1
      stride(2)=1
      stride(3)=1

      imap(1)=1
      imap(2)=imap(1)*(n1+2)
      imap(3)=imap(2)*(n2+1)

      start(1) = 1
      start(2) = 1
      start(3) = 1

      count(1) = n1
      count(2) = n2
      count(3) = n3

      status=nf_put_varm_double(nid,idvx,start,count,stride,imap,vx)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_put_varm_double(nid,idvy,start,count,stride,imap,vy)
      if (status.ne.nf_noerr) call handle_err(status)

      status=nf_put_varm_double(nid,idvz,start,count,stride,imap,vz)
      if (status.ne.nf_noerr) call handle_err(status)


      status = nf_close(nid)

      return
      end
