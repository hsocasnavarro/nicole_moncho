; Take a model atmosphere and/or a spectral profile in IDL format and
; write a binary file in NICOLE's native format

pro idl_to_nicole,model=m_in,i=i,q=q,u=u,v=v,file=filename,extractx=indx,$
                  extracty=indy
  model=0
  prof=0
  if (n_elements(m_in) gt 0) then model=1
  if (n_elements(i) gt 0) then prof=1
  filemod=filename
  fileprof=filename
  if model eq 1 and prof eq 1 then begin
     print,'Adding suffix .model and .prof to filename in order to distinguish'
     print,' model atmospheres from spectral profiles'
     filemod=filemod+'.model'
     fileprof=fileprof+'.prof'
  endif
  if model eq 1 then begin
     m=m_in
     openw,iunit,filemod,/get_lun,/swap_if_big
     if (size(m.t))(0) eq 1 then begin
        nz=long((size(m.t))(1))
        m=new_model(1,1,nz)
        nx=1l
        ny=1l
     endif else begin
        nx=long((size(m.t))(1))
        ny=long((size(m.t))(2))
        nz=long((size(m.t))(3))
     endelse
     if n_elements(indx) eq 0 then indx=indgen(nx)
     if n_elements(indy) eq 0 then indy=indgen(ny)
     npix=long(nx)*long(ny)
     nnx=long(n_elements(indx))
     nny=long(n_elements(indy))
     tmp=lon64arr(22*nz+13+92)
     ; 18.04
     ;tmp(0)=4049129056382445934
     ;tmp(1)=2314885530823504944
     ; spic18.12
     tmp(0)=8313474953548884334
     tmp(1)=3616722794636863856
     tmp(2)=nny*long64(2)^32 + nnx
     tmp(3)=nz
     writeu,iunit,tmp
     tmp=dblarr(22*nz+13+92)
     for iix=0l, nnx-1 do for iiy=0l, nny-1 do begin
        ix=indx(iix)
        iy=indy(iiy)

        tmp(0:nz-1)=m.z(ix,iy,*)
        tmp(nz:2*nz-1)=m.tau(ix,iy,*)
        tmp(2*nz:3*nz-1)=m.t(ix,iy,*)
        tmp(3*nz:4*nz-1)=m.gas_p(ix,iy,*)
        tmp(4*nz:5*nz-1)=m.rho(ix,iy,*)
        tmp(5*nz:6*nz-1)=m.el_p(ix,iy,*)
        tmp(6*nz:7*nz-1)=m.v_los(ix,iy,*)
        tmp(7*nz:8*nz-1)=m.v_mic(ix,iy,*)
        tmp(8*nz:9*nz-1)=m.b_los_z(ix,iy,*)
        tmp(9*nz:10*nz-1)=m.b_los_x(ix,iy,*)
        tmp(10*nz:11*nz-1)=m.b_los_y(ix,iy,*)
        tmp(11*nz:12*nz-1)=m.b_x(ix,iy,*)
        tmp(12*nz:13*nz-1)=m.b_y(ix,iy,*)
        tmp(13*nz:14*nz-1)=m.b_z(ix,iy,*)
        tmp(14*nz:15*nz-1)=m.v_x(ix,iy,*)
        tmp(15*nz:16*nz-1)=m.v_y(ix,iy,*)
        tmp(16*nz:17*nz-1)=m.v_z(ix,iy,*)
        tmp(17*nz:18*nz-1)=m.nH(ix,iy,*)
        tmp(18*nz:19*nz-1)=m.nHminus(ix,iy,*)
        tmp(19*nz:20*nz-1)=m.nHplus(ix,iy,*)
        tmp(20*nz:21*nz-1)=m.nH2(ix,iy,*)
        tmp(21*nz:22*nz-1)=m.nH2plus(ix,iy,*)
        tmp(22*nz)=m.v_mac(ix,iy)
        tmp(22*nz+1)=m.stray_frac(ix,iy)
        tmp(22*nz+2)=m.ffactor(ix,iy)
        tmp(22*nz+3)=m.keep_el_p(ix,iy)
        tmp(22*nz+4)=m.keep_gas_p(ix,iy)
        tmp(22*nz+5)=m.keep_rho(ix,iy)
        tmp(22*nz+6)=m.keep_nH(ix,iy)
        tmp(22*nz+7)=m.keep_nHminus(ix,iy)
        tmp(22*nz+8)=m.keep_nHplus(ix,iy)
        tmp(22*nz+9)=m.keep_nH2(ix,iy)
        tmp(22*nz+10)=m.keep_nH2plus(ix,iy)
        tmp(22*nz+11)=m.chrom_x(ix,iy)
        tmp(22*nz+12)=m.chrom_y(ix,iy)
        tmp(22*nz+12+1:22*nz+12+92)=m.abundance(ix,iy,0:91)
        i=22*nz+92+13
        tmp[i]=m.spic_nz[ix,iy] & i=i+1
        tmp[i]=m.spic_nlambda[ix,iy] & i=i+1
        tmp[i:i+maxspic_nz]=m.spic_z[ix,iy,*] & i=i+maxspic_nz
        tmp[i:i+maxspic_nz*maxspic_nlambda]=m.spic_boundary_int[ix,iy,*] & i=i+maxspic_nz*maxspic_nlambda
        tmp[i]=m.spic_temp[ix,iy] & i=i+1
        tmp[i]=m.spic_dens_factor[ix,iy] & i=i+1
        tmp[i]=m.spic_ds[ix,iy] & i=i+1
        tmp[i]=m.spic_doppler[ix,iy] & i=i+1
        
        writeu,iunit,tmp
     endfor
     free_lun,iunit
  endif     
  if prof eq 1 then begin
     if n_elements(q) eq 0 then q=i*0.
     if n_elements(u) eq 0 then u=i*0.
     if n_elements(v) eq 0 then v=i*0.
     openw,iunit,fileprof,/get_lun,/swap_if_big
     nx=(size(i))(1)
     ny=(size(i))(2)
     if n_elements(indx) eq 0 then indx=indgen(nx)
     if n_elements(indy) eq 0 then indy=indgen(ny)
     npix=long(nx)*long(ny)
     nnx=n_elements(indx)
     nny=n_elements(indy)
     nlam=(size(i))(3)
     npix=long(nx)*long(ny)
     tmp=lon64arr(nlam*4)
     tmp(0)=3328834590979877230
     tmp(1)=2314885530823713331
     tmp(2)=nny*long64(2)^32 + nnx
     tmp(3)=nlam
     writeu,iunit,tmp
     tmp=dblarr(nlam*4)
     for iix=0l, nnx-1 do for iiy=0l, nny-1 do begin
        ix=indx(iix)
        iy=indy(iiy)

        tmp(indgen(nlam)*4)=i(ix,iy,*)
        tmp(indgen(nlam)*4+1)=q(ix,iy,*)
        tmp(indgen(nlam)*4+2)=u(ix,iy,*)
        tmp(indgen(nlam)*4+3)=v(ix,iy,*)
        writeu,iunit,tmp
     endfor
     free_lun,iunit
  endif     

end
