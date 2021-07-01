pro ex_aiaxrt_dem_0920

  ; Combine AIA and XRT to work out the DEM
  ;
  ; Note that this is a quick run through of the code so my region selection
  ; will probably not be optimal - using this for your own event/data do a more
  ; careful region selection of the same size area in AIA and XRT
  ; - this caveat goes some way to explaining why the ssw and py versions of this
  ;   are slightly different and why I also try a dem calc using the py dn values
  ;
  ; 30-Jun-2021 IGH
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;

  ; Load in the XRT data - created in make_xrt_for_python.pro
  fits2map,'XRT_20200912_204028_Be_thin_Open_512.fits',xmap

  ; need to make region smaller than used in python as sub_map here rounds up pixels??
  ; Or something to do with coordinate system sswidl vs sunpy maps?
  xr=[-890,-881]
  yr=[214,220]

  sub_map,xmap,sxmap,xrange=xr,yrange=yr
  ; Python version should be 3.60277
  xdnpxs=mean(sxmap.data)/xmap.dur/16.
  print,xdnpxs
  ;loadct,3,/silent
  ;plot_map,sxmap

  ; Load in AIA submaps prepd in the python notebook
  ff=file_search('','maps_prep_092021*.fits')
  fits2map,ff,amap
  ; Tweak aia submap but same area as XRT one? Need to be careful with the big XRT pixels
  xra=[-890,-881]
  yra=[217,223]
  sub_map,amap,samap,xrange=xra,yrange=yra
  ;  plot_map,samap[2]
  adnpxs=dblarr(6)
  for i=0,5 do begin
    adnpxs[i]=mean(samap[i].data)/samap[i].dur
  endfor

  print,adnpxs

  ; Combine the data together
  dnin=[adnpxs,xdnpxs]

  ;  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Make the response functions
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; AIA
  ; Note that in python version have generic t0 responses and
  ; correct the data for the degradation
  ; would need to use result=aia_bp_get_corrections()
  ; for that appraoch here...
  tresp=aia_get_response(/temperature,/dn,/eve,timedepend_date='12-Sep-2020')

  ids=[0,1,2,3,4,6]
  channels_aia=tresp.channels[ids]
  logt_aia=tresp.logte
  tr_aia=tresp.all[*,ids]
  units_aia=tresp.units

  ; XRT
  ; Using the "default" appraoch here, need to check Kathy's updated one
  wave_resp = make_xrt_wave_resp(contam_time='12-Sep-2020')
  xrt_tresp = make_xrt_temp_resp(wave_resp, /apec_default)

  ; Only need 'Be-thin'
  filts=xrt_tresp.name
  for f=0, n_elements(filts)-1 do begin
    idf=strpos(filts[f],';')
    filts[f]=strmid(filts[f],0,idf)
  endfor
  ;  id1=where(filts eq 'Al-poly')
  id2=where(filts eq 'Be-thin')
  ids=id2;[id1,id2]
  filters=filts[ids]

  print,xrt_tresp[ids].name
  units_xrt=xrt_tresp[ids[0]].temp_resp_units
  gd=where(xrt_tresp[ids[0]].temp gt 0.0)
  logt_xrt=alog10(xrt_tresp[ids[0]].temp[gd])
  tr_xrt=xrt_tresp[ids].temp_resp[gd]

  ; XRT TR is 5.5 - 8.0 in 0.1, AIA TR is 4.0 - 9.0 in 0.05 so need some tweaking before combining
  gd_aia=where(logt_aia ge 5.0 and logt_aia le 8.0)
  logt=logt_aia[gd_aia]
  nt=n_elements(logt)

  tr=dblarr(nt,7)
  tr[*,0:5]=tr_aia[gd_aia,*]
  tr[*,6]=10^interpol(alog10(tr_xrt),logt_xrt,logt)
  ;  plot,[5,8],[1,1],/nodata,/ylog,yrange=[1e-30,1e-23],xrange=[5,8]
  ;  for i=0,6 do oplot,logt,tr[*,i]


  ;  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; So have the data and TR ready to run the DEM code
  ; run the regularization

  ; Choose  T binning for the DEM
  mint=5.6
  maxt=6.8
  ntb=42.
  temps=10^(mint+(maxt-mint)*findgen(ntb)/(ntb-1.))
  logtemps=alog10(temps)
  mlogt=get_edges(logtemps,/mean)


  ; from https://github.com/ianan/demreg/blob/master/idl/dn2dem_pos_nb.pro
  dn2dem_pos_nb, dnin, 0.2*dnin,tr,logt,temps,dem,edem,elogt,chisq,dn_reg,/gloci

  ; What about using the numbers from python directly as the locations slighlty different in the submaps
  ; between ssw and sunpy - so when doing with your own data can just visually check your roi actually
  ; contains the thing you want....
  dninpy=[2.16585241,  24.95020394, 457.63229199, 576.8995456,  204.86727998,  13.59922846,   3.60276666]
  dn2dem_pos_nb, dninpy, 0.2*dninpy,tr,logt,temps,dempy,edempy,elogtpy,chisqpy,dn_regpy,/gloci

  ; Do some plots
  !p.thick=2
  !p.multi=[0,1,2]
  !p.charsize=1.5
  yr=max(dem)*[1e-4,5]
  window,1,xsize=400,ysize=800
  plot,minmax(logtemps),[1,1],/ylog,xtit='Log!D10!N T', ytit='DEM [cm!U-5!N K!U-1!N]',$
    yrange=yr,ystyle=17,xstyle=17,/nodata,tit='ssw roi (red), py roi (blue)'
  loadct,39,/silent
  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dem[i]*[1,1],color=250
  demax=(dem+edem) < yr[1]
  demin=(dem-edem) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[demin[i],demax[i]],color=250

  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dempy[i]*[1,1],color=80
  demaxpy=(dempy+edempy) < yr[1]
  deminpy=(dempy-edempy) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[deminpy[i],demaxpy[i]],color=80

  plot_oo,[1e-1,1e4],[1e-1,1e4],lines=1,xtit='DN_in',ytit='DN_reg',tit='ssw roi (red), py roi (blue)'
  oplot,dnin,dn_reg,psym=1,color=250
  oplot,dninpy,dn_regpy,psym=7,color=80,thick=1
  chans=[channels_aia,filters]
  for i=0, 6 do xyouts, dninpy[i],dn_regpy[i],chans[i],/data,color=80,chars=1
  for i=0, 6 do xyouts, dnin[i],dn_reg[i]*2,chans[i],/data,color=250,chars=1,orien=90


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; What about the errors in AIA
  anumpix=n_elements(samap[0].data)
  aerrper=dblarr(6)
  achan=[94,131,171,193,211,335]
  ; do \loud to get a printout of the error components
  for i=0,5 do aerrper[i]=100.*aia_bp_estimate_error(adnpxs[i]*samap[i].dur,achan[i],n_sample=anumpix)/$
    (adnpxs[i]*samap[i].dur)
  print,'AIA % error: ',aerrper
  ; so still smaller than the systematics we are using


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; What about the errors in XRT....

  stop
end