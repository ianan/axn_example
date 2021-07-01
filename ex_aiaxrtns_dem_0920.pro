pro ex_aiaxrtns_dem_0920

  ; Combine AIA, XRT and NuSTAR to work out the DEM
  ;
  ; Note that this is a quick run through of the code so my region selection of the AIA/XRT maps
  ; will probably not be optimal - using this for your own event/data do a more
  ; careful region selection of the same size area in AIA and XRT
  ; - this caveat goes some way to explaining why the ssw and py versions of this
  ;   are slightly different and why I also try a dem calc using the py dn values
  ;
  ; See ex_aiaxrt_dem_0920.pro for a more details AIA+XRT only one.
  ;
  ; 01-Jul-2021 IGH
  ;
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Load in the XRT data - created in make_xrt_for_python.pro
  fits2map,'XRT_20200912_204028_Be_thin_Open_512.fits',xmap

  ; need to make region smaller than used in python as sub_map here rounds up pixels??
  ; Or something to do with coordinate system sswidl vs sunpy maps?
  xr=[-890,-881]
  yr=[214,220]

  sub_map,xmap,sxmap,xrange=xr,yrange=yr
  ; Python version should be 3.60277
  xdnpxs=mean(sxmap.data)/xmap.dur/16.

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

;  print,xrt_tresp[ids].name
  units_xrt=xrt_tresp[ids[0]].temp_resp_units
  gd=where(xrt_tresp[ids[0]].temp gt 0.0)
  logt_xrt=alog10(xrt_tresp[ids[0]].temp[gd])
  tr_xrt=xrt_tresp[ids].temp_resp[gd]

  ; XRT TR is 5.5 - 8.0 in 0.1, AIA TR is 4.0 - 9.0 in 0.05 so need some tweaking before combining
  gd_aia=where(logt_aia ge 5.0 and logt_aia le 8.0)
  logt=logt_aia[gd_aia]
  nt=n_elements(logt)

  tr=dblarr(nt,8)
  tr[*,0:5]=tr_aia[gd_aia,*]
  tr[*,6]=10^interpol(alog10(tr_xrt),logt_xrt,logt)
  ;  plot,[5,8],[1,1],/nodata,/ylog,yrange=[1e-30,1e-23],xrange=[5,8]
  ;  for i=0,6 do oplot,logt,tr[*,i]

  ; Make the NuSTAR response and get the data all in one from the pha, arf, rmf
  ; uses https://github.com/ianan/nustar_sac/blob/master/idl/make_nstresp.pro

  msdir='/Users/iain/Downloads/sep2020__smallflares_spec/80610208001/event_cl/spec_2039_2045/'
  fname='nu80610208001B06_chu12_N_sr'
  make_nstresp,msdir+fname,ns,ebands=[2.5,4.0]

  print,'NuSTAR ',ns.eid[0],' [count/s]: ',ns.rate[0],' +/- ',ns.erate[0]

  ; combine all data together
  dnin=[adnpxs,xdnpxs,ns.rate[0]]
  ednin=0.2*dnin
  ; as the NuSTAR counting error bigger include that with 20% systematic
  ednin[7]=sqrt(ednin[7]^2+ns.erate[0]^2)

  ; Need to multiply NS response by the AIA area to get DEM in correct units
  area=(xra[1]-xra[0])*(yra[1]-yra[0])*7.25e7*7.25e7
  ; Interpolation the NuSTAR response into the T binning of the other channels
  tr[*,7]=10^interpol(alog10(ns.tresp*area),ns.logt,logt)
  ;  plot,[5,8],[1,1],/nodata,/ylog,yrange=[1e-30,1e-23],xrange=[5,8]
  ;  for i=0,7 do oplot,logt,tr[*,i]


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
  dn2dem_pos_nb, dnin, ednin,tr,logt,temps,dem,edem,elogt,chisq,dn_reg,/gloci

  ; Do some plots
  !p.thick=2
  !p.multi=[0,1,2]
  !p.charsize=1.5
  yr=max(dem)*[1e-4,5]
  clr=225
  window,1,xsize=500,ysize=800
  plot,minmax(logtemps),[1,1],/ylog,xtit='Log!D10!N T', ytit='DEM [cm!U-5!N K!U-1!N]',$
    yrange=yr,ystyle=17,xstyle=17,/nodata,tit='AIA+XRT+NuSTAR: ssw roi (red)'
  loadct,39,/silent
  for i=0,n_elements(mlogt)-1 do oplot,logtemps[i:i+1],dem[i]*[1,1],color=clr
  demax=(dem+edem) < yr[1]
  demin=(dem-edem) > yr[0]
  for i=0,n_elements(mlogt)-1 do oplot, mlogt[i]*[1,1],[demin[i],demax[i]],color=clr


  plot_oo,[1e-1,1e4],[1e-1,1e4],lines=1,xtit='DN_in',ytit='DN_reg',tit='AIA+XRT+NuSTAR: ssw roi (red)'
  oplot,dnin,dn_reg,psym=1,color=250
  chans=[channels_aia,filters,ns.eid[0]]
  for i=0, n_elements(chans)-1 do xyouts, dnin[i]*1.5,dn_reg[i],chans[i],/data,color=clr,chars=1,orien=0

  print,'chisq: ',chisq
  print,'dn_in: ',dnin
  print,'dn_in/dn_reg :', dnin/dn_reg


  stop
end