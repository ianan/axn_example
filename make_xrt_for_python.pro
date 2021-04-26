pro make_xrt_for_python

  ; Prep data from Sep 2018 and Sep 2020 and save out in a format sunpy is happy with.
  ; Save out the data as well as the grade maps
  ;
  ; 26-04-2021  IGH

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Example from Sep 2018 data
  ff18=file_search('','XRT20180928_203102.9.fits')

  read_xrt,ff18,ind,data,/force
  ; grade_type=1 as these are summed data of 512x512 (so each pixel is sum of 4x4)
  xrt_prep,ind,data,indp,datap,/float,grade_map=gm,grade_type=1,/coalign
  ; Change so compatiable with sunpy better
  indp.timesys='UTC'
  filtout=indp.ec_fw1_+'_'+indp.ec_fw2_
  resout=strcompress(string(indp.naxis1),/rem)
  fnameout='XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
  write_xrt,indp,datap,outdir='',outfile=fnameout,/ver
  gfnameout='gm_XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
  write_xrt,indp,gm,outdir='',outfile=gfnameout,/ver


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Example from Sep 2020 data
  ff20=file_search('','XRT20200912_204028.8.fits')
  read_xrt,ff20,ind,data,/force
  xrt_prep,ind,data,indp,datap,/float,grade_map=gm,grade_type=1,/coalign
  indp.timesys='UTC'
  filtout=indp.ec_fw1_+'_'+indp.ec_fw2_
  resout=strcompress(string(indp.naxis1),/rem)
  fnameout='XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
  write_xrt,indp,datap,outdir='',outfile=fnameout,/ver
  gfnameout='gm_XRT_'+break_time(indp.date_obs)+'_'+filtout+'_'+resout+'.fits'
  write_xrt,indp,gm,outdir='',outfile=gfnameout,/ver




  stop
end