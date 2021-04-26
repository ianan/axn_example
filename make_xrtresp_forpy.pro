pro make_xrtresp_forpy

  ; Get the XRT temperature response for the below dates and save in format python likes
  ; Note that this is the default APEC way of getting the response should use the
  ; newer way that uses CHIANTI - but doesn't make a huge difference to the response
  ;
  ; 26-04-2021 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~

  dates=['28-Sep-2018','12-Sep-2020']
  nd=n_elements(dates)

  for d=0, nd-1 do begin
    wave_resp = make_xrt_wave_resp(contam_time=dates[d])
    xrt_tresp = make_xrt_temp_resp(wave_resp, /apec_default)

    print,wave_resp[1].contam.thick_time,wave_resp[1].contam.thick
    ; For this only need 'Al-poly', 'Be-thin'
    filts=xrt_tresp.name
    for f=0, n_elements(filts)-1 do begin
      idf=strpos(filts[f],';')
      filts[f]=strmid(filts[f],0,idf)
    endfor
    id1=where(filts eq 'Al-poly')
    id2=where(filts eq 'Be-thin')
    ids=[id1,id2]
    filters=filts[ids]

    print,xrt_tresp[ids].name
    units=xrt_tresp[ids[0]].temp_resp_units
    gd=where(xrt_tresp[ids[0]].temp gt 0.0)
    logt=alog10(xrt_tresp[ids[0]].temp[gd])
    tr=xrt_tresp[ids].temp_resp[gd]

    date=anytim(dates[d],/ccsds)

    save,file='xrt_tresp_'+strmid(break_time(dates[d]),0,8)+'.dat',filters,logt,tr,units,date
  endfor
  stop
end