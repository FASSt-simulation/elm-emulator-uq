;================================================
;  topo_5.ncl
;================================================
; Concepts illustrated:
;   - Drawing two topographic maps using a subset of 1' data
;   - Using "MeshFill" for faster contouring
;   - Using functions for cleaner code
;   - Using cnFillPalette to assign a color palette to contours
;   - Using two different colormaps on one frame
;----------------------------------------------------------------------
; This script draws subsets of a global 1' topo grid downloaded from:
;
; http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/cell_registered/netcdf/
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;----------------------------------------------------------------------
; This procedure draws subsets of a global 1' topographic map read in
; from a netCDF file. The first plot is of land, the second of ocean.
; Two different colormaps and contour levels are used.
;----------------------------------------------------------------------
begin
;----------------------------------------------------------------------
  slat=(/38.8953,42.71,35.978,38.895,44.452/)
  slon=(/-105.5464,-80.36,-79.09,-120.6,-121.6/)
  labels=(/"US-NR1","CA-TP4","US-Dk3","US-Blo","US-Me2"/)
 
 wks = gsn_open_wks("pdf","fig2_5ENF_sites")   ; send graphics to PNG file
;
; Increase memory for contours. This is necessary if you are
; contouring a large grid. Otherwise, you might get this error:
;
;     fatal:ContourPlotDraw: Workspace reallocation would exceed maximum size 262144
;
  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize" : 100000000000000000
  end setvalues
;---Read data
  topo_file = "ETOPO1_Bed_c_gmt4.grd.nc"
  a         = addfile(topo_file,"r")
  elev_lnd  = a->z({25:60},{-125:70})
  elev_lnd  = where(elev_lnd.lt.-90,elev_lnd@_FillValue,elev_lnd)
printVarSummary(elev_lnd)
;---Set some resources
  res                   = True

  res@gsnDraw           = False    ; will panel later
  res@gsnFrame          = False    ; will panel later

  res@cnFillOn          = True     ; turn on contour fill
  res@lbOrientation     = "Vertical"

  res@cnLinesOn         = False    ; turn off contour lines
  res@cnLineLabelsOn    = False    ; turn off line labels
  res@cnFillMode        = "MeshFill"
  res@lbBoxLinesOn      = False

  res@gsnAddCyclic      = False

;---Zoom in on map
  res@mpMinLatF         = 25
  res@mpMaxLatF         = 60
  res@mpMinLonF         = -125
  res@mpMaxLonF         = -70
  res@mpCenterLonF      = (res@mpMinLonF + res@mpMaxLonF) / 2.

  res@mpDataBaseVersion = "MediumRes"   ; Improve map outlines
  res@mpFillDrawOrder   = "Draw"

  res@gsnLeftString     = ""            ; make sure these
  res@gsnRightString    = ""            ; are blank

;---Get color map for land
  cmap = read_colormap_file("OceanLakeLandSnow")
  res@cnFillPalette        = cmap(3:,:)

;---Set contours for land levels
  mnmxint = nice_mnmxintvl( min(elev_lnd), max(elev_lnd), 18, False)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = mnmxint(0)
  res@cnMaxLevelValF       = mnmxint(1)
  res@cnLevelSpacingF      = mnmxint(2)/8.
  res@lbTitleString              = "m"
  res@lbTitlePosition            = "Right"            ; title position
  res@lbTitleDirection           = "Across"
  res@lbTitleFontHeightF         = 0.025
 res@lbLabelJust        = "CenterLeft"  ; left justify labels
  res@mpOceanFillColor     = "blue"
  res@mpLandFillColor      = -1            ; transparent

  plot_lnd  = gsn_csm_contour_map(wks,elev_lnd,res)

  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFont     = 22
  txres@txFontColor ="black"
  do i = 0, 4
      text= gsn_add_text(wks,plot_lnd,labels(i),slon(i)+1,slat(i)-1.5,txres)
  end do
 
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.


    resP                = True
    ;resP@gsnMaximize    = True
    resP@gsnFrame       = False
    resP@gsnPanelLeft   = 0.04
    resP@gsnPanelRight  = 0.9
    resP@gsnPanelTop    = 0.95
    resP@gsnPanelYWhiteSpacePercent = 5
    resP@gsnPanelXWhiteSpacePercent = 0
    resP@gsnPanelSave   = True
  gsn_panel(wks,plot_lnd,(/1,1/),resP)

 do i = 0, 4
    if (.not.ismissing(slat(i)))
      gsres@gsMarkerColor      ="black"
      gsres@gsMarkerThicknessF = 10
    gsn_polymarker(wks,plot_lnd, slon(i),slat(i),gsres)
    end if
  end do

frame(wks)
end
;----------------------------------------------------------------------
