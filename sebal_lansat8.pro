pro SEBAL_Lansat8

  ; This code is developed for the "Water and electricity authority of Khoozestan province, Iran"
  ; ---------------------------------------------------------------------------------------------
  ; The goal is to estimate evapotranspiration through SEBAL algorithem.
  ; --------------------------------------------------------------------
  
  ; ---------------------------------------------------------------------------------------
  ; Essential inputs: ---------------------------------------------------------------------
  ; 01. Lansat8 Metadata ------------------------------------------------------------------
  ; 02. Spatial Subset Coordinates (ENVI image coordinates) -------------------------------
  ; 03. Average Elevation -----------------------------------------------------------------
  ; 04. Sun Elevation (in metadata file) --------------------------------------------------
  ; 05. Earth sun distance (in metadata file) ---------------------------------------------
  ; 06. Instant temperature (Kelvin) ------------------------------------------------------
  ; 07. Instant temperature (Celsius)------------------------------------------------------
  ; 08. Wind Speed (10m) (From synoptic data) ---------------------------------------------
  ; 09. Maximum temperature (From synoptic data) ------------------------------------------
  ; 10. Minimum temperature (From synoptic data) ------------------------------------------
  ; 11. Julian day (image date) -----------------------------------------------------------
  ; 12. Synoptic station latitude ---------------------------------------------------------
  ; 13. Synoptic station longitude --------------------------------------------------------
  ; 14. Albedo (From synoptic data) -------------------------------------------------------
  ; 15. Humidity (From synoptic data) -----------------------------------------------------
  ; 16. Decimal time (From synoptic data) -------------------------------------------------
  ; 17. Sunshine hours --------------------------------------------------------------------  
  ; 18. Hot & Cold pixel Data (Net radiation - Soil heat flux - Land surface temperature) -
  ; 19. Minimum elevation -----------------------------------------------------------------
  ; ---------------------------------------------------------------------------------------

  COMPILE_OPT idl2


; Adding ENVI Commands to Program
  e = ENVI(/headless)
  IF e EQ !NULL THEN e = ENVI(/headless)
; -------------------------------------------------------------------

; Data Import Section
  Landsat_Image = DIALOG_PICKFILE(TITLE='Import Landsat MTL. File')
; -------------------------------------------------------------------

; Landsat MTL. File to Raster
  Landsat = e.OpenRaster(Landsat_Image)
; -------------------------------------------------------------------
 
; Defining spectral Subsets
  Multi_Spectral = Landsat[0]
  Thermal_Infrared = Landsat[3]
; -------------------------------------------------------------------

; Coordinate System
  Cor = Multi_Spectral.SpatialRef
; -------------------------------------------------------------------

; Defining Spatial Subset Area
  ULX = 0.0
  ULY = 0.0
  LRX = 0.0
  LRY = 0.0

  READ, ULX, PROMPT= 'Upper Left X = '
  READ, ULY, PROMPT= 'Upper Left Y = '
  READ, LRX, PROMPT= 'Lower Right X = '
  READ, LRY, PROMPT= 'Lower Right Y = '

  Sub_Frame = [ULX,ULY,LRX,LRY]
; -------------------------------------------------------------------

; Multispectral Data Spatial Subset
  MS_Sub = ENVITask('SubsetRaster')
  MS_Sub.Input_Raster = Multi_Spectral
  MS_Sub.Sub_Rect = Sub_Frame
  MS_Sub.Execute
; -------------------------------------------------------------------

; Spatial Subset Thermal Data
  TIR_Sub = ENVITask('SubsetRaster')
  TIR_Sub.Input_Raster = Thermal_Infrared
  TIR_Sub.Sub_Rect = Sub_Frame
  TIR_Sub.Execute
; -------------------------------------------------------------------

; Calibration Section

 ; Multispectral Data Radiance Calibration
  MS_Rad = ENVITask('RadiometricCalibration')
  MS_Rad.Input_Raster = MS_Sub.Output_Raster
  MS_Rad.Calibration_Type = 'Radiance'
  MS_Rad.Execute

 ; Multispectral Data TOA Reflectance Calibration
  MS_Ref_TOA = ENVITask('RadiometricCalibration')
  MS_Ref_TOA.Input_Raster = MS_Sub.Output_Raster
  MS_Ref_TOA.Calibration_Type = 'Top-Of-Atmosphere Reflectance'
  MS_Ref_TOA.Execute

 ; Thermal Data Radiance Calibration
  TIR_Rad = ENVITask('RadiometricCalibration')
  TIR_Rad.Input_Raster = TIR_Sub.Output_Raster
  TIR_Rad.Calibration_Type = 'Radiance'
  TIR_Rad.Execute

 ; Thermal Data Brightness Tempreture Calibration
  TIR_BT = ENVITask('RadiometricCalibration')
  TIR_BT.Input_Raster = TIR_Sub.Output_Raster
  TIR_BT.Calibration_Type = 'Brightness Temperature'
  TIR_BT.Execute

 ; Multispectral Data Atmospheric Correction
  MS_Ref_Surf = ENVITask('QUAC')
  MS_Ref_Surf.Input_Raster = MS_Rad.Output_Raster
  MS_Ref_Surf.Sensor = 'Landsat TM/ETM/OLI'
  MS_Ref_Surf.Execute
; -------------------------------------------------------------------

; Spectral Indices : NDVI and LAI
  MS_VEG_Ind = ENVITask('SpectralIndices')
  MS_VEG_Ind.Input_Raster = MS_Ref_Surf.Output_Raster
  MS_VEG_Ind.Index = ['Normalized Difference Vegetation Index','Leaf Area Index']
  MS_VEG_Ind.Execute

  NDVI_LAI = MS_VEG_Ind.Output_Raster
  NDVI = NDVI_LAI.GetData(Bands = [0])
  LAI = NDVI_LAI.GetData(Bands = [1])
; -------------------------------------------------------------------

; Emissivity
  NDVImin = MIN(NDVI)
  NDVImax = MAX(NDVI)
  PV = (NDVI-NDVImin)/(NDVImax-NDVImin)
  Emissivity = 0.004 * PV + 0.986
  Emissivity[WHERE(NDVI LT 0.2,/NULL)] = 0.97
  Emissivity[WHERE(NDVI GT 0.5,/NULL)] = 0.99
; -------------------------------------------------------------------

; Land Surface Temperature (LST)
  TIR_BT = TIR_BT.Output_Raster
  TIR_BT = TIR_BT.GetData(Bands = [0])
  TIR_Rad = TIR_Rad.Output_Raster
  TIR_Rad = TIR_Rad.GetData(Bands = [0])
  LST_Kel = (TIR_BT) / (1 + (TIR_Rad * (TIR_BT / 14380.0) * ALOG( Emissivity)))
  LST_Cel =  LST_Kel - 273.15
  
  PRINT, ''
  PRINT, '**************************************'
  PRINT, 'Define Land Surface Temperature Output'
  PRINT, '**************************************'
  
  LST_CelEx = ENVIRaster(LST_Cel, SPATIALREF=Cor, URI=DIALOG_PICKFILE(TITLE='Land Surface Temperature Output'))
  LST_CelEx.Save
; -------------------------------------------------------------------

; Albedo
 ;ESUN Values for each Landsat band:
  ESUN_1 = 1719.0
  ESUN_2 = 1787.0
  ESUN_3 = 1746.0
  ESUN_4 = 1536.0
  ESUN_5 = 997.0
  ESUN_6 = 811.0
  ESUN_7 = 75.0
  ESUN_SUM = ESUN_1 + ESUN_2 + ESUN_3 + ESUN_4 + ESUN_5 + ESUN_6 + ESUN_7

 ;Omega rate for each band:
  W_2 = ESUN_2 / ESUN_SUM
  W_3 = ESUN_3 / ESUN_SUM
  W_4 = ESUN_4 / ESUN_SUM
  W_5 = ESUN_5 / ESUN_SUM
  W_6 = ESUN_6 / ESUN_SUM
  W_7 = ESUN_7 / ESUN_SUM

 ;TOA Reflectance for each band:
  MS_Ref_TOA = MS_Ref_TOA.Output_Raster
  Ref_TOA_2 = MS_Ref_TOA.GetData(Bands = [1])
  Ref_TOA_3 = MS_Ref_TOA.GetData(Bands = [2])
  Ref_TOA_4 = MS_Ref_TOA.GetData(Bands = [3])
  Ref_TOA_5 = MS_Ref_TOA.GetData(Bands = [4])
  Ref_TOA_6 = MS_Ref_TOA.GetData(Bands = [5])
  Ref_TOA_7 = MS_Ref_TOA.GetData(Bands = [6])

  Albedo_TOA = (Ref_TOA_2 * W_2) + (Ref_TOA_3 * W_3) + (Ref_TOA_4 * W_4) + (Ref_TOA_5 * W_5) + (Ref_TOA_6 * W_6) + (Ref_TOA_7 * W_7)

  Elevation = 0.0
  READ, Elevation, PROMPT= 'Average Elevation : '

  T_sw = 0.75 + (2 * 0.00001) * Elevation
  Albedo_Surf = (Albedo_TOA - 0.03) / (T_sw ^ 2.0)
; -------------------------------------------------------------------
 
; Incoming Shortwave Radiation (RS_Income)
  Sun_Elevation = 0.0
  READ, Sun_Elevation, PROMPT= 'Sun Elevation : '
  
  Teta = 90 - Sun_Elevation
  Teta_d = (Teta * !PI)/180.0
  
  Earth_Sun_Distance = 0.0
  READ, Earth_Sun_Distance , PROMPT= 'Earth Sun Distance : '
  
  dr = 1 / (Earth_Sun_Distance ^ 2.0)
  G_sc = 1367.0
  
  RS_In = G_sc * COS(Teta_d) * dr * T_sw
  
  PRINT, ''
  PRINT, 'RS_In = ', RS_In
; -------------------------------------------------------------------

; Outgoing Longwave Radiation (RL_Outgoing)
  LAI = NDVI_LAI.Getdata(Bands = [1])
  Emissivity_0 = 0.95 + (0.01 * LAI)
  Emissivity_0[WHERE(LAI GE 3.0,/NULL)] = 0.98
  StBlts_Const = ((5.67) * (10 ^ (-8.0)))
  RL_Out = Emissivity_0 * StBlts_Const * (LST_Kel ^ 4.0)
; -------------------------------------------------------------------

; Incoming Longwave Radiation
  Emissivity_A = 0.85 * ((-ALOG(T_sw)) ^ 0.09)
  T_inst = 0.0
  READ, T_inst, PROMPT= 'Instant Temperature : '
  RL_In = Emissivity_A * stblts_const * (T_inst ^ 4.0)
  PRINT, ''
  PRINT, 'RL_In = ', RL_In
; -------------------------------------------------------------------

; Net Radiation : Rn
  Rn = (1 - Albedo_Surf) * RS_In + RL_In - RL_Out - (1 - Emissivity_0) * RL_In
  
  PRINT, ''
  PRINT, '***************************'
  PRINT, 'Define Net Radiation Output'
  PRINT, '***************************'
  
  RnEx = ENVIRaster(Rn, SPATIALREF=Cor, URI= DIALOG_PICKFILE(TITLE='Net Radiation Output'))
  RnEx.Save
; -------------------------------------------------------------------

; Soil Heat Flux : G
  G_Rn = (LST_Cel / Albedo_Surf) * (((0.0038 * Albedo_Surf) + (0.0074 * (Albedo_Surf ^ 2.0)) * (1-(0.98 * (NDVI ^ 4.0)))))
  G = Rn * G_Rn
  
  PRINT, ''
  PRINT, '****************************'
  PRINT, 'Define Soil Heat Flux Output'
  PRINT, '****************************'
  
  GEx = ENVIRaster(G, SPATIALREF=Cor, URI= DIALOG_PICKFILE(TITLE='Soil Heat Flux Output'))
  GEx.Save
; -------------------------------------------------------------------

; Sensible Heat Flux : H : rah
  K = 0.41
  Zx = 10.0
  Ux = 0.0
  READ, Ux, PROMPT= 'Wind_Speed: '
  
  Z_om = 0.018 * LAI
  U_Star_10m = (K * Ux) / (ALOG(Zx/Z_om))
  U_200 = U_Star_10m * ((ALOG(200.0/Z_om) / K))
  U_Star_200m = ((K * U_200) / (ALOG(200.0/Z_om)))
  
  Z1 = 0.1
  Z2 = 2.0
  
  Rah = (ALOG(Z2 / Z1)) / (U_Star_200m * K)
; -------------------------------------------------------------------

; Hot and Cold Pixel Mapping
  Hot_Cold_Map = (NDVI GT 0.85 AND LST_Kel LT 309) * 1 + (NDVI GE 0 AND NDVI LT 0.2 AND LST_Kel GT 320) * 2
  
  PRINT, ''
  PRINT, '****************************'
  PRINT, 'Define Hot & Cold Map Output'
  PRINT, '****************************'
  
  Hot_Cold_MapEx = ENVIRaster(Hot_Cold_Map, SPATIALREF=Cor, URI= DIALOG_PICKFILE(TITLE='Hot & Cold Map Output'))
  Hot_Cold_MapEx.Save
; -------------------------------------------------------------------

; ET Refernce : Penman Montith
  T_ins_Cel = 0.0
  T_min = 0.0
  T_max = 0.0
  Julian = 0.0
  Lat = 0.0
  Lon = 0.0
  Albedo = 0.0
  Humidity = 0.0
  Time = 0.0
  Sunshine = 0.0

  READ, T_ins_Cel , PROMPT= 'Instant temperature (Celsius) : '
  READ, T_min , PROMPT= 'Minimum temperature : '
  READ, T_max , PROMPT= 'Maximum temperature : '
  READ, Julian , PROMPT= 'Julian Date : '
  READ, Lat , PROMPT= 'Station Latitude : '
  READ, Lon , PROMPT= 'Station Longitude : '
  READ, Albedo , PROMPT= 'Station Albedo : '
  READ, Humidity , PROMPT= 'Rate of Humidity : '
  READ, Time , PROMPT= 'Time (Decimal) : '
  READ, Sunshine , PROMPT= 'Sunshine Hours : '

  phi = (lat*!PI)/180.0
  delt_1 = (2503.0*exp((17.27*T_ins_Cel)/(T_ins_Cel+237.3)))/(T_ins_Cel+237.3)^2.0
  delt_2 = 0.409*sin(((2*!PI*julian)/365.0)-1.39)
  omega_s = ACOS(-TAN(phi)*TAN(delt_2))
  n = (24.0*omega_s)/(!PI)
  b = (2*!PI*(julian-81.0))/364.0
  Sc = 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
  Lz = 51.43
  w = (!PI/12.0)*((time+(lz-lon)/15.0)+Sc-12.0)
  w1 = w-(!PI/24.0)
  w2 = w+(!PI/24.0)
  ra = ((12.0*4.92*0.96)/!PI)*((w2-w1)*sin(phi)*sin(delt_2)+cos(phi)*cos(delt_2)*(sin(w2)-sin(w1)))
  rs = (0.25+0.5*(sunshine/n))*ra
  rn_s = (1-albedo)*rs
  e0_t = 0.6108*exp((17.27*T_ins_Cel)/(T_ins_Cel+237.3))
  e_a = (humidity/100.0)*e0_t
  rs_o = (0.75+2*(0.00001)*elevation)*ra
  fcd = 1.35*(rs/rs_o)-0.35
  rn_l = (2.042*0.0000000001)*fcd*(0.34-0.14*(e_a^0.5))*(T_ins_Cel+273.15)^4.0
  rn_ref = (rn_s)-(rn_l)
  g_ref = 0.04*rn_ref
  p = 101.3*((293-0.0065*elevation)/293.0)^5.26
  gamma = 0.000665*p
  u2 = (Ux*4.87)/(alog(67.8*10.0-5.42))
  e0_t_min = 0.6108*exp((17.27*t_min)/(t_min+237.3))
  e0_t_max = 0.6108*exp((17.27*t_max)/(t_max+237.3))
  e_s = (e0_t_max+e0_t_min)/2.0
  cn = 66.0
  cd = 0.25
  et_reference = ((0.408*delt_1*(rn_ref-g_ref)+(gamma*cn)*u2*(e_s-e_a))/(T_ins_Cel+273.15))/(delt_1+gamma*(1+cd*u2))
  PRINT, ''
  PRINT, 'ET_Refrence = ', et_reference
; -------------------------------------------------------------------

; Latent Heat Flux : Landa
  Landa = (2.501 - (0.00236 * LST_Cel)) * (10 ^ 6.0)
; -------------------------------------------------------------------

; Soil Heat Flux For Hot & Cold Pixels 
  Gr = 9.81
  cp = 1004.0

  Rn_Cold = 0.0
  G_Cold = 0.0
  Rn_Hot = 0.0
  G_Hot = 0.0
  LST_Cold = 0.0
  LST_Hot = 0.0
  MinElevation = 0.0
  
  READ, Rn_Cold , PROMPT= 'Cold Pixel Net radiation : '
  READ, G_Cold , PROMPT= 'Cold Pixel Soil Heat Flux : '
  READ, LST_Cold, PROMPT= 'Cold Pixel Land Surface Temp. : '
  READ, Rn_Hot , PROMPT= 'Hot Pixel Net radiation : '
  READ, G_Hot , PROMPT= 'Hot pixel Soil Heat Flux : '  
  READ, LST_Hot, PROMPT= 'Hot Pixel Land Surface Temp. : '
  READ, MinElevation, PROMPT= 'Minimum Elevation : '
  
  Temp_Rah = Rah
  Temp_U_Star = U_Star_200m
  
 ; Iteration Loop (10 times)
  FOR i = 1, 10 do begin
    H_Cold = (Rn_Cold - G_Cold) - 1.05 * 1000.0 * ((et_reference * 0.001) / (3600.0)) * Landa
    H_Hot = Rn_Hot - G_Hot

   ; Near Surface Temperature
    Air_Pressure = 101.3 * (((293 - 0.0065 * Elevation)/293.0) ^ 5.26)
    Rho =  (1000.0 * Air_Pressure) / (287.0 * 1.01 * LST_Kel)
    dT_Cold = (H_Cold * Temp_Rah) / (Rho * cp)
    dT_Hot = (H_Hot * Temp_Rah) / (Rho * cp)

    Dif_dT = dT_Hot - dT_Cold
    Dif_LST = LST_Hot - LST_Cold

    a = Dif_dT / Dif_LST
    b = (-a) * LST_Hot + dT_Hot

    LST_DEM = LST_Kel + 0.0065 * (Elevation - MinElevation)

    dT_Total = a * LST_DEM + b

   ; H : Initial Sensible Heat Flux
    H = (rho * cp) * (dT_Total / Temp_Rah)

   ; L : Monine Abokhov
    L = (-((Rho * cp * (Temp_U_Star ^ 3.0) * LST_Kel) / (K * Gr * H )))

   ; Sai Parameters
    Sai_200m = FLTARR(SIZE(L,/DIMENSIONS))
    Sai_2m = Sai_200m
    Sai_01m = Sai_200m

    Unstability = WHERE(L LT 0,/NULL, count)

    IF count NE 0 THEN BEGIN

      X_200m = (1.0 - 16.0 * (200.0 / L[Unstability])) ^ 0.25
      Sai_200m[Unstability] = 2 * ALOG((1 + X_200m) / (2.0)) + ALOG( (1 + (X_200m ^ 2.0)) / 2.0) - 2 * ATAN(X_200m) + 0.5 * !PI

      X_2m = (1.0 - 16.0 * (2.0 / L[Unstability])) ^ 0.25
      Sai_2m[Unstability] = 2 * ALOG( (1 + (X_2m ^ 2.0)) / 2.0)

      X_01m = (1.0 - 16.0 * (0.1 / L[Unstability])) ^ 0.25
      Sai_01m[Unstability] = 2 * ALOG( (1 + (X_01m ^ 2.0)) / 2.0)

    ENDIF

    Stability = WHERE(L GT 0,/NULL, count)

    IF count NE 0 THEN BEGIN

      Sai_200m[Stability] = (-5) * (2.0 / L[Stability])
      Sai_2m[Stability] = (-5) * (2.0 / L[Stability])
      Sai_01m[Stability] = (-5) * (0.1 / L[Stability])

    ENDIF

   ; Corrected Rah
    U_Star_200m_Cor = (U_200 * K) / (ALOG( 200.0 / Z_om) - Sai_200m)
    Rah_Cor = (ALOG(Z2/Z1) - Sai_2m + Sai_01m) / (U_Star_200m_Cor * K)

    Temp_Rah = Rah_Cor
    Temp_U_Star = U_Star_200m_Cor

   ; testing the loop
    print, 'Iteration num. ', i

  ENDFOR
; -------------------------------------------------------------------

; Final Instant Evapotranspiration
  LET = Rn - G - H
  ET_ins = (3600.0 * LET) / ((2.501 - 0.00236 * (LST_Kel - 273.15)) * (10 ^ 6.0))

  PRINT, ''
  PRINT, '****************************************'
  PRINT, 'Define Instant Evapotranspiration Output'
  PRINT, '****************************************'

  ET_insEx = ENVIRaster(ET_ins, SPATIALREF=Cor, URI=DIALOG_PICKFILE(TITLE='Instant Evapotranspiration Output'))
  ET_insEx.SAVE

  PRINT, ''
  PRINT, '*******                 *******'
  PRINT, '*******************************'
  PRINT, '**** End of code execution ****'
  PRINT, '*******************************'
  PRINT, '*******                 *******'
  PRINT, ''

end