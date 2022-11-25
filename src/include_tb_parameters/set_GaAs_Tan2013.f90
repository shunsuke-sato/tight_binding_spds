! GaAs parameters from J Comput Electron 12, 56 (2013)


!  energy_unit_tb = 4.7065*ev
  energy_unit_tb = 1d0*ev

  Ea_s   = -4.5863d0*energy_unit_tb
  Ec_s   = -1.3323d0*energy_unit_tb
  Ea_p   = 1.4694d0*energy_unit_tb
  Ec_p   = 9.5885d0*energy_unit_tb
  Ea_d   = 11.2878d0*energy_unit_tb !*0d0+10*ev
  Ec_d   = 35.2863d0*energy_unit_tb !*0d0+10*ev
  Ea_s2  = 10.0480d0*energy_unit_tb !*0d0+10*ev
  Ec_s2  = 25.6752d0*energy_unit_tb !*0d0+10*ev


  ss_sigma      = -1.7615d0*energy_unit_tb
  s2s2_sigma    = -0.8374d0*energy_unit_tb !*0d0
  s2a_sc_sigma  = -1.1173d0*energy_unit_tb !*0d0
  sa_s2c_sigma  = -2.9313d0*energy_unit_tb !*0d0

  sapc_sigma    = 2.1768d0*energy_unit_tb
  scpa_sigma    = 3.6705d0*energy_unit_tb
  s2apc_sigma   = 2.6877d0*energy_unit_tb !*0d0
  s2cpa_sigma   = 1.8335d0*energy_unit_tb !*0d0

  sadc_sigma    = -2.1172d0*energy_unit_tb !*0d0
  scda_sigma    = -2.9128d0*energy_unit_tb !*0d0
  s2adc_sigma   = -0.4974d0*energy_unit_tb !*0d0
  s2cda_sigma   = -2.9971d0*energy_unit_tb !*0d0

  pp_sigma   =  3.8065d0*energy_unit_tb
  pp_pi      = -1.5010d0*energy_unit_tb

  padc_sigma  =  -1.2077d0*energy_unit_tb !*0d0
  pcda_sigma  =  -1.9855d0*energy_unit_tb !*0d0
  padc_pi     =   3.1547d0*energy_unit_tb !*0d0
  pcda_pi     =   2.3234d0*energy_unit_tb !*0d0

  dd_sigma  =  -1.9986d0*energy_unit_tb !*0d0
  dd_pi     =   3.1681d0*energy_unit_tb !*0d0
  dd_delta  =  -2.3137d0*energy_unit_tb !*0d0

  delta_a  = 0.1259d0*energy_unit_tb
  delta_c  = 0.1235d0*energy_unit_tb

! note: the values should be renormalized so as to reproduce the original paper results.
! Most probably, there were some typos in their definision
  delta_a = 3d0*delta_a
  delta_c = 3d0*delta_c
