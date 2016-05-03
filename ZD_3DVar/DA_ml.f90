module DA_ml
implicit none
logical, parameter ::     &
  DEBUG_DA=.false.,       &   ! general purpose debug messages
  DEBUG_DA_OBS=.false.,   &   ! observation info
  DEBUG_DA_3DV=.false.,   &   ! 3DVar module
  DEBUG_DA_1STEP=.false., &   ! run only 1 DA step (no adv/chem)
  DEBUG_DA_OUTPUT=.true.     ! hourly output before/after DA step
character(len=*), parameter   ::                &
  DA_NAMELIST="namelist.nml",                   &
  DA_FMT_DEF ="('3DVar@PPP YYYY-MM-DD hh: ',A,'.')",&
  NMC_FMT_DEF="('B-NMC@PPP YYYY-MM-DD hh: ',A,'.')"
character(len=len(DA_FMT_DEF)):: da_fmt_msg=''
character(len=128)            :: da_msg=''
real, save                    :: datim_before,datim_after
endmodule DA_ml
