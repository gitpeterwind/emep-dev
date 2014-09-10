module DA_ml
implicit none
logical, parameter            :: DA_DEBUG=.false.
character(len=*), parameter   ::                &
  DA_NAMELIST="namelist.nml",                   &
  DA_FMT_DEF ="('3DVar YYYY-MM-DD hh: ',A,'.')",&
  NMC_FMT_DEF="('B-NMC YYYY-MM-DD hh: ',A,'.')"
character(len=len(DA_FMT_DEF)):: da_fmt_msg=''
character(len=128)            :: da_msg=''
real, save                    :: datim_before,datim_after
endmodule DA_ml
