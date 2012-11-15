module DA_ml
logical, parameter            :: DA_DEBUG=.false.
character(len=*), parameter   :: DA_FMT_DEF="('3DVar YYYY-MM-DD hh: ',A,'.')",&
                                NMC_FMT_DEF="('B-NMC YYYY-MM-DD hh: ',A,'.')"
character(len=83)             :: da_fmt_msg='', da_msg=''
real, save                    :: datim_before,datim_after
end module DA_ml
