!Reads the local fractions type emissions, and rewrite them in "new" format
!uses about 30GB memory for a global file (0.5 deg resolution)
!mpif90 -r8 -O2 -o RWEmis Country_mod.f90 ReadWriteEmis.f90  `nc-config --fflags` `nc-config --flibs`
!debug:
!mpif90 -r8  -o RWEmis -debug-parameters all -traceback -ftrapuv -g -fpe0 -O0 -fp-stack-check Country_mod.f90 ReadWriteEmis.f90  `nc-config --fflags` `nc-config --flibs`
program rwemis
  use netcdf
  use Country_mod
  implicit none
  character(len = 200) ::path,fileName,filenameOut,varname,speciesnames(100),projection,periodicity,SECTORS_NAME
  integer :: i,j, imax, jmax, ic, ix
  integer :: itime, ntime, isec, nsectors, ispec, nspecies, found, size, iland,nlandfound
  real  :: lonstart,latstart,dlon,dlat
  integer :: NCMAX = 10 !Max number of countries in one gridcell
  integer :: countrycodes(MAXNLAND),Countryicode(MAXNLAND)
  real, allocatable :: cdfemis(:,:),cdfemisOut(:,:),emisSecC(:,:,:,:),fractions(:,:,:)
  real, allocatable :: landcode(:,:,:),nlandcode(:,:),speciessum(:,:)
  integer :: ncFileID, ncFileIDOut
  logical :: createVariableOnly
  
  call init_Country()
  filename='/nobackup/forsk/sm_petwi/Data/Emis_ECLIPSEv6b_GNFR/Emis_ECLIPSEv6b_CLE_GNFR_05deg_2015.nc'
  filenameOut='testNew.nc'
  path='notset'
  read(*,'(A)')path
  read(*,'(A)')filename
  filenameOut=trim(filename)
  filename=trim(path)//'/'//trim(filename)
  write(*,*)'convert '//trim(filename)
  write(*,*)' into '//trim(filenameOut)

  if(trim(path)=='./' .or. trim(path)=='.' .or. trim(path)=='')then
    write(*,*)'Do not want to overwrite!'
    stop
  end if

  !1) read metadata
  nsectors = 13 !modify if different
  call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))
  call ReadMetadata(ncFileID, imax, jmax, ntime, nspecies, lonstart, dlon, latstart, dlat,projection, speciesnames, SECTORS_NAME)
  if(SECTORS_NAME/='GNFR')then
    write(*,*)'Only GNFR sectors fully implemented'
    stop
  end if
  if (ntime==1) then
     periodicity='yearly'
  else if(ntime==12) then
     write(*,*)'only yearly files implemented'
     stop
     periodicity='monthly'
  else
     write(*,*)'only yearly and monthly files implemented'
     stop
     periodicity='time'
  end if
  
  write(*,*)'grid sizes',imax, jmax, ntime
  write(*,*)'projection  and number of time stamps ',trim(projection),ntime
  write(*,*)'periodicity ',trim(periodicity)
  write(*,*)'number of species and sectors',nspecies,nsectors
  write(*,*)'longitudes start and step: ',lonstart, dlon
  write(*,*)'latitudes start and step: ',latstart, dlat
  
  do i=1,MAXNLAND
     Countryicode(i)=Country(i)%icode
  end do
  allocate(cdfemis(imax,jmax))
  allocate(speciessum(imax,jmax))
  allocate(cdfemisOut(imax,jmax))
  allocate(emisSecC(imax,jmax,nsectors,MAXNLAND))!NB: very large array
  allocate(fractions(imax,jmax,NCMAX))
  allocate(landcode(imax,jmax,NCMAX))
  allocate(nlandcode(imax,jmax))
  nlandcode = 0
  landcode = 0
  size=imax*jmax
  !2) create emisfile with metadata
  call createNetCDF(filenameOut,imax,jmax,nsectors,ntime,projection,periodicity,SECTORS_NAME,lonstart,latstart,dlon,dlat)
  call check(nf90_open(trim(fileNameOut),nf90_share+nf90_write,ncFileIDOut)  )
  !first collect a list of all countries defined
  !we assume that the countries lsit does not change with time 
  do itime = 1, 1
     call readCDF(ncFileID, 'Codes', landcode, size*NCMAX, itime, found)
     call readCDF(ncFileID, 'NCodes', nlandcode, size, itime, found)
  end do
  countrycodes = 0
  nlandfound = 0
  do i=1,imax
     do j=1,jmax
        do ic=1,nint(nlandcode(i,j))
           call find_index(nint(landcode(i,j,ic)),Countryicode, MAXNLAND, ix)
           if (ix<1) then
              write(*,*)'Country code not found ',nint(landcode(i,j,ic)),i,j
              stop
           end if
           if(ix==IC_DUMMY)then
              write(*,*)'WARNING: cannot handle undefined country yet (N/A gives invalid character for NetCDF)'
              stop
           end if
           if (countrycodes(ix)==0) nlandfound = nlandfound+1
           countrycodes(ix)=countrycodes(ix)+1
        end do
     end do
  end do
  
  write(*,*)'found emissions defined for ', nlandfound,' countries'
  !first we only create all the variables. This makes the file much smaller.
  createVariableOnly = .true.
  do itime = 1, ntime
     do ispec = 1, nspecies
!        if(trim(speciesnames(ispec))/='nox' .and.trim(speciesnames(ispec))/='sox' .and.trim(speciesnames(ispec))/='co' .and.trim(speciesnames(ispec))/='pm25' .and.trim(speciesnames(ispec))/='pmco' .and.trim(speciesnames(ispec))/='nh3' .and.trim(speciesnames(ispec))/='voc' )cycle
        do ic=1,MAXNLAND
           if (countrycodes(ic)<=0) cycle !we only write defined countries              
           write(varname,"(A)")trim(country(ic)%code)
           call writeCDFsector(ncFileIDOut, trim(speciesnames(ispec)), trim(varname), emisSecC(1,1,1,ic), imax,jmax,nsectors, itime, createVariableOnly)
        end do
        !4b) write totals for all countries and sectors
        varname=''
        call writeCDFsector(ncFileIDOut, trim(speciesnames(ispec)), trim(varname), speciessum, imax,jmax,0, itime, createVariableOnly)
     end do
  end do
  write(*,*)'all variables created'
  call check(nf90_close(ncFileIDOut))
  call check(nf90_open(trim(fileNameOut),nf90_share+nf90_write,ncFileIDOut)  )
  
  createVariableOnly = .false.
  do itime = 1, ntime
     do ispec = 1, nspecies
!        if(trim(speciesnames(ispec))/='nox' .and.trim(speciesnames(ispec))/='sox' .and.trim(speciesnames(ispec))/='co' .and.trim(speciesnames(ispec))/='pm25' .and.trim(speciesnames(ispec))/='pmco' .and.trim(speciesnames(ispec))/='nh3' .and.trim(speciesnames(ispec))/='voc' )cycle
        write(*,*)' Reading ',trim(speciesnames(ispec))

        emisSecC=0.0
        speciessum=0.0
        do isec = 1, nsectors
           !3) read data and country fractions
           write(varname,"(A,I2.2)")trim(speciesnames(ispec))//'_sec',isec
!           write(*,*)'Reading ',trim(varname)
           call readCDF(ncFileID, 'fractions_'//varname, fractions, size*NCMAX, itime, found)
           cdfemis=0.0
           call readCDF(ncFileID, varname, cdfemis, size, itime, found)
           !convert into country and sector array
           do i=1,imax
              do j=1,jmax
                 do ic=1,nint(nlandcode(i,j))
                    call find_index(nint(landcode(i,j,ic)),Countryicode, MAXNLAND, ix)
                    emisSecC(i,j,isec,ix)=emisSecC(i,j,isec,ix) + cdfemis(i,j)*fractions(i,j,ic)
                    speciessum(i,j)=speciessum(i,j) + cdfemis(i,j)*fractions(i,j,ic)
                 end do
              end do
           end do
        end do
        !4) write data country by country
        !note that we write zero in fields that are not defined           
        do ic=1,MAXNLAND
           if (countrycodes(ic)<=0) cycle !we only write defined countries              
           write(varname,"(A)")trim(country(ic)%code)
           call writeCDFsector(ncFileIDOut, trim(speciesnames(ispec)), trim(varname), emisSecC(1,1,1,ic), imax,jmax,nsectors, itime,createVariableOnly)
        end do
        !4b) write totals for all countries and sectors
        varname=''
        call writeCDFsector(ncFileIDOut, trim(speciesnames(ispec)), trim(varname), speciessum, imax,jmax,0, itime,createVariableOnly)
     end do
  end do
  call check(nf90_close(ncFileID))
  call check(nf90_close(ncFileIDOut))
end program rwemis

subroutine Check(status,errmsg)
  use netcdf
  implicit none
  integer, intent ( in) :: status
  character(len=*), intent(in), optional :: errmsg

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    if(present(errmsg)) print *, "ERRMSG: ", trim(errmsg)
    write(*,*)'NetCDF error?'
    stop
  end if
end subroutine Check

subroutine createNetCDF(filename,imax,jmax,nsectors,ntime,projection,periodicity,SECTORS_NAME,lonstart,latstart,dlon,dlat)
  use netcdf
  implicit none
  character(len=*),  intent(in)  :: fileName,projection,periodicity,SECTORS_NAME
  integer ,  intent(in)  :: imax,jmax,nsectors,ntime
  real ,  intent(in) :: lonstart,latstart,dlon,dlat
  integer :: i,j,timeDimID, iDimID, jDimID, varID, secDimID,ivarID, jvarID ,secvarID,timevarID,ncFileID
  character(len=10) :: created_date,created_hour
  real , allocatable :: xcoord(:),ycoord(:)
  integer, parameter :: NSECTORS_GNFR=13
  integer :: sectornb(NSECTORS_GNFR)
  character (len=1) ::sector1(NSECTORS_GNFR)
  character (len=40) ::sectorname(NSECTORS_GNFR),varname
  call check(nf90_create(fileName,nf90_hdf5,ncFileID),"create:"//trim(fileName))
  call check(nf90_def_dim(ncFileID,"lon",imax,iDimID),"dim:lon")
  call check(nf90_def_dim(ncFileID,"lat",jmax,jDimID),"dim:lat")
  call check(nf90_def_dim(ncFileID,"sector",nsectors,secDimID),"dim:sec")
  !NB: if the time dimension is defined as unlimited, the reading is VERY slow if many variables.
  !seems to increase as the square of the number of variables in the file, even if the variables are not read(!!!)
!  call check(nf90_def_dim(ncFileID,"time",nf90_unlimited,timeDimID),"dim:time")
  call check(nf90_def_dim(ncFileID,"time",ntime,timeDimID),"dim:time")

  ! define coordinate variables
  call check(nf90_def_var(ncFileID,"lon",nf90_float,[idimID],ivarID))
  call check(nf90_put_att(ncFileID,ivarID,"standard_name","longitude"))
  call check(nf90_put_att(ncFileID,ivarID,"long_name","longitude"))
  call check(nf90_put_att(ncFileID,ivarID,"units","degrees_east"))
  call check(nf90_def_var(ncFileID,"lat",nf90_float,[jdimID],jvarID))
  call check(nf90_put_att(ncFileID,jvarID,"long_name","latitude"))
  call check(nf90_put_att(ncFileID,jvarID,"units","degrees_north"))
  call check(nf90_put_att(ncFileID,jvarID,"standard_name","latitude"))
  call check(nf90_def_var(ncFileID,"sector",nf90_int,[secdimID],secvarID))
  call check(nf90_put_att(ncFileID,secvarID,"long_name",trim(SECTORS_NAME)//" sector index"))
  call check(nf90_def_var(ncFileID,"time",nf90_float,[timedimID],timevarID))
 
  ! Write global attributes
  call check(nf90_put_att(ncFileID,nf90_global,"Conventions", "CF-1.6 for coordinates" ))
  call Date_And_Time(date=created_date,time=created_hour)
  call check(nf90_put_att(ncFileID,nf90_global,"created_date", created_date))
  call check(nf90_put_att(ncFileID,nf90_global,"created_hour", created_hour))
  call check(nf90_put_att(ncFileID,nf90_global,"projection",trim(projection)))
  call check(nf90_put_att(ncFileID,nf90_global,"periodicity",trim(periodicity)))
  call check(nf90_put_att(ncFileID,nf90_global,"SECTORS_NAME",trim(SECTORS_NAME)))
  if(trim(SECTORS_NAME)=='GNFR')then
       sector1(1)='A'
       sectorname(1)='PublicPower'
       sector1(2)='B'
       sectorname(2)='Industry'
       sector1(3)='C'
       sectorname(3)='OtherStationaryComb'
       sector1(4)='D'
       sectorname(4)='Fugitive'
       sector1(5)='E'
       sectorname(5)='Solvents'
       sector1(6)='F'
       sectorname(6)='RoadTransport'
       sector1(7)='G'
       sectorname(7)='Shipping'
       sector1(8)='H'
       sectorname(8)='Aviation'
       sector1(9)='I'
       sectorname(9)='Offroad'
       sector1(10)='J'
       sectorname(10)='Waste'
       sector1(11)='K'
       sectorname(11)='AgriLivestock'
       sector1(12)='L'
       sectorname(12)='AgriOther'
       sector1(13)='M'
       sectorname(13)='Other'
       do i = 1,NSECTORS_GNFR
          write(varname,"(A,I2.2)")'sec',i
          call check(nf90_put_att(ncFileID,nf90_global,trim(varname),sector1(i)//'_'//trim(sectorname(i))))
       end do
    end if
  call check(nf90_enddef(ncFileID))

  allocate(xcoord(imax))
  do i=1,imax
     xcoord(i)=(i-1)*dlon + lonstart
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord) )

  allocate(ycoord(jmax))
  do j=1,jmax
     ycoord(j)=(j-1)*dlat + latstart
  enddo
  call check(nf90_put_var(ncFileID, jVarID, ycoord) )

  do i=1,NSECTORS_GNFR
     sectornb(i)=i
  enddo
  call check(nf90_put_var(ncFileID, secVarID, sectornb) )
    
  call check(nf90_close(ncFileID))

  deallocate(xcoord)
  deallocate(ycoord)

end subroutine createNetCDF


subroutine find_index(wanted, list, lsize, Index)
  integer, intent(in) :: wanted
  integer, dimension(1:lsize), intent(in) :: list
  integer ::   n
  Index = -1
  do n = 1, lsize
     if (wanted == list(n)) then
        Index = n
        exit
     end if
  end do
end subroutine find_index
  
subroutine find_indexc(wanted, list, lsize, Index)
  character(len=*), intent(in) :: wanted
  character(len=*), dimension(1:lsize), intent(in) :: list
  integer ::   n
  Index = -1
  do n = 1, lsize
     if (wanted == list(n)) then
        Index = n
        exit
     end if
  end do
end subroutine find_indexc
  
subroutine  createnewsectorvariable(ncFileID,varname,speciesname,countrycode,units)
  use netcdf
  use Country_mod
  implicit none
  integer ,intent(in) ::ncFileID
  character (len = *),intent(in) ::varname,speciesname,countrycode, units
  integer :: idimID, jdimID,secdimID,timeDimID,varID,xtype,ndims
  character*100 ::name
  integer :: imax, jmax, nsectors, ix
  
  call check(nf90_redef(ncid = ncFileID))
  call check(nf90_inq_dimid(ncFileID,"lon"  ,idimID),"dim:lon")
  call check(nf90_inq_dimid(ncFileID,"lat"  ,jdimID),"dim:lat")
  call check(nf90_inq_dimid(ncFileID,"sector"  ,secdimID),"dim:sec")
  call check(nf90_inq_dimid(ncFileID,"time"  ,timedimID),"dim:time")

  call check(nf90_inquire_dimension(ncFileID,idimID,len=imax))
  call check(nf90_inquire_dimension(ncFileID,jdimID,len=jmax))
  call check(nf90_inquire_dimension(ncFileID,secdimID,len=nsectors))
  call check(nf90_def_var(ncFileID,varname,nf90_float,&
             [iDimID,jDimID,secDimID,timeDimID],varID),"defvar:"//trim(varname))
!define variable as to be compressed
  call check(nf90_def_var_deflate(ncFileid,varID,shuffle=0,deflate=1,&
       deflate_level=4),"compress:"//trim(varname))
  call check(nf90_def_var_chunking(ncFileID,varID,NF90_CHUNKED,&
       !small chunks give faster read for emep model
       (/min(IMAX,40),min(JMAX,40),nsectors,1/)),"chunk3D:"//trim(varname))
  
  call check(nf90_put_att(ncFileID, varID, "units", units))
  call check(nf90_put_att(ncFileID, varID, "species", trim(speciesname)))
  if(trim(speciesname)=='nox')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 46.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else  if(trim(speciesname)=='sox')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 64.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else   if(trim(speciesname)=='nh3')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 17.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else   if(trim(speciesname)=='co')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 28.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  end if
  call check(nf90_put_att(ncFileID, varID, "country_ISO", trim(countrycode)))
  call find_indexc( trim(countrycode),Country(1:MAXNLAND)%code, MAXNLAND, ix)
  call check(nf90_put_att(ncFileID, varID, "countrycode", Country(ix)%icode))

  call check(nf90_enddef(ncid = ncFileID))
  
end subroutine createnewsectorvariable
subroutine  createnewvariable(ncFileID,varname,units)
  use netcdf
  implicit none
  integer ,intent(in) ::ncFileID
  character (len = *),intent(in) ::varname, units
  integer :: idimID, jdimID,timeDimID,varID,xtype,ndims
  character*100 ::name
  integer :: imax, jmax
  
  call check(nf90_redef(ncid = ncFileID))
  call check(nf90_inq_dimid(ncFileID,"lon"  ,idimID),"dim:lon")
  call check(nf90_inq_dimid(ncFileID,"lat"  ,jdimID),"dim:lat")
  call check(nf90_inq_dimid(ncFileID,"time"  ,timedimID),"dim:time")

  call check(nf90_inquire_dimension(ncFileID,idimID,len=imax))
  call check(nf90_inquire_dimension(ncFileID,jdimID,len=jmax))
  call check(nf90_def_var(ncFileID,varname,nf90_float,&
             [iDimID,jDimID,timeDimID],varID),"defvar:"//trim(varname))
!define variable as to be compressed
  call check(nf90_def_var_deflate(ncFileid,varID,shuffle=0,deflate=1,&
       deflate_level=4),"compress:"//trim(varname))
  call check(nf90_def_var_chunking(ncFileID,varID,NF90_CHUNKED,&
       (/min(IMAX,100),min(JMAX,100),1/)),"chunk2D:"//trim(varname))
  
  call check(nf90_put_att(ncFileID, varID, "units", units))
  if(trim(varname)=='nox')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 46.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else  if(trim(varname)=='sox')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 64.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else   if(trim(varname)=='nh3')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 17.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  else   if(trim(varname)=='co')then
    call check(nf90_put_att(ncFileID, varID, "molecular_weight", 28.0))
    call check(nf90_put_att(ncFileID, varID, "molecular_weight_units", "g mole-1"))
  end if

  call check(nf90_enddef(ncid = ncFileID))
  
end subroutine createnewvariable

subroutine readCDF(ncFileID, varname, Rvar, size, itime, found)
  use netcdf
  implicit none
  character(len = *),intent(in) :: varname
  integer ,intent(in) :: ncFileID,size,itime
  real, intent(out) :: Rvar(size)
  integer ,intent(out) :: found
  integer :: i,j,xtype, ndims,varID,dimids(NF90_MAX_VAR_DIMS),nAtts,status
  character*100 ::name
  integer :: startvec(NF90_MAX_VAR_DIMS),count(NF90_MAX_VAR_DIMS)
  found=0
  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  
  if(status == nf90_noerr) then     
     call check(nf90_Inquire_Variable(ncFileID,VarID,name,&
          xtype,ndims,dimids,nAtts),"GetDimsId")
     startvec=1
     count=0
     do i=1,ndims
        call check(nf90_inquire_dimension(ncid=ncFileID, dimID=dimids(i), &
             len=count(i)),"GetDims")
     end do
     if (count(ndims)<itime) then
        print *, 'time stamp does not exist: ',trim(varname)
        return
     end if
     !only fetch one time (last dimension) 
     startvec(ndims) = itime
     count(ndims) = 1
     call check(nf90_get_var(ncFileID, VarID, Rvar,start=startvec,count=count))
     found = 1
  else
     print *, 'variable does not exist: ',trim(varname)
     stop
  endif

end subroutine readCDF

subroutine ReadMetadata(ncFileID, imax, jmax, ntime, nspecies, lonstart, dlon, latstart, dlat,projection, speciesnames, SECTORS_NAME)
  use netcdf
  implicit none
  integer ,intent(in) :: ncFileID
  character(len=200), intent(out) :: speciesnames(100),projection, SECTORS_NAME
  integer ,intent(out) :: imax, jmax, ntime, nspecies
  real, intent(out) :: lonstart, dlon, latstart, dlat
  integer :: i,j,n,xtype, varID,dimids(NF90_MAX_VAR_DIMS),nAtts,status
  character*100 ::name
  integer :: startvec(NF90_MAX_VAR_DIMS),count(NF90_MAX_VAR_DIMS)
  integer :: ndims, nVariables, nAttributes, unlimitedDimId
  real :: rvar(2)
  
  status = nf90_get_att(ncFileID, nf90_global, "SECTORS_NAME", SECTORS_NAME )
  if(status /= nf90_noerr) then     
    write(*,*)'Did not find global attribute SECTORS_NAME. Assuming GNFR'
    SECTORS_NAME='GNFR'
  end if
  call check(nf90_get_att(ncFileID, nf90_global, "projection", projection))
  if(trim(projection) /= 'lon lat')then
     write(*,*)'Sorry, can only handle lon lat projection!'
     stop
  end if
  status = nf90_inq_varid(ncid = ncFileID, name = 'Codes', varID = VarID)
  if(status == nf90_noerr) then     
     call check(nf90_Inquire_Variable(ncFileID,VarID,name,&
          xtype,ndims,dimids,nAtts),"GetDimsId")
     count=0
     do i=1,ndims
        call check(nf90_inquire_dimension(ncFileID,dimids(i),len=count(i)),"len:dimN")
     end do
     imax = count(1)
     jmax = count(2)
     ntime = count(ndims)
  else
     write(*,*)'expected a file with Codes as variable '
     stop
  end if
  call check(nf90_inq_varid(ncid = ncFileID, name = 'lon', varID = VarID))
  call check(nf90_get_var(ncFileID, VarID, Rvar,count=[2]))
  lonstart = nint(Rvar(1)*100000)/100000.
  dlon = nint((Rvar(2)-Rvar(1))*100000)/100000.
  call check(nf90_inq_varid(ncid = ncFileID, name = 'lat', varID = VarID))
  call check(nf90_get_var(ncFileID, VarID, Rvar,count=[2]))
  latstart = nint(Rvar(1)*100000)/100000.
  dlat = nint((Rvar(2)-Rvar(1))*100000)/100000.
  
  call check(nf90_inquire(ncFileID, ndims, nVariables, nAttributes))
  nspecies = 0
  do i=1, nVariables
     call check(nf90_Inquire_Variable(ncFileID,i,name,&
          xtype,ndims,dimids,nAtts),"InqVar")
     !find name that start with "fractions_" and end with "_sec01"
     n=len_trim(name)
     if(n>5) then
        if(name(n-5:n)=='_sec01' .and. name(1:10)=='fractions_')then
           nspecies = nspecies + 1
           speciesnames(nspecies) = name(11:n-6)
           write(*,*)nspecies,' found ',trim(speciesnames(nspecies))
        end if
     end if
  end do

end subroutine ReadMetadata
subroutine writeCDFsector(ncFileID, speciesname, countrycode, Rvar, imax, jmax,nsectors, itime, createVariableOnly)
  use netcdf
  implicit none
  character(len = *),intent(in) :: speciesname, countrycode
  integer ,intent(in) :: ncFileID, imax, jmax, itime, nsectors
  real, intent(in) :: Rvar(imax,jmax,nsectors)
  logical, intent(in) :: createVariableOnly !NB: bus error if made optional !?
  integer :: i,j,xtype, ndims,varID,dimids(NF90_MAX_VAR_DIMS),nAtts,status
  character*600 ::name,varname,units
  integer :: startvec(NF90_MAX_VAR_DIMS),count(NF90_MAX_VAR_DIMS)
  
  if(nsectors>0)then
     write(varname,"(A)")trim(speciesname)//'_'//trim(countrycode)
  else
     write(varname,"(A)")trim(speciesname)
  end if
  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID)
  if(status /= nf90_noerr) then     
    !make new variable
     units = 'tonnes/year'
     if(nsectors>0)then
        call createnewsectorvariable(ncFileID,varname,speciesname,countrycode,units)            
     else
        call createnewvariable(ncFileID,varname,units)
     end if
  endif
  if(createVariableOnly) return

  call check(nf90_inq_varid(ncid = ncFileID, name = trim(varname), varID = VarID))
  if (nsectors>0) then
     startvec=1
     startvec(4)=itime
     count=0
     count(1)=imax
     count(2)=jmax
     count(3)=nsectors
     count(4)=1
  else
     startvec=1
     startvec(3)=itime
     count=0
     count(1)=imax
     count(2)=jmax
     count(3)=1
  end if
  call check(nf90_put_var(ncFileID, VarID, Rvar,start=startvec,count=count) )
  call check(nf90_inq_varid(ncid = ncFileID, name = 'time', varID = VarID))
  call check(nf90_put_var(ncFileID, VarID, [itime],start=[itime],count=[1]) )!overwrite even if already defined for simplicity

end subroutine writeCDFsector
