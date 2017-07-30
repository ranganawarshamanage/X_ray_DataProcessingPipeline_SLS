MODULE SUB_ROUTINES
IMPLICIT NONE

PUBLIC :: createinputfile
PUBLIC :: ranker
PUBLIC :: isa_cutoff
PUBLIC :: cc_set
PUBLIC :: dataset_remover
PUBLIC :: cc_frames
!!$PUBLIC :: frame_remover

CONTAINS
  subroutine createinputfile(i,dat_fmt,xtal,run,tot_rotation_angle,a,b,c,alf,bet,gam,spgr,FDL,ABSOPT,&
       resol_ll,resol_ul,strng_pxl,min_pixel_spot,frac_index,REFDATA,ref_dataset,head_file_name)
    implicit none
    integer, intent(in) :: i,run,spgr,strng_pxl,min_pixel_spot
    real, intent(in) :: a,b,c,alf,bet,gam,resol_ll,resol_ul,frac_index,tot_rotation_angle
    character*10, intent(in) :: dat_fmt,FDL,ABSOPT
    character*20, intent(in),optional :: ref_dataset
    character*10, intent(in) :: xtal
    logical, intent(in) :: REFDATA
    character(len=255) :: line,frame_name,basename,img_path,refdata_path,head_file_name
    real :: exposure_time,wavelength,detector_distance,beamx,beamy,oscillation_angle
    integer :: char_length_frame,char_length_base,char_length_headfilename,nframes,dx,j
    integer :: nframe_strt,nframe_end,nframe_use
    logical :: file_exist

1000 FORMAT(A)
1001 FORMAT(A,A)
1002 FORMAT(A,6F8.2)  
1003 FORMAT(A,F7.2)
1004 FORMAT(A,I5)
1005 FORMAT(A,F4.1,A,F4.1)
1006 FORMAT(A,F7.5)
1007 FORMAT(A,F8.4)
1008 FORMAT(A,F4.2)
1009 FORMAT(A,I5,A,I5)

    file_exist = .false.

    if(dat_fmt=='CBF')then
       open(unit=12,file='frame1.head',status='old')
       do j = 1, 35
          read(12,'(a)') line
          if(j==3)  read(line(1:),*)    frame_name
          if(line(3:15)=='Exposure_time') read(line(17:25),*) exposure_time
          if(line(3:12)=='Wavelength') read(line(14:20),*) wavelength
          if(line(3:19)=='Detector_distance') read(line(21:27),*) detector_distance
          if(line(3:9)=='Beam_xy')then
             read(line(12:18),*) beamx
             read(line(21:27),*) beamy
          end if
          if(line(3:17)=='Angle_increment') read(line(19:24),*) oscillation_angle
       end do
       char_length_frame = LEN_TRIM(frame_name)
       char_length_base  = char_length_frame - 5
       read(frame_name(6:char_length_base),*) basename
       close(12)
       img_path = '../slink/'//trim(xtal)//'/'//trim(basename)//'?????.cbf CBF'
       if(REFDATA) refdata_path = '../slink/wd_path/'//trim(ref_dataset) ! 25.01.2016
       if(tot_rotation_angle /= 0)then
          nframes = int(tot_rotation_angle/oscillation_angle)
       else
          call system("ls ../slink/"//trim(xtal)//'/*.cbf | wc -l > temp.tmp')
          open(11,file='temp.tmp',status='old')
          read(11,*) nframes
          close(11);call system("rm temp.tmp")
       end if
       dx = int(detector_distance*1000)

       open(unit=20,file='XDS.INP',status='new')
       if(run==1) write(20,1000) 'JOB=XYCORR INIT COLSPOT IDXREF'
       if(run==2) write(20,1000) 'JOB=DEFPIX INTEGRATE CORRECT'
       write(20,1000) 'MAXIMUM_NUMBER_OF_JOBS=4'
       write(20,1000) 'MAXIMUM_NUMBER_OF_PROCESSORS=12'
       write(20,1000)
       write(20,1000) '! for this experiment:'
       write(20,1003) 'ORGX=              ',beamx
       write(20,1003) 'ORGY=              ',beamy
       write(20,1004) 'DETECTOR_DISTANCE= ',dx
       write(20,1007) 'OSCILLATION_RANGE= ',oscillation_angle
       write(20,1006) 'X-RAY_WAVELENGTH=  ',wavelength
       write(20,1001) 'NAME_TEMPLATE_OF_DATA_FRAMES=',img_path
       write(20,1004) 'DATA_RANGE=        1 ',nframes 
       write(20,1004) 'SPOT_RANGE=        1 ',nframes
       write(20,1000) 'BACKGROUND_RANGE=  2 10'
       write(20,1000)
       write(20,1004) 'SPACE_GROUP_NUMBER=',spgr
       write(20,1002) 'UNIT_CELL_CONSTANTS= ',a,b,c,alf,bet,gam
       write(20,1000)
       if(REFDATA) write(20,1001) 'REFERENCE_DATA_SET=',trim(adjustl(refdata_path)) !25.01.2016
       write(20,1000) 'REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !DISTANCE'
       write(20,1000) 'REFINE(INTEGRATE)=BEAM ORIENTATION CELL !AXIS DISTANCE'
       write(20,1000) 'REFINE(CORRECT)=BEAM ORIENTATION CELL AXIS !DISTANCE'
       write(20,1000)
       write(20,1001) 'FRIEDEL''S_LAW=',FDL
       write(20,1001) 'STRICT_ABSORPTION_CORRECTION=',ABSOPT
       write(20,1000)
       write(20,1000) '! parameters with changes wrt default values:'
       write(20,1000) 'TRUSTED_REGION=0.00 1.15'
       write(20,1000) 'VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=8000. 30000.'
       write(20,1000) '!MINIMUM_ZETA=0.05'
       write(20,1000) 'CORRECTIONS=DECAY MODULATION ABSORP !default value'
       write(20,1000) '!CORRECTIONS= ! for minisets'
       write(20,1000) '!MINIMUM_I/SIGMA=50 !prevent CORRECT from scaling'
       write(20,1005) 'INCLUDE_RESOLUTION_RANGE=',resol_ll,' ',resol_ul
       write(20,1004) 'STRONG_PIXEL=',strng_pxl
       write(20,1004) 'MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=',min_pixel_spot
       write(20,1008) 'MINIMUM_FRACTION_OF_INDEXED_SPOTS=',frac_index
       write(20,1000) '! parameters specifically for this detector and beamline:'
       write(20,1000) 'DETECTOR=PILATUS NX=2463 NY=2527 QX=0.172  QY=0.172 !PILATUS 6M'
       write(20,1000) 'MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500'
       write(20,1000) 'DIRECTION_OF_DETECTOR_X-AXIS=1 0 0'
       write(20,1000) 'DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0'
       write(20,1000) 'INCIDENT_BEAM_DIRECTION=0 0 1'
       write(20,1000) 'ROTATION_AXIS=1 0 0'
       write(20,1000) 'FRACTION_OF_POLARIZATION=0.99'
       write(20,1000) 'POLARIZATION_PLANE_NORMAL=0 1 0'
       write(20,1000) 'SENSOR_THICKNESS=0.32'
       write(20,1000) '!UNTRUSTED_RECTANGLE= 487  495     1 2527'
       write(20,1000) '!UNTRUSTED_RECTANGLE= 981  989     1 2527'
       write(20,1000) '!UNTRUSTED_RECTANGLE=1475 1483     1 2527'
       write(20,1000) '!UNTRUSTED_RECTANGLE=1969 1977     1 2527'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463   195  213'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463   407  425'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463   619  637'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463   831  849'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  1043 1061'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  1255 1273' 
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  1467 1485'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  1679 1697'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  1891 1909'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  2103 2121'
       write(20,1000) '!UNTRUSTED_RECTANGLE=   1 2463  2315 2333'
       write(20,1000) 'NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13 !used by: INTEGRATE'
       close(20)

    elseif(trim(adjustl(dat_fmt))=='H5')then
       inquire(file='../slink/h5header.tmp',exist=file_exist)
       if(file_exist)then
          open(unit=12,file='../slink/h5header.tmp',status='old')
       else
          print*, 'h5header.tmp doest not exist'
          stop
       end if

       ! write h5header.tmp reading
       do
          read(12,'(a)') line
          if(line(1:19)=='incident_wavelength') read(line(26:33),*) wavelength
          if(line(1:17)=='detector_distance') read(line(26:31),*) detector_distance
          if(line(1:13)=='beam_center_y') read(line(26:32),*) beamy
          if(line(1:13)=='beam_center_x') read(line(26:32),*) beamx
          if(line(1:19)=='omega_range_average') read(line(26:30),*) oscillation_angle
          if(line(1:15)=='frames_per_file')then
             read(line(16:),*) nframes
             exit
          end if
       end do
       close(12)
       ! data path construction
       char_length_headfilename = LEN_TRIM(head_file_name)
       print*, 'head_file_name = ',trim(head_file_name),' characters= ', char_length_headfilename
       char_length_base  = char_length_headfilename - 9
       read(head_file_name(1:char_length_base),*) basename
       img_path = '../slink/'//trim(basename)//'??????.h5'

       ! reference data path
       if(REFDATA) refdata_path = '../slink/wd_path/'//trim(ref_dataset) ! 25.01.2016

       ! first and last frame numbers
       nframe_strt = i*nframes + 1 - nframes      
       nframe_use = int(tot_rotation_angle/oscillation_angle)
       nframe_end = nframe_strt + nframe_use - 1

       ! detector distance
       dx = int(detector_distance*1000)

       open(unit=20,file='XDS.INP',status='new')
       if(run==1) write(20,1000) 'JOB=XYCORR INIT COLSPOT IDXREF'
       if(run==2) write(20,1000) 'JOB=DEFPIX INTEGRATE CORRECT'
       write(20,1000) 'MAXIMUM_NUMBER_OF_JOBS=8'
       write(20,1000) 'MAXIMUM_NUMBER_OF_PROCESSORS=12'
       write(20,1000)
       write(20,1000) '! for this experiment:'
       write(20,1003) 'ORGX=              ',beamx
       write(20,1003) 'ORGY=              ',beamy
       write(20,1004) 'DETECTOR_DISTANCE= ',dx
       write(20,1007) 'OSCILLATION_RANGE= ',oscillation_angle
       write(20,1006) 'X-RAY_WAVELENGTH=  ',wavelength
       write(20,1001) 'NAME_TEMPLATE_OF_DATA_FRAMES=',img_path
       write(20,1009) 'DATA_RANGE=        ',nframe_strt,' ',nframe_end 
       write(20,1009) 'SPOT_RANGE=        ',nframe_strt,' ',nframe_end
       write(20,1009) 'BACKGROUND_RANGE=  ',nframe_strt,' ',nframe_strt+9
       write(20,1000)
       write(20,1004) 'SPACE_GROUP_NUMBER=',spgr
       write(20,1002) 'UNIT_CELL_CONSTANTS= ',a,b,c,alf,bet,gam
       write(20,1000)
       if(REFDATA) write(20,1001) 'REFERENCE_DATA_SET=',trim(adjustl(refdata_path)) !25.01.2016
       write(20,1000) 'REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION'
       write(20,1000) 'REFINE(INTEGRATE)=POSITION BEAM ORIENTATION CELL DISTANCE !AXIS'
       write(20,1000) 'REFINE(CORRECT)=BEAM ORIENTATION CELL AXIS DISTANCE'
       write(20,1000)
       write(20,1001) 'FRIEDEL''S_LAW=',FDL
       write(20,1001) 'STRICT_ABSORPTION_CORRECTION=',ABSOPT
       write(20,1000)
       write(20,1000) '! parameters with changes wrt default values:'
       write(20,1000) 'TRUSTED_REGION=0.00 1.15'
       write(20,1000) 'VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=8000. 30000.'
       write(20,1000) '!MINIMUM_ZETA=0.05'
       write(20,1000) 'CORRECTIONS=DECAY MODULATION ABSORP !default value'
       write(20,1000) '!CORRECTIONS= ! for minisets'
       write(20,1000) '!MINIMUM_I/SIGMA=50 !prevent CORRECT from scaling'
       write(20,1005) 'INCLUDE_RESOLUTION_RANGE=',resol_ll,' ',resol_ul
       write(20,1004) 'STRONG_PIXEL=',strng_pxl
       write(20,1004) 'MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=',min_pixel_spot
       write(20,1008) 'MINIMUM_FRACTION_OF_INDEXED_SPOTS=',frac_index
       write(20,1000) 'CLUSTER_RADIUS=2'
       write(20,1000) 'SEPMIN=4'
       write(20,1000) '! parameters specifically for this detector and beamline:'
       write(20,1000) 'DETECTOR=EIGER NX=4150 NY=4371 QX=0.075  QY=0.075 ! mm'
       write(20,1000) 'MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=3000000'
       write(20,1000) 'SENSOR_THICKNESS=0.32 ! mm'
       write(20,1000)
       write(20,1000) 'DIRECTION_OF_DETECTOR_X-AXIS=1 0 0'
       write(20,1000) 'DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0'
       write(20,1000) 'INCIDENT_BEAM_DIRECTION=0 0 1'
       write(20,1000) 'ROTATION_AXIS=1 0 0'
       write(20,1000) 'FRACTION_OF_POLARIZATION=0.99'
       write(20,1000) 'POLARIZATION_PLANE_NORMAL=0 1 0'
       write(20,1000)
       write(20,1000) 'Detector mask date: 2016.02.29 ??'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151    514  552'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   1065 1103'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   1616 1654'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   2167 2205'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   2718 2756'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   3269 3307'
       write(20,1000) 'UNTRUSTED_RECTANGLE=    0 4151   3820 3858'
       write(20,1000)  
       write(20,1000) 'NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13 !used by: INTEGRATE'
       close(20)
    end if

  end subroutine createinputfile

  SUBROUTINE ranker(FDL,ABSOPT,resol_ll,resol_ul)
    USE SORT
    IMPLICIT NONE
    REAL, INTENT(IN) :: resol_ll,resol_ul
    INTEGER nwell,nn,i,link_len
    REAL Rmeas,rmi
    CHARACTER*10, INTENT(in) :: FDL,ABSOPT
    CHARACTER(LEN=132) :: link_CORRECT,line,command,link_temp,link_xds_ascii
    CHARACTER(LEN=132), DIMENSION(:), ALLOCATABLE :: LINK
    REAL, DIMENSION(:), ALLOCATABLE :: ARR
    INTEGER, DIMENSION(:), ALLOCATABLE :: order

1000 FORMAT(A)
1005 FORMAT(A,F4.1,A,F4.1)
1001 FORMAT(A,A)
    Rmeas = 0

    OPEN(UNIT=10,FILE='linkcorrectfile.txt',STATUS='old')
    CALL SYSTEM("cat linkcorrectfile.txt | wc -l > x.tmp")
    OPEN(UNIT=11,FILE='x.tmp',STATUS='old')
    READ(11,'(I4)') nwell
    CLOSE(11);CALL SYSTEM("rm x.tmp")
    ALLOCATE(LINK(nwell));ALLOCATE(ARR(nwell));ALLOCATE(order(nwell))
    DO nn=1, nwell
       READ(10,1000) link_CORRECT
       link_CORRECT = trim(link_CORRECT)
       CALL SYSTEM("tac " // link_CORRECT // " >  " // "correct.tmp")
       OPEN(UNIT=12,FILE='correct.tmp',STATUS='old')
       OPEN(UNIT=20,FILE='stat.txt',STATUS='replace')
       DO
          READ(12,1000) line
          IF(line(5:34)=='WILSON STATISTICS OF DATA SET') EXIT
       END DO
       DO i=1,10
          READ(12,*)
       END DO
       DO
          READ(12,1000) line
          WRITE(20,1000) line
          IF(line(2:20)=='SUBSET OF INTENSITY')THEN
             EXIT
          END IF
       END DO
       CALL SYSTEM("tac " // "stat.txt" // " > " // "info.txt")
       CALL SYSTEM("echo " // link_CORRECT // " >> " // "outfile.out")
       CALL SYSTEM("rm stat.txt")
       CALL SYSTEM("rm correct.tmp")
       CALL SYSTEM("cat info.txt >> outfile.out")
       OPEN(UNIT=14,FILE='info.txt',STATUS='old')
       DO i=1,4
          READ(14,*)
       END DO
       DO i=1,3
          READ(14,1000) line
          READ(line(93:98),*) rmi
          Rmeas = Rmeas + abs(rmi)
       END DO
       WRITE(21,*) TRIM(ADJUSTL(link_CORRECT)), Rmeas/3
       LINK(nn)=link_CORRECT;ARR(nn)=Rmeas/3
       Rmeas = 0
       CLOSE(14)
    END DO
    CALL quick_sort(ARR,order)
    OPEN(UNIT=22,FILE='XSCALE_GBL_SRT.INP',STATUS='REPLACE')
    WRITE(22,1000) 'OUTPUT_FILE=XSCALE_GBL_SRT.ahkl'
    WRITE(22,1001) 'FRIEDEL''S_LAW=',FDL
    WRITE(22,1000) 'SAVE_CORRECTION_IMAGES=FALSE'
    DO i=1,nwell
       link_len = LEN_TRIM(LINK(order(i)))
       READ(LINK(order(i))(1:link_len-10),1000) link_temp
       link_xds_ascii = TRIM(ADJUSTL(link_temp))//'XDS_ASCII.HKL'
       WRITE(22,1001) 'INPUT_FILE=',trim(adjustl(link_xds_ascii))
       WRITE(22,1005) 'INCLUDE_RESOLUTION_RANGE=',resol_ll,' ',resol_ul
       WRITE(22,1000) 'MINIMUM_I/SIGMA=0'
    END DO
    CLOSE(10);CLOSE(11);CLOSE(12);CLOSE(20);CLOSE(21);CLOSE(22)

  END SUBROUTINE ranker

  subroutine isa_cutoff(FDL,ABSOPT,resol_ll,resol_ul,xscale_isa_template)
    implicit none

    character*10, intent(in) :: FDL,ABSOPT
    real, intent(in)   :: resol_ll,resol_ul
    character(len=255),intent(out) :: xscale_isa_template
    !    character(len=255),intent(out) :: xscale_isa_hkl,xscale_isa_inp
    character(len=120) :: line,dset !,xscale_isa_hkl,xscale_isa_inp
    character*3        :: strng4
    real               :: ISa
    logical            :: isempty
    integer            :: ISa_Th,i,j,nset

    isempty = .false.

1000 FORMAT(A)
1001 FORMAT(A,A)
1004 FORMAT(A,I5)
1005 FORMAT(A,F4.1,A,F4.1)
1009 FORMAT(F5.2)

    nset = 0
    i = 0
    j = 0

    open(unit=13,file='XSCALE.LP',status='old')
    do
       read(13,1000) line
       if(trim(adjustl(line))=='a        b          ISa    ISa0   INPUT DATA SET')then
          do
             i = i + 1  !this i carrys the total number of datasets
             read(13,1000) line
             print*, trim(line)
             write(23,1000) trim(adjustl(line))
             if(len_trim(line)== 0)then
                isempty = .true.
                exit
             end if
          end do
       end if
       if(isempty) exit
    end do
    isempty = .false.
    close(13);close(23)

    write(6,*) 'Enter your ISa Threshold:'
    read(*,*) ISa_Th
    write(strng4,'(I3)') ISa_Th
    xscale_isa_template = 'XSCALE_ISa'//trim(adjustl(strng4))
    !   xscale_isa_hkl = 'XSCALE_ISa'//trim(adjustl(strng4))//'.ahkl'
    !   xscale_isa_inp = 'XSCALE_ISa'//trim(adjustl(strng4))//'.INP'
    open(unit=14,file='fort.23',status='old')
    open(unit=24,file=trim(adjustl(xscale_isa_template))//'.INP',status='replace')
    write(24,1001) 'OUTPUT_FILE=',trim(adjustl(xscale_isa_template))//'.ahkl'
    write(24,1001) 'FRIEDEL''S_LAW=',FDL
    write(24,1000) 'SAVE_CORRECTION_IMAGES=FALSE'
    print*, i
    do j = 1,i-1 
       read(14,1000) line
       print*,line
       read(line(24:28),1009) ISa
       read(line(38:),1000)   dset
       if(ISa > ISa_Th)then
          nset = nset + 1  !this nset carrys the ISa number of datasets
          write(24,1001) 'INPUT_FILE=',trim(dset)
          write(24,1005) 'INCLUDE_RESOLUTION_RANGE=',resol_ll,' ',resol_ul
          write(24,1004) 'MINIMUM_I/SIGMA=',ISa_Th
       end if
    end do
    close(14);close(24)
  end subroutine isa_cutoff

  subroutine cc_set(resol_ll,resol_ul,nset,filename_template,filename_template_ccd)

    USE symops_mod
    USE resolution_mod
    !USE rc_mod   !29072017
    !USE hash_mod !29072017
    USE gnufor2

    IMPLICIT NONE

    ! character*255, intent(in) :: xscale_isa_template
    real, intent(in)   :: resol_ll,resol_ul
    character*255, intent(out) :: filename_template,filename_template_ccd
    integer, intent(out) :: nset
    integer :: unit,ier1
    integer,dimension(:), allocatable :: dataset_i
    real, dimension(:), allocatable :: cc_i

    INTEGER :: ih,ik,il,oldhkl,jh,jk,jl,ihkl,iop,i,j,k,n,ispgr,ncheck, &
         ier2,iset,nuniq,nobs,oldset,oldnuniq, &
         ii,i1,i2,nx,ny,npairs
    REAL s2,resol,fsq,sigma,dum4(4),CC2,SE
    CHARACTER (LEN=132) :: xscale_isa_hkl
    CHARACTER (LEN=80) :: line!,filename_template
    CHARACTER*1 W
    REAL :: sumyy,sumX,sumY,sumXY,sumX2,sumY2,AA,BB,CC,DD,CorC,start,finish
    LOGICAL last

    TYPE :: uniq
       integer uindex
       integer setnum
       real    intensity
       type(uniq), POINTER :: nextuniq
    END TYPE uniq

    TYPE(uniq), POINTER :: ob, firstuniq

    INTEGER, DIMENSION(:), ALLOCATABLE :: uindex, setnum
    REAL,    DIMENSION(:), ALLOCATABLE :: intensity, x, y
    REAL,    DIMENSION(:,:), ALLOCATABLE :: sumI
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: nfreq
    CHARACTER, DIMENSION(:), ALLOCATABLE :: filename_arr

1000 FORMAT(A)

    CALL SYSTEM("ls -l *.ahkl")
    
    unit=10
    WRITE(6,*)
    WRITE(6,*) 'Enter the name of XSCALE.HKL type file?'
    READ(5,'(a)') xscale_isa_hkl

    do i=1,len_trim(xscale_isa_hkl)
       if(xscale_isa_hkl(i:i)=='.') exit
    end do
    read(xscale_isa_hkl(1:i-1),1000) filename_template

!!$    ALLOCATE(filename_arr(len_trim(xscale_isa_hkl)))
!!$    DO i=1,len_trim(xscale_isa_hkl)
!!$       read(xscale_isa_hkl(i:i),'(a)') W
!!$       if(W == '.') exit
!!$       filename_arr(i:i) = W
!!$    END DO
!!$    filename_template = TRANSFER(filename_arr,filename_template)

    print*,trim(adjustl(filename_template))

    OPEN(UNIT=unit,FILE=trim(adjustl(xscale_isa_hkl)),STATUS='old',ACTION='READ',iostat=ier1)
    IF (ier1/=0) THEN
       PRINT*,'ier=',ier1
       STOP 'could not open '!, trim(adjustl(xscale_isa_hkl)) 
    END IF

    filename_template_ccd = trim(adjustl(filename_template))//'_ccd'
    print*,trim(adjustl(filename_template_ccd))
    CALL SYSTEM("cp " // trim(adjustl(filename_template)) // ".INP " &
         // trim(adjustl(filename_template_ccd)) // ".INP" )

    OPEN(UNIT=22,FILE=trim(adjustl(filename_template_ccd))//'.txt',status='replace')

    WRITE(22,*) '% 1st column- dataset number' 
    WRITE(22,*) '% 2nd column- correlation coefficent against all other sets'
    WRITE(22,*) '% 3rd column- standard error in CC'
    WRITE(22,*) '% 4th column- number of common pairs'

    ncheck=0

    ! Open input file; get spacegroup and cell constants
    DO
       READ(unit,1000,iostat=ier2) line
       IF (ier2/=0) THEN
          PRINT*,'ier=',ier2
          PRINT*,'line=',line(:LEN_TRIM(line))
          STOP 'read error in xds_file'
       END IF
       IF (line(1:20)=='!SPACE_GROUP_NUMBER=') then
          READ(line(21:),*) ispgr
          CALL getops(ispgr,op)
          WRITE(*,'(a)'),line(:LEN_TRIM(line))
       ELSE IF (line(1:21).eq.'!UNIT_CELL_CONSTANTS=') then 
          READ(line(22:),*) a,b,c,alf,bet,gam
          WRITE(*,'(a)'),line(:LEN_TRIM(line))
       ELSE IF (line(1:9) == '!ITEM_H=1') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:9) == '!ITEM_K=2') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:9) == '!ITEM_L=3') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:12) == '!ITEM_IOBS=4') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:19) == '!ITEM_SIGMA(IOBS)=5') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:13) == '!ITEM_ISET=10') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:7) == '! ISET=') THEN
          READ(line(8:),*) iset
          nset=max(nset,iset)
       ELSE IF (line(1:14) == '!END_OF_HEADER') THEN
          EXIT
       END IF
    END DO
    IF (ncheck /= 6) stop 'header items not found'
    print*,'nset=',nset

    allocate(dataset_i(nset));allocate(cc_i(nset))
    
    nuniq=0   ! Number of unique reflections
    nobs=0    ! Number of observations

    last=.FALSE.
    ihkl=0
    i = 0

    ALLOCATE(ob); NULLIFY(ob%nextuniq); firstuniq => ob

    DO
       READ(unit,1000) line
       IF (line(1:12) == '!END_OF_DATA') last=.TRUE. 
       IF (.NOT.last) THEN
          READ(line,*)ih,ik,il,fsq,sigma,dum4,iset
          IF (sigma<0.) CYCLE  ! Filter out reflections with sigma <0
          CALL getres(ih,ik,il,s2,resol)
          IF (resol>resol_ll.OR.resol<resol_ul) CYCLE  ! Filter out reflections not within given resolution limit
          nobs=nobs+1  ! Counts the number of observations
          CALL asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop) ! Generate a unique identifier (ihkl) for each reflection.
          ! This identifier is used to identify unique reflections.
!!$        write(20,'(4I10,F12.1,I10)') jh,jk,jl,ihkl,fsq,iset ! For testing
          IF (nobs==1) oldhkl=ihkl
          i = i+1  ! Just a counter and is ranging from 1 to number of observations (nobs)
          IF (i == 1) nuniq = 1  ! First reflection is always unique
          IF (ihkl /= oldhkl.OR.last) THEN  ! process unique reflection
             nuniq = nuniq+1	          ! update # of unique refs
          END IF
          IF (i == 1) THEN
             ob%uindex    = nuniq
             ob%setnum    = iset
             ob%intensity = fsq
!!$           WRITE(21,'(2I10,F12.1)') ob%setnum,ob%uindex,ob%intensity ! For testing
          ELSE
             ob%uindex    = nuniq
             ob%setnum    = iset
             ob%intensity = fsq
!!$           WRITE(21,'(2I10,F12.1)') ob%setnum,ob%uindex,ob%intensity ! For testing
          END IF
          ALLOCATE(ob%nextuniq); ob => ob%nextuniq; NULLIFY(ob%nextuniq)
       END IF
       IF (last) EXIT  ! Exit from the loop on hit "!END_OF_DATA"
       oldhkl=ihkl
    END DO  ! Loop over all reflections in the XSCALE.HKL file
    PRINT*, 'READING COMPLETE.'
    PRINT*,'# obs (excluding misfits), # unique =',nobs,nuniq

    PRINT*,'ALLOCATING ARRAYS FOR DATA...'
    IF(nobs > 0)THEN
       ALLOCATE(uindex(nobs)); ALLOCATE(setnum(nobs)); ALLOCATE(intensity(nobs))
       ALLOCATE(x(nuniq)); ALLOCATE(y(nuniq))
       ALLOCATE(sumI(nuniq,nset)); ALLOCATE(nfreq(nuniq,nset))
    END IF
    PRINT*,'ALLOCATED'

!!!!!!!!! Initializing arrays !!!!!!!!!!!!!!!!!!!!!
    DO i = 1, nobs
       uindex(i) = 0
       setnum(i) = 0
       intensity(i) = 0.
    END DO
    print*, 'initialization 1 passe'
    DO j = 1, nuniq
       x(j) = 0.
       y(j) = 0.
       DO i = 1, nset
          sumI(j,i) = 0.
          nfreq(j,i) = 0
       END DO
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    PRINT*, 'TRANSFERRING DATA TO ARRAYS...'
    i = 0
    ob => firstuniq
    DO WHILE(ASSOCIATED(ob%nextuniq))
       i = i + 1
       uindex(i) = ob%uindex
       setnum(i) = ob%setnum
       sumI(uindex(i),setnum(i)) = sumI(uindex(i),setnum(i)) + ob%intensity
       nfreq(uindex(i),setnum(i)) = nfreq(uindex(i),setnum(i)) + 1
       ob => ob%nextuniq
    END DO
    ! I assume that I want to calculate the CC between set_ii and the others. In this case
    ! X values are calculated using set ii and Y values are calculated using all other sets.

    ! Select the suitable subset of unique reflections from X and Y populations
    n = 0
    npairs = 0  ! This variable carries the number of pairs between two subsets (X and Y)
    SE = 0.

    PRINT*, 'CALCULATING CC. PLEASE WAIT...'
    CALL cpu_time(start)
    DO ii = 1, nset
       DO j = 1, nuniq
          IF (nfreq(j,ii) > 0) THEN
             nx = nx +1  ! Not used
             DO i = 1, nset
                IF ((i /= ii).AND.(nfreq(j,i) > 0)) THEN
                   ny = ny + nfreq(j,i)  
                   sumyy = sumyy + sumI(j,i)    
                END IF
             END DO
             IF(ny /= 0)THEN
                npairs = npairs + 1
                x(npairs) = sumI(j,ii)/nfreq(j,ii) 
                y(npairs) = sumyy/ny ! Kay method of averaging
!!$              WRITE(25,'(2I5,2ES15.5,I10)') ii,j,x(npairs),y(npairs),npairs ! Fort testing
             END IF
          END IF
          ny = 0
          sumyy = 0.
       END DO  ! Loop over all unique reflections

       ! Calculation of correlation coefficeint CC between two sets
       ! CorrelationCoefficient = [{SUM(XY)-(SUMX.SUMY)/n}/sqrt[{SUMX^2-(SUMX)^2/n}{SUMY^2-(SUMY)^2/n}]
       sumX = 0.
       sumY = 0.
       sumXY = 0.
       sumX2 = 0.
       sumY2 = 0.
       AA = 0.
       BB = 0.
       CC = 0.
       DD = 0.
       CorC = 0. ! Correlation Coefficeint
       DO i = 1, npairs
          sumX = sumX + x(i)        ! SUM(Xi)
          sumY = sumY + y(i)        ! SUM(Yi)
          sumXY = sumXY + x(i)*y(i) ! SUM(XiYi)

          sumX2 = sumX2 + x(i)**2   ! SUM(X^2)
          sumY2 = sumY2 + y(i)**2   ! SUM(Y^2)
       END DO
       AA = sumXY - (sumX * sumY)/npairs
       BB = sumX2 - ((sumX)**2)/npairs
       CC = sumY2 - ((sumY)**2)/npairs
       DD = sqrt(BB*CC)
       CorC = AA/DD
       CC2 = (CorC**2)
       SE = sqrt((1-CC2)/(npairs-2))
       dataset_i(ii) = ii
       cc_i(ii) = CorC
       WRITE(22,'(I5,2F8.3,I10)') ii,CorC,SE,npairs ! Writing to cc_datasets.txt
       WRITE(30,'(I5,2F8.3,I10)') ii,CorC,SE,npairs
       npairs = 0
    END DO  ! Loop over all datasets
    CALL cpu_time(finish)
    PRINT*, ' CPU time for CC calculation=',finish-start,'  seconds'
    CLOSE(unit); CLOSE(22), CLOSE(30)
    call plot(dataset_i,CorC_i)
  end subroutine cc_set

  subroutine dataset_remover(filename_template,filename_template_ccd,FDL,nset,nsetccd,xscale_ccd_template)
    IMPLICIT NONE
    character*255,intent(in) :: filename_template,filename_template_ccd
    integer,intent(in) :: nset
    character*10, intent(in) :: FDL
    integer,intent(out) :: nsetccd
    CHARACTER*255, intent(out) :: xscale_ccd_template
    CHARACTER (LEN=132) :: line=''
    character*3 :: strng
    INTEGER :: ier1,ier2,iset,np,i,j,cutoff2,cc2
    REAL cutoff,cc,err

1000 FORMAT(A)
1001 FORMAT(A,A)
1002 FORMAT(A,A,A)

    nsetccd = 0
!!$    PRINT*, 'How many datasets do you have?'
!!$    READ(*,*) nset
    print*
    print*, trim(adjustl(filename_template))
    print*, trim(adjustl(filename_template_ccd))
    print*, nset
    print*
    PRINT*, 'What is your CC cutoff (0 - 1)?'
    READ(*,*) cutoff
    cutoff2 = int(cutoff*100)
    WRITE(strng,'(I3)') cutoff2
    PRINT*, cutoff2,strng
    xscale_ccd_template = trim(adjustl(filename_template_ccd))//trim(adjustl(strng))
    print*, trim(adjustl(xscale_ccd_template))

    OPEN(UNIT=10,FILE=trim(adjustl(filename_template_ccd))//'.txt',STATUS='old',ACTION='READ',iostat=ier1)
    OPEN(UNIT=11,FILE=trim(adjustl(filename_template_ccd))//'.INP',STATUS='old',ACTION='READ',iostat=ier2)
    OPEN(UNIT=20,FILE=trim(adjustl(xscale_ccd_template))//'.INP',STATUS='replace',ACTION='WRITE')

    WRITE(20,1002) 'OUTPUT_FILE=',trim(adjustl(xscale_ccd_template))//'.ahkl'
    WRITE(20,1001) 'FRIEDEL''S_LAW=',FDL
    WRITE(20,1000) 'SAVE_CORRECTION_IMAGES=FALSE'
    IF (ier1/=0 .OR. ier2/=0) THEN
       PRINT*,'ier1=',ier1
       PRINT*,'ier2=',ier2
       STOP 'could not open file'
    END IF
    DO i = 1, 4
       READ(10,*) line
    END DO
    DO i = 1, 3
       READ(11,*) line
    END DO

    DO i = 1, nset
       READ(10,*) iset, cc, err, np
       IF(cc > cutoff)THEN
          nsetccd = nsetccd + 1
          DO j = 1, 3
             READ(11,1000) line
             WRITE(20,1000) line
          END DO
       ELSE
          DO j = 1, 3
             READ(11,1000) line
          END DO
       END IF
    END DO
    CLOSE(10);CLOSE(11);CLOSE(20)

  end subroutine dataset_remover

  subroutine cc_frames(xscale_isa_hkl,resol_ll,resol_ul)

    USE symops_mod
    USE resolution_mod
    !USE rc_mod
    !USE hash_mod

    IMPLICIT NONE

    character*255, intent(in) :: xscale_isa_hkl
    real, intent(in)   :: resol_ll,resol_ul
    integer :: unit,ier1

    INTEGER, PARAMETER  :: maxref=50000
    INTEGER ih,ik,il,oldhkl,jh,jk,jl,ihkl,iop,i,j,k,n,ispgr,ncheck, &
         ier2,iset,nset,nuniq,nobs,oldset,oldnuniq,ii,i1,i2,frm2,oldfrm2, &
         jj,l,npairs,nx,ny,nfrms,max_ii,ncc,is,olduindx,setn,oldframe
    INTEGER :: ifq(maxref)
    REAL s2,resol,fsq,sigma,dum2(2),frm,dum4,sumCC,start,finish
    CHARACTER (LEN=80) :: line
    REAL :: sumxx,sumyy,sumX,sumY,sumXY,sumX2,sumY2,AA,BB,CC,DD,CC2
    LOGICAL last,present

    TYPE :: dataset
       integer iiset
       TYPE(dataset),POINTER :: nextset
    END TYPE dataset

    TYPE :: dataframe
       integer iframe
       TYPE(dataframe),POINTER :: nextframe
    END TYPE dataframe

    TYPE :: dataint
       real iintensity
       TYPE(dataint),POINTER :: nextint
    END TYPE dataint

    TYPE ::uniq
       integer uindex,ufrq
       type(dataset) setnum
       type(dataint) intensity
       type(dataframe) framenum
       TYPE(uniq), POINTER :: nextuniq
    END TYPE uniq

    TYPE(uniq), POINTER :: ob, firstuniq
    TYPE(dataset), POINTER :: setnum, firstset
    TYPE(dataframe), POINTER :: framenum, firstframe
    TYPE(dataint), POINTER :: intensity, firstint

    INTEGER, DIMENSION(:), ALLOCATABLE :: set,uindex,ufrq,iiset,iframe
    REAL, DIMENSION(:), ALLOCATABLE :: iintensity,xx,yy

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: dset,dframe,np
    REAL, DIMENSION(:,:), ALLOCATABLE :: dint,CorC,SE

1000 FORMAT(A)

    !CALL hash_init(2*maxref)
    ncheck=0
    nset = 0

    unit=10
    OPEN(UNIT=unit,FILE=xscale_isa_hkl,STATUS='old',ACTION='READ',iostat=ier1)
    IF (ier1/=0) THEN
       PRINT*,'ier=',ier1
       STOP 'could not open file '!,xscale_isa_hkl
    END IF

    OPEN(UNIT=28,FILE='cc_frames.txt',STATUS='new')
    WRITE(28,*) '1st column- dataset number'
    WRITE(28,*) '2nd column- frame number'
    WRITE(28,*) '3rd column- CC between the frame and the rest of all data'
    WRITE(28,*) '4th column- number of pairs'
    WRITE(28,*)

    OPEN(UNIT=29,FILE='ccf_averaged.txt',STATUS='new')

    npairs = 0
    nx = 0
    ny = 0
    sumX = 0.
    sumY = 0.
    sumXY = 0.
    sumX2 = 0.
    sumY2 = 0.
    AA = 0.
    BB = 0.
    CC = 0.
    DD = 0.
    nuniq = 0
    last=.FALSE.
    present = .FALSE.
    ihkl=0
    i = 0
    ii = 1
    jj = 1
    max_ii = 1 ! This variable carrys the number of times a given uniq reflection be found in the dataset.

    ! Open input file; get spacegroup and cell constants
    DO
       READ(unit,1000,iostat=ier2) line
       IF (ier2/=0) THEN
          PRINT*,'ier=',ier2
          PRINT*,'line=',line(:LEN_TRIM(line))
          STOP 'read error in xds_file'
       END IF
       IF (line(1:20)=='!SPACE_GROUP_NUMBER=') then
          READ(line(21:),*) ispgr
          CALL getops(ispgr,op)
          WRITE(*,'(a)'),line(:LEN_TRIM(line))
       ELSE IF (line(1:21).eq.'!UNIT_CELL_CONSTANTS=') then 
          READ(line(22:),*) a,b,c,alf,bet,gam
          WRITE(*,'(a)'),line(:LEN_TRIM(line))
       ELSE IF (line(1:9) == '!ITEM_H=1') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:9) == '!ITEM_K=2') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:9) == '!ITEM_L=3') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:12) == '!ITEM_IOBS=4') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:19) == '!ITEM_SIGMA(IOBS)=5') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:13) == '!ITEM_ISET=10') THEN
          ncheck=ncheck+1
       ELSE IF (line(1:7) == '! ISET=') THEN
          READ(line(8:),*) iset
          nset=max(nset,iset)
       ELSE IF (line(1:14) == '!END_OF_HEADER') THEN
          EXIT
       END IF
    END DO
    IF (ncheck /= 6) stop 'header items not found'
    print*,'nset=',nset

    ALLOCATE(set(nset))
    ALLOCATE(ob); NULLIFY(ob%nextuniq)
    firstuniq => ob; nobs = 0
    ALLOCATE(setnum); NULLIFY(setnum%nextset)
    firstset => setnum
    ALLOCATE(framenum); NULLIFY(framenum%nextframe)
    firstframe => framenum
    ALLOCATE(intensity); NULLIFY(intensity%nextint)
    firstint => intensity

    PRINT*, 'READING DATA. PLEASE WAIT...'

    DO
       READ(unit,1000) line
       IF (line(1:12) == '!END_OF_DATA') last=.TRUE. 
       IF (.NOT.last) THEN
          READ(line,*)ih,ik,il,fsq,sigma,dum2,frm,dum4,iset
          If(abs(frm-NINT(frm))<0.05) CYCLE  ! After discussion with Kay
          frm2 = INT(frm)+1                  ! After discussion with Kay
          IF (sigma<0.) CYCLE                ! Filter out reflections with sigma <0
          CALL getres(ih,ik,il,s2,resol)
          IF (resol>resol_ll.OR.resol<resol_ul) CYCLE  ! Filter out reflections not within given resolution limit
          nobs=nobs+1  ! Counts the number of observations
          CALL asuput(ih,ik,il,jh,jk,jl,ihkl,iop,op,size(op,3)/ncenop) ! Generate a unique identifier (ihkl) for each reflection.
!!$        write(21,'(4I10,F12.1,2I10)') jh,jk,jl,ihkl,fsq,frm2,iset ! For testing
          IF (nobs==1) THEN
             oldhkl=ihkl
             nuniq = 1
             oldset = iset
          END IF
          IF (ihkl == oldhkl) THEN
             IF (nobs == 1) ii = 0
             i = nuniq
             ii = ii + jj
             IF(ii > max_ii) max_ii = ii
             ifq(i) = ii ! Choice for better performance
             ob%uindex = nuniq
             !           ob%ufrq = ii  
             ob%setnum%iiset = iset
             ob%framenum%iframe = frm2
             ob%intensity%iintensity = fsq
!!$           WRITE(22,*) i,ii,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
          ELSE IF (ihkl /= oldhkl) THEN
             nuniq = nuniq+1
             j = nuniq
             jj = 1
             ifq(j) = jj ! Choice for better performance         
             !           ob%ufrq = jj
             ob%uindex = nuniq
             ob%setnum%iiset = iset
             ob%framenum%iframe = frm2
             ob%intensity%iintensity = fsq
!!$           WRITE(22,*) j,jj,ob%uindex,ob%setnum%iiset,ob%framenum%iframe,ob%intensity%iintensity !For testing
             ii = 1
          END IF
          ALLOCATE(setnum%nextset); setnum => setnum%nextset; NULLIFY(setnum%nextset)
          ALLOCATE(framenum%nextframe); framenum => framenum%nextframe; NULLIFY(framenum%nextframe)
          ALLOCATE(intensity%nextint); intensity => intensity%nextint; NULLIFY(intensity%nextint)
          ALLOCATE(ob%nextuniq); ob => ob%nextuniq; NULLIFY(ob%nextuniq)
       END IF
       IF (last) EXIT  !Exit from the loop on hit "!END_OF_DATA"
       oldhkl=ihkl
       oldset = iset
    END DO  !Loop over all reflections in the XSCALE.HKL file
    PRINT*, 'READING COMPLETE.'
    PRINT*,'max_ii=',max_ii !For testing
    PRINT*,'# obs (excluding misfits), # unique =',nobs,nuniq

    PRINT*,'ALLOCATING ARRAYS FOR DATA...'
    IF(nobs > 0)THEN
       ALLOCATE(uindex(nuniq)); ALLOCATE(iiset(nobs)); ALLOCATE(iframe(nobs))
       ALLOCATE(xx(nobs)); ALLOCATE(yy(nobs))
       ALLOCATE(dset(nuniq,max_ii)); ALLOCATE(dframe(nuniq,max_ii)); ALLOCATE(dint(nuniq,max_ii)) 
    END IF
    PRINT*,'ALLOCATED'

    PRINT*, 'TRANSFERRING DATA TO ARRAYS...'
    i = 0
    j = 0
    ob => firstuniq
    setnum => firstset
    framenum => firstframe
    intensity => firstint
    DO WHILE(ASSOCIATED(ob%nextuniq))
       i = i + 1
       iiset(i)  = ob%setnum%iiset
       iframe(i) = ob%framenum%iframe
       IF(i == 1) olduindx = ob%uindex
       IF(olduindx == ob%uindex)THEN
          j = j + 1
          k = ob%uindex
          uindex(k)   = k
          dset(k,j)   = ob%setnum%iiset
          dframe(k,j) = ob%framenum%iframe
          dint(k,j)   = ob%intensity%iintensity
!!$        WRITE(23,*) i,j,uindex(k),dset(k,j),dframe(k,j),dint(k,j) !For testing
       ELSEIF(olduindx /= ob%uindex)THEN
          j = 1
          k = ob%uindex
          uindex(k)   = k
          dset(k,j)   = ob%setnum%iiset
          dframe(k,j) = ob%framenum%iframe
          dint(k,j)   = ob%intensity%iintensity
!!$        WRITE(23,*) i,j,uindex(k),dset(k,j),dframe(k,j),dint(k,j) !For testing
       END IF
       olduindx = ob%uindex
       setnum => setnum%nextset
       framenum => framenum%nextframe
       intensity => intensity%nextint
       ob => ob%nextuniq
    END DO

    oldframe = 0
    DO is = 1, nset
       DO j = 1, nobs
          IF(is == iiset(j))THEN
             IF(oldframe < iframe(j)) oldframe = iframe(j)
          END IF
       END DO
       set(is) = oldframe
       WRITE(24,*) is, set(is)
       oldframe = 0
    END DO

    nfrms = MAXVAL(set); PRINT*,'Highest frame-number of all datasets=', nfrms

    ALLOCATE(np(nfrms,nset)); ALLOCATE(CorC(nfrms,nset));ALLOCATE(SE(nfrms,nset))
    DO i = 1, nset
       DO j = 1, nfrms
          CorC(j,i) = 0.
          SE(j,i)   = 0.
          np(j,i)   = 0
       END DO
    END DO

    PRINT*, 'CALCULATING CC. PLEASE WAIT...'
    CALL cpu_time(start)
    DO i = 1, nset
       DO j = 1, set(i)  !Maximum number of frames in ith dataset
          DO k = 1, nuniq
             DO l = 1, ifq(k)  !Occurence of kth unique reflection in all data
                IF((dset(k,l)==i .AND. dframe(k,l)==j) .AND. (dint(k,l) /= 0))THEN
                   nx = nx + 1
                   sumxx = sumxx + dint(k,l)
                   present = .TRUE.
                END IF
             END DO

             IF(present)THEN
                DO l = 1, ifq(k)  !Occurence of kth unique reflection in all data
                   IF((.NOT.(dset(k,l)==i .AND. dframe(k,l)==j)).AND.(dint(k,l) /= 0))THEN 
                      ny = ny + 1
                      sumyy = sumyy + dint(k,l)
                   END IF
                END DO
                IF(ny /= 0)THEN
                   xx(k) = sumxx/nx
                   yy(k) = sumyy/ny
                   npairs = npairs + 1

!!$                 WRITE(25,*) nx,ny,npairs,i,j,k,xx(k),yy(k) !For testing
!!$                 WRITE(26,*) nx,sumxx,xx(k),ny,sumyy,yy(k)  !For testing

                   sumX = sumX + xx(k)         ! SUM(Xi)
                   sumY = sumY + yy(k)         ! SUM(Yi)
                   sumXY = sumXY + xx(k)*yy(k) ! SUM(XiYi)

                   sumX2 = sumX2 + xx(k)**2    ! SUM(X^2)
                   sumY2 = sumY2 + yy(k)**2    ! SUM(Y^2)
                   ! PRINT*, sumX,sumY,sumXY,sumX2,sumY2
                END IF
             END IF
             present = .FALSE.
             nx = 0
             ny = 0
             sumxx = 0.0
             sumyy = 0.0
          END DO
          ! Here do the CC calculation.
          ! PRINT*, npairs
          IF(sumX /= 0 .OR. sumY /= 0)THEN
             AA = sumXY - (sumX * sumY)/npairs
             BB = sumX2 - ((sumX)**2)/npairs
             CC = sumY2 - ((sumY)**2)/npairs
             DD = sqrt(BB*CC)
             CorC(j,i) = AA/DD
             np(j,i) = npairs
!!$           WRITE(27,*) sumX,sumY,sumXY,sumX2,sumY2,AA,DD,CorC(j,i),np(j,i),i,j !For testing
             npairs = 0
             sumX = 0.
             sumY = 0.
             sumXY = 0.
             sumX2 = 0.
             sumY2 = 0.
          END IF
          ! Standard error in CC
          CC2 = (CorC(j,i)**2)
          SE(j,i) = sqrt((1-CC2)/(np(j,i)-2))
          IF(np(j,i) > 3) WRITE(6,'(2I5,2F8.3,I5)') i,j,CorC(j,i),SE(j,i),np(j,i) !Screen output
       END DO
    END DO
    CALL cpu_time(finish)
    PRINT*, ' CPU time =',finish-start,'  seconds'

    PRINT*, 'PRINTING CC VALUES TO ARRAYS...'
    DO i = 1, nset
       DO j = 1, nfrms
          IF(np(j,i) < 3)THEN  ! Filter out not-defined situations for the CC calculation
             CorC(j,i) = 0.0
             SE(j,i)   = 0.0
             np(j,i)   = 0
          END IF
          WRITE(28,'(2I5,2F8.3,I5)') i,j,CorC(j,i),SE(j,i),np(j,i) !Writing to cc_frames.txt
       END DO
    END DO

    ! Here calculates the averaged CC for each frame aveaged over all datasets.
    sumCC = 0.
    ncc = 0  
    DO j = 1, nfrms
       DO i = 1, nset
          IF(CorC(j,i) /= 0)THEN
             ncc = ncc + 1
             sumCC = sumCC + CorC(j,i)
          END IF
       END DO
       IF(ncc /= 0)THEN
          WRITE(29,'(I5,F8.3)') j,sumCC/ncc !Writing to cc_averaged.txt
          sumCC = 0.
          ncc = 0
       END IF
    END DO

    PRINT*, 'EVERYTHING IS SUCCESSFUL.' 
    PRINT*, 'Please look at cc_frames.txt AND cc_averaged.txt in your working directory. Good luck!'
    CLOSE(unit);CLOSE(24);CLOSE(28);CLOSE(29)
  end subroutine cc_frames

!!$  subroutine frame_remover(xscale_inp,ccf_remove_in)
!!$    IMPLICIT NONE
!!$    character*255,intent(in) :: xscale_inp,ccf_remove_in
!!$    ! eventually ccf_remove_in variable should carry the filename of frames_remove.txt
!!$    ! and
!!$    ! xscale_inp should carry the filename of XSCALE_ISa_ccd.INP
!!$    INTEGER :: unit,ier,nset,i,iset,f1,f2,j
!!$    CHARACTER (LEN=132) :: line
!!$
!!$    nset=100 !hard coded
!!$    unit=10
!!$
!!$    OPEN(UNIT=11,FILE='frames_remove.txt',STATUS='old',ACTION='READ')
!!$    OPEN(UNIT=12,FILE='XSCALE_cc_cutoff0p9.INP',STATUS='old',ACTION='READ')
!!$    OPEN(UNIT=20,FILE=xscale_in,STATUS='new',ACTION='WRITE')
!!$    OPEN(UNIT=21,FILE=xscale_log,STATUS='new',ACTION='WRITE')
!!$
!!$    DO i=1,3
!!$       READ(12,*) line
!!$       WRITE(20,*) line
!!$    END DO
!!$    DO i=1,nset
!!$       READ(11,*) iset,f1,f2
!!$       IF(i==1)THEN
!!$          READ(12,'(a)') line
!!$       ELSE
!!$          DO j=1,3
!!$             READ(12,'(a)') line
!!$          END DO
!!$       END IF
!!$       IF((f2-f1)==24)THEN
!!$          WRITE(21,*) 'Dataset # =',iset,' removed'
!!$          CYCLE
!!$       ELSEIF(f2==0)THEN
!!$          WRITE(20,*) line
!!$          WRITE(20,*) 'INCLUDE_RESOLUTION_RANGE=30.0 2.8'
!!$          WRITE(20,*) 'MINIMUM_I/SIGMA=3'
!!$          WRITE(21,*) 'Dataset # =',iset,' put in; no modification'
!!$       ELSE
!!$          CALL ascii_file(unit,iset,f1,f2,line)
!!$          CLOSE(unit)
!!$       END IF
!!$    END DO
!!$
!!$
!!$  end subroutine frame_remover
  
END MODULE SUB_ROUTINES

