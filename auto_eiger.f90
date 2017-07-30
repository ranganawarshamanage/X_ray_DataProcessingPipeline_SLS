!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   MAIN PROGRAM              !!!!!!!!
!!!!!!!           Author: Rangana Warshamanage      !!!!!!!!
!!!!!!!                 Date: 22.09.2015            !!!!!!!!
!!!!!!!           Modified date: 29.02.2016         !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE copyright_module
  CHARACTER (LEN=*), PARAMETER :: copyright_text = '(&
       &"________________________________________________________________________________"/&
       &"                                                                                "/&
       &"                 M X Automator                                                  "/&
       &"                 _______________________________________________________________"/&
       &"                 version 1.0 / February 2016                                    "/&
       &"                                                                                "/&
       &"                 Copyright (c) Rangana Warshamanage, PSI, Switzerland, 2016.    "/&
       &"                 All rights reserved.                                           "/&
       &"                                                                                "/&
       &"                 Data processing pipeline for multi-crystals using              "/&
       &"                 XDS and rotation data collected on PILATUS 6M nad EIGER 16M    "/&
       &"                 at beamlines X06SA/X10SA at the Paul Scherrer Institut.        "/&
       &"________________________________________________________________________________"/&
       &) '
END MODULE copyright_module

MODULE parameter_module

  CHARACTER (LEN=1) :: AA,isa_again,response
  CHARACTER (LEN=10) :: FDL,ABSOPT,strng1,strng2,xtal,dat_fmt
  CHARACTER (LEN=50) :: ref_dataset
  CHARACTER(LEN=255) :: data_dir,dir_well,link,link_dir_well,img_path,XSCALE_local,&
       XGBL_path,XLCL_path,G_path,head_file,cwd_path,line_inp,&
       xscale_isa_template,filename_template,filename_template_ccd,&
       xscale_ccd_template
  CHARACTER, DIMENSION(:), ALLOCATABLE :: strng3
  INTEGER, DIMENSION(:), ALLOCATABLE :: mps
  INTEGER :: spgr,length,i,nn,j,run,nset,nsetccd,nwell,&
       strng_pxl,min_pixel_spot,xl_n,npxl_lines,STAT1,STAT2
  REAL :: a,b,c,alf,bet,gam,resol_ll,resol_ul,frac_index,tot_rotation_angle

  CHARACTER (LEN=20) :: pathfile = 'pathfile.txt'
  CHARACTER (LEN=20) :: inpfile  = 'auto.inp'
  CHARACTER (LEN=20) :: pxlfile = 'pxlfile.txt'
  CHARACTER (LEN=20) :: XSCALE_global= 'XSCALE_global.INP'
  LOGICAL :: pathfile_exists = .FALSE.
  LOGICAL :: inpfile_exists  = .FALSE.
  LOGICAL :: pxlfile_exists  = .FALSE.
  LOGICAL :: file1_exists    = .FALSE.
  LOGICAL :: file2_exists    = .FALSE.
  LOGICAL :: file3_exists    = .FALSE.
  LOGICAL :: linkcorr_exists = .FALSE.
  LOGICAL :: xscgl_exists    = .FALSE.
  LOGICAL :: REFDATA         = .FALSE.
  LOGICAL :: dir_e  = .FALSE.
  LOGICAL :: dir_e2 = .FALSE.
  LOGICAL :: dir_e3 = .FALSE.
  INTEGER :: nxtals = 0
  INTEGER :: unitpxl = 12
  INTEGER :: imps = 0

END MODULE parameter_module

MODULE global_module

  USE SUB_ROUTINES
  USE copyright_module
  USE parameter_module
  USE gnufor2

END MODULE global_module


PROGRAM automation

  USE global_module
  IMPLICIT NONE

1000 FORMAT(A)
1001 FORMAT(A,A)
1005 FORMAT(A,F4.1,A,F4.1)


  WRITE(*,copyright_text)

  call getcwd(cwd_path) !Get the path of working directory

  PRINT*
  WRITE(6,*) 'YOUR WORKING DIRECTORY:'
  PRINT*, trim(cwd_path)

  call file_check(pathfile); call file_check(inpfile)
  call pxlfile_check(pxlfile,pxlfile_exists,mps)

  WRITE(6,*) 'Enter your name of the data processing directory?'
  READ(5,*) data_dir

  call read_input(pxlfile_exists)

  call file_open('linkcorrectfile.txt')
  call file_open(XSCALE_global,FDL,ABSOPT)


  CALL SYSTEM("wc -l < " // pathfile // " > nwell.tmp")
  call system_var('nwell.tmp',nwell)

  CALL SYSTEM("rm nwell.tmp"); CALL SYSTEM("rm img_*")

  OPEN(UNIT=10,FILE=pathfile,STATUS='old')
  DO nn = 1,nwell 
     READ(10,1000) dir_well      !Name of the current well
     print*, dir_well
     dir_well = TRIM(ADJUSTL(dir_well))
     link_dir_well = TRIM(ADJUSTL(dir_well))//TRIM(ADJUSTL(data_dir))
     link_dir_well = TRIM(ADJUSTL(link_dir_well))
     print*, link_dir_well
     WRITE(strng1,'(I5)') nn
     strng1=TRIM(ADJUSTL(strng1))
     XGBL_path = 'img_'//strng1
     XGBL_path = TRIM(ADJUSTL(XGBL_path))
     print*, XGBL_path
     CALL SYSTEM("ln -s " // link_dir_well // " " // XGBL_path)
     CALL SYSTEM("pwd")
     CALL CHDIR(dir_well,STAT1)
     IF(STAT1 /= 0)THEN
        WRITE(6,*) 'navigtion to ', dir_well,' FAILED'
        STOP
     ELSE
        CALL SYSTEM("ln -s " // cwd_path // "wd_path") ! 25.01.2016      
        if(trim(adjustl(dat_fmt))=='H5')then
           call system('ls *_master.h5 > master.name')
           open(unit=10,file='master.name',status='old')
           read(10,'(a)') head_file
           close(10)
           call system("rm master.name")

           call system('h5header.sh '//trim(head_file)//' > h5header.tmp')
           call sleep(2)
           call system('ls -l *_data_*.h5 | wc -l > nxtals.tmp ')
           call system_var('nxtals.tmp',nxtals)
        end if

        dir_e2 = .FALSE.
        INQUIRE(FILE=trim(adjustl(data_dir))//'/.',EXIST=dir_e2)

        IF(dir_e2)THEN ! if data_dir already exists
           CALL CHDIR(trim(adjustl(data_dir)),STAT2)
           IF(STAT2 /= 0)THEN
              WRITE(6,*) 'navigation to ',trim(adjustl(data_dir)),' folder FAILED'
              STOP
           ELSE
              CALL SYSTEM("rm slink")
              CALL SYSTEM("ln -s " // dir_well // "slink") 
              call file_open('XSCALE_local.INP',FDL,ABSOPT)
              i = 0;

              DO
                 i = i + 1
                 if(trim(adjustl(dat_fmt))=='CBF')then
                    IF((i==1).OR.(i==nxtals))THEN
                       CALL SYSTEM("ls -l " // dir_well // " | grep ^d | grep xtal_ | wc -l > xx.temp")
                       call system_var('xx.temp',nxtals)
                    END IF
                    print*, i
                 end if

                 IF(i<=nxtals)THEN
!!$              IF(pxlfile_exists)THEN
!!$                 imps = imps + 1
!!$                 min_pixel_spot = mps(imps)
!!$              END IF
                    WRITE(strng2,'(I6)') i
                    strng2 =TRIM( ADJUSTL(strng2))
                    xtal='xtal_'//strng2             !xtal_?????
                    xtal = trim(adjustl(xtal))
                    print*, xtal                     !for testing
                    dir_e3 = .FALSE.
                    INQUIRE(FILE=trim(adjustl(xtal))//'/.',EXIST=dir_e3) 
                    IF(dir_e3)THEN                   !if xtal_n exists
                       CALL CHDIR(xtal)
                       file1_exists = .FALSE.
                       INQUIRE(FILE='GXPARM.XDS',EXIST=file1_exists)
                       IF(file1_exists)THEN
                          WRITE(6,*) xtal, 'indexing and integration has been successful'
                          write(41,*)  xtal, 'indexing and integration has been successful'
                          XLCL_path = '../'//trim(xtal)//'/XDS_ASCII.HKL'
                          G_path = trim(XGBL_path)//'/'//trim(xtal)//'/XDS_ASCII.HKL'
                          call file_write(21,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                          call file_write(22,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                          call file_write(23,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                          PRINT*,trim(XGBL_path)//'/'//trim(xtal)//'/CORRECT.LP'
                       END IF
                       CALL CHDIR("../")
                       CYCLE
                    ELSE
                       write(40,1000) xtal
                       print*, LEN_TRIM(xtal) 
                       CALL SYSTEM("mkdir " // xtal)
                       CALL CHDIR(xtal)
                       CALL SYSTEM("pwd")
                       if(trim(adjustl(dat_fmt))=='CBF')then
                          head_file = '../../'//trim(xtal)//'/*00001.cbf'
                          print*, head_file
                          CALL SYSTEM("head -35 " // trim(head_file) // " > frame1.head")
                       end if
                       do run = 1, 2
                          CALL createinputfile(i,dat_fmt,xtal,run,tot_rotation_angle,a,b,c,alf,bet,gam,spgr,FDL,ABSOPT,&
                               resol_ll,resol_ul,strng_pxl,min_pixel_spot,frac_index,REFDATA,ref_dataset,head_file)
                          CALL SYSTEM("xds_par XDS.INP")
                          if(run==1) CALL SYSTEM("mv XDS.INP XDS_IDX.INP")
                       end do
                       file2_exists = .FALSE.
                       INQUIRE(FILE='GXPARM.XDS',EXIST=file2_exists)
                       IF(file2_exists)THEN
                          WRITE(6,*) xtal, 'indexing and integration is successful'
                          write(41,*)  xtal, 'indexing and integration is successful'
                          XLCL_path = '../'//trim(xtal)//'/XDS_ASCII.HKL'
                          G_path = trim(XGBL_path)//'/'//trim(xtal)//'/XDS_ASCII.HKL'
                          call file_write(21,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                          call file_write(22,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                          call file_write(23,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                       ELSE
                          WRITE(6,*) xtal, 'indexing and integration failed'
                       END IF
                       CALL CHDIR("../") !Back to automation directory
                    END IF
                    CYCLE
                 ELSE
                    EXIT
                 END IF
              END DO ! over nxtals
              CLOSE(21)
              dir_e = .FALSE.
              INQUIRE(FILE='./autoXSCALE/.',EXIST=dir_e)
              IF(dir_e) CALL SYSTEM("rm -r autoXSCALE")
              CALL SYSTEM("mkdir autoXSCALE")
              CALL CHDIR("autoXSCALE")
              CALL SYSTEM("cp ../XSCALE_local.INP XSCALE.INP" )
              CALL SYSTEM("xscale XSCALE.INP")
              CALL CHDIR(trim(cwd_path)) !Back to working directory
           END IF
        ELSE ! if 'data_dir' does not exist in the well.
           CALL SYSTEM("mkdir " // trim(adjustl(data_dir)))   !Directory for data processing in this well
           CALL CHDIR(trim(adjustl(data_dir)),STAT2)
           print*, STAT2
           call system("pwd")
           IF(STAT2 /= 0)THEN
              WRITE(6,*) 'navigation to ',trim(adjustl(data_dir)),' folder failed'
              STOP
           ELSE
	      CALL SYSTEM("rm slink")
              CALL SYSTEM("ln -s " // dir_well // "slink")
              call file_open('XSCALE_local.INP',FDL,ABSOPT)
              i = 0;
              DO
                 i = i + 1
                 if(trim(adjustl(dat_fmt))=='CBF')then
                    IF((i==1).OR.(i==nxtals))THEN
                       CALL SYSTEM("ls -l " // dir_well // " | grep ^d | grep xtal_ | wc -l > xx.temp")
                       call system_var('xx.temp',nxtals)
                    END IF
                    print*, i
                 end if

                 IF(i<=nxtals)THEN
                    IF(pxlfile_exists)THEN
                       imps = imps + 1
                       min_pixel_spot = mps(imps)
                    END IF
                    WRITE(strng2,'(I6)') i
                    strng2 =TRIM( ADJUSTL(strng2))
                    xtal='xtal_'//strng2 
                    xtal = trim(adjustl(xtal))
                    write(40,1000) xtal
                    print*, LEN_TRIM(xtal) 
                    CALL SYSTEM("mkdir " // xtal)
                    CALL CHDIR(xtal)
                    CALL SYSTEM("pwd")
                    if(trim(adjustl(dat_fmt))=='CBF')then
                       head_file = '../../'//trim(xtal)//'/*00001.cbf'
                       print*, head_file
                       CALL SYSTEM("head -35 " // trim(head_file) // " > frame1.head")
                    end if
                    do run = 1, 2
                       CALL createinputfile(i,dat_fmt,xtal,run,tot_rotation_angle,a,b,c,alf,bet,gam,spgr,FDL,ABSOPT,&
                            resol_ll,resol_ul,strng_pxl,min_pixel_spot,frac_index,REFDATA,ref_dataset,head_file)
                       CALL SYSTEM("xds_par XDS.INP")
                       if(run==1) CALL SYSTEM("mv XDS.INP XDS_IDX.INP")
                    end do
                    file3_exists = .FALSE.
                    INQUIRE(FILE='GXPARM.XDS',EXIST=file3_exists)
                    IF(file3_exists)THEN
                       WRITE(6,*) xtal, 'indexing and integration is successful'
                       write(41,*)  xtal, 'indexing and integration is successful'
                       XLCL_path = '../'//trim(xtal)//'/XDS_ASCII.HKL'
                       G_path = trim(XGBL_path)//'/'//trim(xtal)//'/XDS_ASCII.HKL'
                       call file_write(21,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                       call file_write(22,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                       call file_write(23,trim(XLCL_path),trim(G_path),resol_ll,resol_ul)
                    ELSE
                       WRITE(6,*) xtal, 'indexing and integration failed'
                    END IF
                    CALL CHDIR("../") !Back to automation directory
		    CYCLE
                 ELSE
                    EXIT
                 END IF
              END DO ! over nxtals
              CLOSE(21)
              CALL SYSTEM("mkdir autoXSCALE")
              CALL CHDIR("autoXSCALE")
              CALL SYSTEM("cp ../XSCALE_local.INP XSCALE.INP" )
              CALL SYSTEM("xscale XSCALE.INP")
              CALL CHDIR(trim(cwd_path)) !Back to working directory
           END IF
        END IF
     END IF
  END DO ! over nwells
  CLOSE(10);CLOSE(22);CLOSE(23)

  CALL CHDIR(trim(cwd_path))
  CALL ranker(FDL,ABSOPT,resol_ll,resol_ul)
  CALL SYSTEM("cp XSCALE_GBL_SRT.INP XSCALE.INP")
  CALL SYSTEM("xscale XSCALE.INP")
  CALL SYSTEM("cp XSCALE.LP XSCALE_GBL_SRT.LP")
  WRITE(6,*) 'Processing and scaline done!'
  WRITE(6,*)
  WRITE(6,*) 'Now starts dataset selection procedure'
  WRITE(6,*)
  WRITE(6,*) '***********************************'
  WRITE(6,*) 'Test 1). Based on dataset ISa value'
  WRITE(6,*) '***********************************'
  WRITE(6,*)
  WRITE(6,*) 'List of ISa values in XSCALE file'
  DO
     CALL isa_cutoff(FDL,ABSOPT,resol_ll,resol_ul,xscale_isa_template)
     WRITE(6,*) 'XSCALE with new input file...'
     CALL SYSTEM("cp " // trim(adjustl(xscale_isa_template)) // ".INP XSCALE.INP")
     CALL SYSTEM("xscale XSCALE.INP")
     CALL SYSTEM("cp XSCALE.LP " // trim(adjustl(xscale_isa_template)) // ".LP")
     WRITE(6,*) 'Do you want to use a different ISa-cutoff? (y/n)'
     READ(5,'(a)') isa_again
     IF(isa_again=='y')THEN
        ! CALL SYSTEM("cp XSCALE.LP " // trim(adjustl(xscale_isa_template)) // ".LP")
        CYCLE
     ELSE
        EXIT
     END IF
  END DO
  WRITE(6,*)
  WRITE(6,*) '**********************************'
  WRITE(6,*) 'Test 2). Based on dataset CC value'
  WRITE(6,*) '**********************************'
  WRITE(6,*)
  CALL cc_set(resol_ll,resol_ul,nset,filename_template,filename_template_ccd)
  CALL dataset_remover(filename_template,filename_template_ccd,FDL,nset,nsetccd,xscale_ccd_template)
  CALL SYSTEM("cp " // trim(adjustl(xscale_ccd_template)) // ".INP XSCALE.INP")
  CALL SYSTEM("xscale XSCALE.INP")
  CALL SYSTEM("cp XSCALE.LP " // trim(adjustl(xscale_ccd_template)) // ".LP" )
  WRITE(6,*)
  WRITE(6,*)
  WRITE(6,*) '************ ALL DONE ************'
  WRITE(6,*)


CONTAINS

  SUBROUTINE file_check(file_in)
    USE global_module

    CHARACTER (LEN=*), INTENT(IN) :: file_in
    LOGICAL :: presence

    INQUIRE(FILE=trim(adjustl(file_in)),EXIST=presence)
    IF(presence)THEN
       WRITE(*,*) file_in, 'exists'
    ELSE
       WRITE(6,*) file_in,' is MISSING'
       STOP
    END IF

  END SUBROUTINE file_check


  SUBROUTINE read_input(pxlfile_present)
    USE global_module

    LOGICAL, INTENT(IN) :: pxlfile_present

    OPEN(UNIT=10,FILE=inpfile,STATUS='OLD',ACTION='READ')
    DO
       READ(10,'(a)') line_inp
       line_inp = trim(adjustl(line_inp))
       do i=1, len_trim(line_inp)
          if(line_inp(i:i)=='!') exit
       end do
       IF(line_inp(1:7)=='FORMAT=') READ(line_inp(8:i),*) dat_fmt
       IF(line_inp(1:20)=='UNIT_CELL_CONSTANTS=') READ(line_inp(21:i),*) a,b,c,alf,bet,gam  
       IF(line_inp(1:20)=='SPACE_GROUP_NUMBER=') READ(line_inp(21:i),*) spgr
       IF(line_inp(1:7)=='FRIEDEL') READ(line_inp(16:i),*) FDL
       IF(line_inp(1:29)=='STRICT_ABSORPTION_CORRECTION=') READ(line_inp(30:i),*) ABSOPT
       IF(line_inp(1:25)=='INCLUDE_RESOLUTION_RANGE=') READ(line_inp(26:i),*) resol_ll,resol_ul
       IF(line_inp(1:13)=='STRONG_PIXEL=') READ(line_inp(14:i),*) strng_pxl
       IF(.not.pxlfile_present)THEN
          IF(line_inp(1:35)=='MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=') READ(line_inp(36:i),*) min_pixel_spot
       END IF
       IF(line_inp(1:34)=='MINIMUM_FRACTION_OF_INDEXED_SPOTS=') READ(line_inp(35:i),*) frac_index
       IF(line_inp(1:21)=='TOTAL_ROTATION_ANGLE=') READ(line_inp(22:i),*) tot_rotation_angle
       IF(line_inp(1:19)=='REFERENCE_DATA_SET=')THEN
          READ(line_inp(20:i),*) ref_dataset
          REFDATA = .TRUE.
       END IF
       IF(line_inp(1:14)=='AUTO_INPUT_END') EXIT
    END DO
    CLOSE(10)
    WRITE(*,*) 'FORMAT ',dat_fmt
    WRITE(*,*) 'UNIT_CELL_CONSTANTS=',a,b,c,alf,bet,gam
    WRITE(*,*) 'SPACE_GROUP_NUMBER=', spgr
    WRITE(*,*) 'FRIEDEL''S_LAW=', FDL
    WRITE(*,*) 'STRICT_ABSORPTION_CORRECTION=', ABSOPT
    WRITE(*,*) 'INCLUDE_RESOLUTION_RANGE=', resol_ll,resol_ul
    WRITE(*,*) 'STRONG_PIXEL=', strng_pxl
    IF(.not.pxlfile_present)THEN
       WRITE(*,*) 'MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=', min_pixel_spot
    END IF
    WRITE(*,*) 'MINIMUM_FRACTION_OF_INDEXED_SPOTS=', frac_index
    WRITE(*,*) 'TOTAL_ROTATION_ANGLE=', tot_rotation_angle
    IF(REFDATA) WRITE(*,*) 'REFERENCE_DATA_SET=', ref_dataset

    WRITE(6,*) 'Is everything alright (y/n)?'
    READ(5,*) response
    IF(((TRIM(ADJUSTL(response)))=='Y').OR.((TRIM(ADJUSTL(response)))=='y'))THEN
       continue
    ELSE
       stop 'Edit your auto.inp and rerun'
    END IF

  END SUBROUTINE read_input

  SUBROUTINE pxlfile_check(file_in,pxlfile_present,smps)
    USE global_module

1000 FORMAT(A)
1001 FORMAT(A,A)
1005 FORMAT(A,F4.1,A,F4.1)

    CHARACTER (LEN=*), INTENT(IN) :: file_in
    LOGICAL, INTENT(OUT) :: pxlfile_present
    INTEGER :: sxl_n,snpxl_lines,sunitpxl,si
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: smps

    pxlfile_present  = .FALSE.
    sunitpxl = 12

    INQUIRE(FILE=file_in,EXIST=pxlfile_present)
    IF(pxlfile_exists)THEN
       WRITE(*,*) file_in, 'exists'
       CALL SYSTEM("wc -l < " // file_in // " > pxl.tmp")
       call system_var('pxl.tmp',snpxl_lines)
!!$       OPEN(UNIT=10,FILE='pxl.tmp',STATUS='OLD')
!!$       READ(10,*) snpxl_lines
!!$       CLOSE(10)
!!$       CALL SYSTEM("rm pxl.tmp")
       ALLOCATE(smps(snpxl_lines-1))
       OPEN(UNIT=sunitpxl,FILE=file_in,STATUS='OLD')
       READ(sunitpxl,1000)
       do i = 1,snpxl_lines-1 
          READ(sunitpxl,*) sxl_n,smps(i)
       end do
       CLOSE(sunitpxl)
    ELSE
       WRITE(6,*) 'pxlfile.txt file is MISSING. auto.inp parameters used.'
    END IF

  END SUBROUTINE pxlfile_check


  SUBROUTINE file_open(file_in,sFDL,sABSOPT)
    USE global_module

1000 FORMAT(A)
1001 FORMAT(A,A)
1005 FORMAT(A,F4.1,A,F4.1)

    CHARACTER (LEN=*), INTENT(IN) :: file_in
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: sFDL,sABSOPT

    IF(file_in=='XSCALE_local.INP')THEN
       OPEN(UNIT=21,FILE=file_in,STATUS='REPLACE',ACTION='WRITE')
       WRITE(21,1000) 'OUTPUT_FILE=XSCALE_local.ahkl'
       WRITE(21,1001) 'FRIEDEL''S_LAW=',sFDL
       WRITE(21,1000) 'SAVE_CORRECTION_IMAGES=',sABSOPT
    END IF

    IF(file_in==XSCALE_global)THEN
       OPEN(UNIT=22,FILE=file_in,STATUS='REPLACE',ACTION='WRITE')
       WRITE(22,1000) 'OUTPUT_FILE=XSCALE_global.ahkl'
       WRITE(22,1001) 'FRIEDEL''S_LAW=',sFDL
       WRITE(22,1001) 'SAVE_CORRECTION_IMAGES=',sABSOPT
    END IF

    IF(file_in=='linkcorrectfile.txt')THEN
       OPEN(UNIT=23,FILE=file_in,STATUS='REPLACE',ACTION='WRITE')
    END IF

  END SUBROUTINE file_open


  SUBROUTINE file_write(file_unit,sXLCL_path,sG_path,sresol_ll,sresol_ul)
    USE global_module

1000 FORMAT(A)
1001 FORMAT(A,A)
1005 FORMAT(A,F4.1,A,F4.1)

    CHARACTER (LEN=*), INTENT(IN), OPTIONAL ::sXLCL_path,sG_path
    CHARACTER (LEN=255) :: string1
    INTEGER, INTENT(IN) :: file_unit
    INTEGER :: si
    REAL, INTENT(IN), OPTIONAL :: sresol_ll, sresol_ul

    IF(file_unit==21)THEN
       WRITE(21,1001) 'INPUT_FILE=',sXLCL_path
       WRITE(21,1005) 'INCLUDE_RESOLUTION_RANGE=',sresol_ll,' ',sresol_ul
       WRITE(21,1000) 'MINIMUM_I/SIGMA=0'
    END IF
    IF(file_unit==22)THEN
       WRITE(22,1001) 'INPUT_FILE=',trim(sG_path)
       WRITE(22,1005) 'INCLUDE_RESOLUTION_RANGE=',sresol_ll,' ',sresol_ul
       WRITE(22,1000) 'MINIMUM_I/SIGMA=0'
    END IF
    IF(file_unit==23)THEN
       si = len_trim(sG_path)-13
       read(sG_path(1:si),1000) string1
       WRITE(23,1001) trim(string1)//'CORRECT.LP'
    END IF

  END SUBROUTINE file_write

  subroutine system_var(strng_in,para_out)
    character(len=*),intent(in) :: strng_in
    integer,intent(out) :: para_out

    open(UNIT=10,FILE=strng_in,STATUS='old')
    read(10,*) para_out; close(10)
    call system("rm "//strng_in)
    
  end subroutine system_var


END PROGRAM automation

