!******************************************************
    subroutine CMS_Tools_Dialog    
! Perform one of a selection of internal tools
! written by Mitchell Brown, USACE-CHL
! 1 and 2 - 04/14/2021
! 3 - 9/30/2021
! 4 - 2/09/2022
! 5 - 3/16/2022
!******************************************************    
    use diag_lib, only: diag_print_message
    implicit none
    integer :: ichoice 
    
50	continue
    
#ifdef _WIN32        
    write(*,*) '*****************************************************************'
    write(*,*) 'Make a choice from the following CMS Tools'
    write(*,*) '1 - Print RT_JULIAN/REFTIME corresponding to actual date'
    write(*,*) '2 - Print actual date corresponding to RT_JULIAN/REFTIME'
    write(*,*) '3 - Read WSE Dataset and output a Max WSE dataset for all times'
    write(*,*) '4 - Take multiple solution datasets and merge into one file'
    write(*,*) '5 - Extract data for one cell from mult. datasets to text file'
    write(*,*) '9 - Exit'
    write(*,*) '*****************************************************************'
    write(*,*) ''
100 read (*,*) ichoice
    
    select case (ichoice)
    case (1) 
      call CMS_date2reftime
      call clear_screen       !This will present a fresh menu after it returns
    case (2) 
      call CMS_reftime2date
      call clear_screen
    case (3)
      call CMS_MaxWSE
      call clear_screen
    case (4)                  !Added 02/09/2022  MEB
      call MultiDatasets_toOne
      call clear_screen
    case (5)                  !Added 03/16/2022  MEB
      call ExtractMultiDatasets_toText
      call clear_screen
    case (9)
      STOP
    case default
      call diag_print_message ('Selection not available, make another selection')
      goto 100
    end select
#else    
    write(*,*) '*****************************************************************'
    write(*,*) 'Make a choice from the following CMS Tools (limited on Linux)'
    write(*,*) '3 - Read WSE Dataset and output a Max WSE dataset for all times'
    write(*,*) '9 - Exit'
    write(*,*) '*****************************************************************'
    write(*,*) ''
200 read (*,*) ichoice
    
    select case (ichoice)    
    case (3)
      call CMS_MaxWSE
      call clear_screen
    case (9)
      STOP
    case default
      call diag_print_message ('Selection not available, make another selection')
      goto 200
    end select
#endif

    goto 50
    
    contains
    
!******************************************************
      subroutine clear_screen
      
#ifdef _WIN32
        call system('cls')
#else
        call system('clear')
#endif

      return
      end subroutine clear_screen

!******************************************************
      subroutine get_correct_fraction_dset_name (in_str,out_str)
        implicit none
        character(len=*), intent(in)  :: in_str
        character(len=*), intent(out) :: out_str
        character(100) :: first, middle, last
        integer :: iloc1,iloc2
        
        iloc1 = index(in_str,'_')
        first = in_str(1:iloc1-1)
        middle = in_str(iloc1:iloc1+2)
        
        iloc1 = index(in_str,'(')
        iloc2 = index(in_str,')')
        last = in_str(iloc1:iloc2)
        
        out_str = trim(first)//trim(middle)//' '//trim(last)
        
      return
      end subroutine get_correct_fraction_dset_name

      
!******************************************************
      subroutine find_in_output_strings (a_str, member)
	    use out_def, only: noutputs, aoutputs
	    implicit none
		character(len=*), intent(in) :: a_str
		integer, intent(out)         :: member
		
		integer  :: i,iloc
		logical  :: found = .false.
		
		do i=1, noutputs
		  iloc = index(a_str,trim(aoutputs(i)))
		  if (iloc .gt. 0) then
			found = .true.
			exit
		  endif
		enddo
		if (found) then
		  member = i
		else
		  member = -1
		endif
		
		return
      end subroutine find_in_output_strings
		  
#ifdef _WIN32		
!******************************************************
! Read input from the keyboard to choose number of files/datasets to extract data for one cell ID into a text file
! written by Mitchell Brown, USACE-CHL  3/16/2022
!****************************************************** 
      subroutine ExtractMultiDatasets_toText
        use diag_lib,  only: diag_print_message
        use geo_def,   only: mapid,idmap
        use out_lib,   only: OPEN_CREATE_DATASET   
		use out_def,   only: noutputs, aoutputs, aoutput_lengths, simlabel
        use const_def, only: rad2deg
        use xmdf
        
        implicit none        
        integer            :: cellID, a_Number, nDsets, sMaxPathLength, vMaxPathLength, aLen
        integer            :: numFiles, i, j, k, ierr, fid, gid, nfid, ndid, ncells, a_NumTimes
        integer            :: nSDsets, nVDsets, complete, iloc, member, nDim, labelLength
        real(4)            :: num, wDir, wAbs, uVal, vVal
        logical            :: exists
        character(len=100) :: newFile, a_Path, cardfile, aString
        character(len=500) :: sPath(1),vPath(1), cardname
        character(len=5)   :: which
        character          :: achoice
        
        character(len=100), allocatable :: filenames(:), dsets(:), paths(:)
        integer, allocatable :: Dims(:), lengths(:)
        real(4),allocatable  :: SValues(:), VValues(:,:)
        real(8),allocatable  :: Times(:)
        
       !Determine output lengths
		do k=1,noutputs
          aoutput_lengths(k) = len_trim(aoutputs(k))
        enddo
      
        !Read card file to get geometry information to map ID correctly.
10      write(*,*)' '
        write(*,'(A)',advance='no') 'Enter name of parameter file (*.cmcards): '
        read(*,*) cardfile
        inquire(FILE=trim(cardfile),EXIST=exists)
        if (.not. exists) then
          call diag_print_message('File not found, try again',' ')
          goto 10
        endif
        open(77,file=cardfile,status='unknown')
        do
          read(77,*) cardname
          if (trim(cardname)=='GRID_FILE') exit
        enddo
        call geo_cards(cardname, exists, .false.)
        do 
          read(77,*) cardname
          if (trim(cardname)=='BATHYMETRY_DATASET') exit
        enddo
        call geo_cards(cardname, exists, .false.)
        do 
          read(77,*) cardname
          if (trim(cardname)=='GRID_CELL_TYPES') exit
        enddo
        call geo_cards(cardname, exists, .false.)
        call read_grid_xmdf
        close(77)
        
        write(*,*)' '
        write(*,'(A)',advance='no') 'Enter Cell ID for which to extract data: '
        read(*,*) cellID
        !map SMS ID to internal ID
        cellID = mapid(cellID)
        
        write(*,*)' '
        write(*,'(A)',advance='no') 'How many datasets to extract information for? '
        read(*,*) numFiles
        allocate(filenames(numFiles),dsets(numFiles),Dims(numFiles))
        
        !Loop through individual files, reading all datasets included in each
		i = 1
        do
          which = ''; aLen = 0
          if (numFiles .gt. 1) then
            which = 'first '
            if (i .ne. 1) which = 'next '
            aLen = len_trim(which)+1
          endif
          write(*,*)' '
          write(*,'(A)',advance='no') 'Enter the '//adjustl(which(1:aLen))//'filename which holds a dataset to merge: '
          read(*,*) filenames(i)
          inquire(FILE=trim(filenames(i)),EXIST=exists)
          if (.not.exists) then
            write(*,'(A)') 'File not found, try again'
            cycle
		  endif
		  
          call XF_OPEN_FILE(trim(filenames(i)),READONLY,fid,ierr)
          call XF_GET_SCALAR_DATASETS_INFO (fid, nSDsets, sMaxPathLength, ierr)
          call XF_GET_VECTOR_DATASETS_INFO (fid, nVDsets, vMaxPathLength, ierr)
          nDsets = nSDsets + nVDsets
          allocate(Paths(nDsets),lengths(nDsets))
          Paths=' '; sPath = ' '; vPath = ' ';
          
          !If no available datasets, try looking under the 'Datasets' folder (this file probably written by SMS, not CMS)
          gid = -1
          if (nDsets .eq. 0) then 
            call XF_OPEN_GROUP(fid, 'Datasets', gid, ierr)
            call XF_GET_SCALAR_DATASETS_INFO (gid, nSDsets, sMaxPathLength, ierr)
            call XF_GET_VECTOR_DATASETS_INFO (gid, nVDsets, vMaxPathLength, ierr)
            nDsets = nSDsets + nVDsets
            if (nDsets .eq. 0) then
              call diag_print_message('Error: no datasets found within the selected file. Preset try again')
              write(*,*) 'Press the Enter key to continue'
              read(*,*) 
              return
            endif
          endif          

          !Get paths
          if (gid .eq. -1) gid = fid       !If still equal to -1, then the first attempt worked
          call XF_GET_SCALAR_DATASET_PATHS (gid, nSDsets, sMaxPathLength, sPath, ierr)
          call XF_GET_VECTOR_DATASET_PATHS (gid, nVDsets, vMaxPathLength, vPath, ierr)
		  
		  !Find simlabel in path to help determine which dataset and what the proper length is
		  iloc = index(sPath(1),'/',.false.)
		  simlabel = sPath(1)(1:iloc)
          labelLength = len_trim(simlabel)
		  
          !List scalar datasets in file
          write(*,*) ' '
          write(*,*) ' Datasets contained in '//trim(filenames(i))//': '
		  k = 1
          do j=1,nSDsets
            aString=''
            if (sPath(1)(1:8) .eq. 'Fraction') then                                        !MEB 08/11/22
              call get_correct_fraction_dset_name(sPath(1)(k:k+sMaxPathLength),aString)
              Paths(j) = 'Datasets/'//trim(aString)
            else
	          call find_in_output_strings (sPath(1)(k:k+sMaxPathLength), member)         !Pass in the string to determine which of the output dataset names it is.
		      lengths(j) = aoutput_lengths(member) + len_trim(simlabel) - 1          !Set the proper length (there is a lot of garbage characters in the array of paths)
              Paths(j) = sPath(1)(k:k+lengths(j))
            endif
            write(*,'(2x,i2,A)') j,' - '//trim(Paths(j))
			      k = k + sMaxPathLength                                                   !Have to parse through multiple paths in one long string to find each.
          enddo
          k = 1
          do j=1,nVDsets
			call find_in_output_strings (vPath(1)(k:k+vMaxPathLength), member)
			lengths(j+nSDsets) = aoutput_lengths(member) + labelLength - 1
            Paths(j+nSDsets) = vPath(1)(k:k+lengths(j+nSDsets))
            write(*,'(2x,i2,A)') j+nSDsets,' - '//trim(Paths(j+nSDsets))
			k = k + vMaxPathLength
          enddo
          write(*,*) ' '
          write(*,'(A)',advance='no') 'Enter number of the dataset to merge: '
          read(*,*) a_Number
          
		  !Add proper path name to the list of datasets to merge.  Also set the number of dimensions based on the type.
          dsets(i) = Paths(a_Number)
          if (a_Number .le. nSDsets) then 
            Dims(i) = 1
          else
            write(*,'(A)',advance='no') '  Do you want to compute direction from U/V components? Y/n: '
            read(*,*) achoice
            if (achoice == 'N' .or. achoice == 'n') then
              Dims(i) = 2
            else
              Dims(i) = 3
            endif
          endif
          deallocate(Paths,lengths)

          i = i + 1                        !Move on to next file
          if (i .gt. numFiles) then
            exit
          endif
        enddo
        
        !Now open datasets and copy to new file
        write(*,*)' '
        write(*,'(A)',advance='no') 'Enter name of text file to export information to: '
        read(*,*) newFile
        inquire(FILE=trim(newFile),EXIST=exists)
        if (exists) then               !If it exists, delete and write over it
          open (11,file=newFile)
          close (11,status='delete')
        endif
        open(11, file=newfile)
        write(11,'(i0)') numFiles
        close(11)
        
        !Iterate through the files to extract the data from XMDF files and then write information to the text file.
        do i=1,numFiles
          call XF_OPEN_FILE (trim(filenames(i)),READONLY,fid,ierr)
          call XF_OPEN_GROUP(fid,trim(dsets(i)),gid,ierr)
          
          call XF_GET_DATASET_NUM_TIMES(gid, a_NumTimes, ierr)
          if(allocated(Times)) deallocate(Times)
          allocate(Times(a_NumTimes))
          Times = 0.0
          
          open(11, file=newfile, ACCESS='APPEND')
          select case (Dims(i))
          case (1)
            write(11,'(2(i0,2x),A)') i, a_NumTimes, trim(dsets(i))
          case (2)
            write(11,'(2(i0,2x),"U and V",3x,A)') i, a_NumTimes, trim(dsets(i))
          case default
            write(11,'(2(i0,2x),"U, V, and Dir",3x,A)') i, a_NumTimes, trim(dsets(i))
          end select
          
          call XF_GET_DATASET_TIMES(gid, a_NumTimes, Times, ierr)
          call XF_GET_DATASET_NUMVALS(gid, nCells, ierr)
          if(.not.allocated(Svalues)) allocate (Svalues(nCells),Vvalues(2,nCells))
          Svalues = 0.0; Vvalues = 0.0
          nDim = Dims(i)

          complete = 0
          write(*,*) ' '
          write(*,'(A,i0,A)') 'Extracting ',a_NumTimes, ' times from dataset: '//trim(dsets(i))
          do j=1,a_NumTimes
            select case (nDim)
            case (1)
              call XF_READ_SCALAR_VALUES_TIMESTEP(gid, j, nCells, Svalues, ierr)
              write(11,'(F13.3,2x,F13.8)') Times(j), Svalues(cellID)
            case (2) 
              call XF_READ_VECTOR_VALUES_TIMESTEP(gid, j, nCells, 2, VValues, ierr)
              write(11,'(F13.3,2(2x,F13.8))') Times(j), Vvalues(1,cellID), Vvalues(2,cellID)
            case default
              call XF_READ_VECTOR_VALUES_TIMESTEP(gid, j, nCells, 2, VValues, ierr)
              uVal = Vvalues(1,cellID) ; vVal = Vvalues(2,cellID)
              wAbs = uVal*uVal + vVal*vVal
              if (wAbs .ne. 0.0) then
                wDir = atan2(uVal/wAbs, vVal/wAbs)
                wDir = wDir * rad2deg
              else 
                wDir = 0
              endif
              write(11,'(F13.3,3(2x,F13.8))') Times(j), Vvalues(1,cellID), Vvalues(2,cellID), wDir
            end select
            
            num = real(j)/real(a_NumTimes)*100
            if (num .ge. complete) then
              write(*,'(2x,I3,x,A)') complete, 'Percent Complete'
              complete = complete + 5
            endif
          enddo
          call XF_CLOSE_GROUP(gid, ierr)  !Close old group
          call XF_CLOSE_FILE (fid, ierr)  !Close old file
          close(11)
        enddo

        write(*,*) 'Saved file: '//trim(newFile)
        write(*,*) ' '
        write(*,*) 'Press the Enter key to continue'
        read(*,*) 
        

      return
      end subroutine ExtractMultiDatasets_toText
      
!******************************************************
! Read input from the keyboard to choose number of files/datasets to merge from individual files into a common file
! - needed by Honghai for combining newer solutions into one file for PTM.
! written by Mitchell Brown, USACE-CHL  2/09/2022
!****************************************************** 
	  subroutine MultiDatasets_toOne
        use diag_lib, only: diag_print_message
        use out_lib,  only: OPEN_CREATE_DATASET   
		use out_def,  only: noutputs, aoutputs, aoutput_lengths, simlabel
        use comvarbl, only: reftime
        use xmdf
        
        implicit none
        integer   :: numFiles, i, j, k, ierr, fid, gid, nfid, ndid, ncells, a_NumTimes
        integer   :: a_Number, nDsets, sMaxPathLength, vMaxPathLength, nDim
        integer   :: nSDsets, nVDsets, complete, iloc, member
        real(4)   :: Num
        logical   :: exists

        integer,allocatable :: Dims(:), lengths(:)
        real(4),allocatable :: SValues(:), VValues(:,:)
        real(8),allocatable :: Times(:)
        character(len=100),allocatable :: filenames(:), dsets(:)
        character(len=100),allocatable :: paths(:)
        character(len=100) :: a_Path, newFile, units, which, aString
        character(len=500) :: sPath(1),vPath(1)
        character(len=5)   :: aExt

        !Determine output lengths
		do k=1,noutputs
          aoutput_lengths(k) = len_trim(aoutputs(k))
        enddo

        write(*,*)' '
        write(*,'(A)',advance='no') 'How many solution datasets to merge into one file? '
        read(*,*) numFiles
        allocate(filenames(numFiles),dsets(numFiles),Dims(numFiles))
        
        !Loop through individual files, reading all datasets included in each
		i = 1
        do
          which = 'first'
          if (i .ne. 1) which = 'next'
          write(*,*)' '
          write(*,'(A)',advance='no') 'Enter '//trim(which)//' filename which holds a dataset to merge: '
          read(*,*) filenames(i)
          inquire(FILE=trim(filenames(i)),EXIST=exists)
          if (.not.exists) then
            write(*,'(A)') 'File not found, try again'
            cycle
		  endif
		  
          call XF_OPEN_FILE(trim(filenames(i)),READONLY,fid,ierr)
          call XF_GET_SCALAR_DATASETS_INFO (fid, nSDsets, sMaxPathLength, ierr)
          call XF_GET_VECTOR_DATASETS_INFO (fid, nVDsets, vMaxPathLength, ierr)
          nDsets = nSDsets + nVDsets
          
          !If no available datasets, try looking under the 'Datasets' folder (this file probably written by SMS, not CMS)  MEB 08/11/22
          gid = -1
          if (nDsets .eq. 0) then 
            call XF_OPEN_GROUP(fid, 'Datasets', gid, ierr)
            call XF_GET_SCALAR_DATASETS_INFO (gid, nSDsets, sMaxPathLength, ierr)
            call XF_GET_VECTOR_DATASETS_INFO (gid, nVDsets, vMaxPathLength, ierr)
            nDsets = nSDsets + nVDsets
            if (nDsets .eq. 0) then
              call diag_print_message('Error: no datasets found within the selected file. Preset try again')
              write(*,*) 'Press the Enter key to continue'
              read(*,*) 
              return
            endif
          endif
          
          if (gid .eq. -1) gid = fid       !If still equal to -1, then the first attempt worked  MEB 08/11/22
          allocate(Paths(nDsets),lengths(nDsets))
          Paths=' '; sPath = ' '; vPath = ' ';
          
          !Get paths
          call XF_GET_SCALAR_DATASET_PATHS (gid, nSDsets, sMaxPathLength, sPath, ierr)
          call XF_GET_VECTOR_DATASET_PATHS (gid, nVDsets, vMaxPathLength, vPath, ierr)
          call XF_CLOSE_FILE(fid,ierr)
		  
		  !Find simlabel in path to help determine which dataset and what the proper length is
		  iloc = index(sPath(1),'/',.false.)
		  simlabel = sPath(1)(1:iloc)
		  
          !List scalar datasets in file
          write(*,*) ' '
          write(*,*) ' Datasets contained in '//trim(filenames(i))//': '
		  k = 1
          do j=1,nSDsets
            aString=''
            if (sPath(1)(1:8) .eq. 'Fraction') then                                        !MEB 08/11/22
              call get_correct_fraction_dset_name(sPath(1)(k:k+sMaxPathLength),aString)
              Paths(j) = 'Datasets/'//trim(aString)
            else
			  call find_in_output_strings (sPath(1)(k:k+sMaxPathLength), member)           !Pass in the string to determine which of the output dataset names it is.
			  lengths(j) = aoutput_lengths(member) + len_trim(simlabel) - 1                !Set the proper length (there is a lot of garbage characters in the array of paths)
              Paths(j) = sPath(1)(k:k+lengths(j))
            endif
            write(*,'(2x,i2,A)') j,' - '//trim(Paths(j))
			k = k + sMaxPathLength                                                         !Have to parse through multiple paths in one long string to find each.
          enddo
		  k = 1
          do j=1,nVDsets
			call find_in_output_strings (vPath(1)(k:k+vMaxPathLength), member)
			lengths(j+nSDsets) = aoutput_lengths(member) + len_trim(simlabel) - 1
            Paths(j+nSDsets) = vPath(1)(k:k+lengths(j+nSDsets))
            write(*,'(2x,i2,A)') j+nSDsets,' - '//trim(Paths(j+nSDsets))
			k = k + vMaxPathLength
          enddo
          write(*,*) ' '
          write(*,'(A)',advance='no') 'Enter number of the dataset to merge: '
          read(*,*) a_Number
          
		  !Add proper path name to the list of datasets to merge.  Also set the number of dimensions based on the type.
          dsets(i) = Paths(a_Number)
          if (a_Number .le. nSDsets) then 
            Dims(i)= 1
          else
            Dims(i)= 2
          endif
          
          deallocate(Paths,lengths)

          i = i + 1                        !Move on to next file
          if (i .gt. numFiles) then
            exit
          endif
        enddo
        
        !Now open datasets and copy to new file
        write(*,*)' '
        write(*,'(A)',advance='no') 'Enter name of new file to create with merged datasets: '
        read(*,*) newFile
        call fileext (newFile,aext)
        if (aext(1:1) .eq. ' ') newFile = trim(newFile)//'.h5'
        
        inquire(FILE=trim(newFile),EXIST=exists)
        if (exists) then
          open (11,file=newFile)
          close (11,status='delete')
        endif
        call XF_CREATE_FILE(trim(newFile),readwrite,nfid,ierr)
        call XF_CLOSE_FILE(nfid,ierr)
        
        do i=1,numFiles
          call XF_OPEN_FILE (trim(newFile),readwrite,nfid,ierr)
          call XF_OPEN_FILE (trim(filenames(i)),READONLY,fid,ierr)
          call XF_OPEN_GROUP(fid,trim(dsets(i)),gid,ierr)
          
          call XF_READ_REFTIME(gid,Reftime,ierr)
          call XF_GET_DATASET_NUM_TIMES(gid, a_NumTimes, ierr)
          if(.not.allocated(Times)) allocate(Times(a_NumTimes))
          Times = 0.0
          
          call XF_GET_DATASET_TIMES(gid, a_NumTimes, Times, ierr)
          call XF_GET_DATASET_NUMVALS(gid, nCells, ierr)
          if(.not.allocated(Svalues)) allocate (Svalues(nCells),Vvalues(2,nCells))
          Svalues = 0.0; Vvalues = 0.0
          nDim = Dims(i)
          call OPEN_CREATE_DATASET(nfid,trim(dsets(i)),ndid,nDim,'',ierr)

          complete = 0
          write(*,*) ' '
          write(*,'(A)') ' Copying dataset: '//trim(dsets(i))
          do j=1,a_NumTimes
            if(nDim == 1) then
              call XF_READ_SCALAR_VALUES_TIMESTEP(gid, j, nCells, Svalues, ierr)
              call XF_WRITE_SCALAR_TIMESTEP(ndid,Times(j),nCells, Svalues, ierr)
            else
              call XF_READ_VECTOR_VALUES_TIMESTEP(gid, j, nCells, nDim, VValues, ierr)
              call XF_WRITE_VECTOR_TIMESTEP(ndid,Times(j),nCells, nDim, VValues, ierr)
            endif
            num = real(j)/real(a_NumTimes)*100
            if (num .ge. complete) then
              write(*,'(2x,I3,x,A)') complete, 'Percent Complete'
              complete = complete + 10
            endif
          enddo
          call XF_CLOSE_GROUP(ndid, ierr) !Close new group
          call XF_CLOSE_FILE (nfid, ierr) !Close new file
          call XF_CLOSE_GROUP(gid, ierr)  !Close old group
          call XF_CLOSE_FILE (fid, ierr)  !Close old file
        enddo

        write(*,*) 'Saved file: '//trim(newFile)
        write(*,*) ' '
        write(*,*) 'Press the Enter key to continue'
        read(*,*) 
        
      return
      end subroutine MultiDatasets_toOne
#endif            
      
!******************************************************
      subroutine CMS_MaxWSE
        use diag_lib, only: diag_print_message
        
        implicit none
        character :: filename*100, aext*10
        logical   :: exists

        do
          write(*,'(A)') 'Enter filename of WSE solution file (*_wse.h5, *_eta.dat, etc)'
          read(*,*) filename
          inquire(FILE=trim(filename),EXIST=exists)
          if (.not.exists) then
            write(*,'(A)') 'File not found, try again'
          else
            exit
          endif
        enddo

        call fileext (filename,aext)
        select case (aext)
          case('h5')
#ifdef _WIN32            
            call CMS_MaxWSE_h5(trim(filename))
#else
            call diag_print_message("XMDF/H5 file functionality is missing for Linux")
#endif
          case('dat')
            call CMS_MaxWSE_dat(trim(filename))
          case default
            call diag_print_message("This tool presently only works on '.dat' and '.h5' files")
        end select
        write(*,*) ' '
        write(*,*) 'Press the Enter key to continue'
        read(*,*) 
          
        return
      end subroutine CMS_MaxWSE

!*************************************************************
      subroutine CMS_MaxWSE_dat (filename)
        use in_def,   only: scaldattype,vecdattype 
        use in_lib,   only: read_dat
        use out_lib,  only: write_scal_dat_file
        use comvarbl, only: reftime
        use size_def, only: ncells
        
        implicit none
        character(len=*), intent(in) :: filename
        real(4), allocatable :: max_array(:)
        
        integer :: nscal, nvec, nTimes, wseRec = -1
        type(scaldattype), pointer :: scaldat(:)
        type(vecdattype),  pointer :: vecdat(:)
        integer :: i, j, iloc
        logical :: found = .false.
        character :: newfile*100, aName*100
        
        call read_dat (filename,nscal,scaldat,nvec,vecdat)
        do i=1,nscal
          if(scaldat(i)%NAME == 'Water_Elevation') then
            nTimes  = scaldat(i)%NT
            nCells  = scaldat(i)%NC
            reftime = scaldat(i)%reftime
            allocate(max_array(nCells))
            max_array = -1000.0
            found = .true.
            wseRec = i
            exit
          else
            cycle
          endif
        enddo
        
        if (.not.found) then
          call diag_print_message ("Error: 'Water_Elevation' dataset not found in file")
          write(*,*) 'Press the Enter key to continue'
          read(*,*) 
          return
        endif
        write(*,*)' '
        write(*,'(A)') ' Parsing existing Water Elevations'
        do i=1,nTimes
          do j=1,nCells
            max_array(j) = max(scaldat(wseRec)%VAL(j,i),max_array(j))
          enddo
        enddo
        
        !Create new file for the Max Water Elevations based on original filename
        iloc    = index(filename,'_eta')
        aName   = filename(1:iloc-1)
        newfile = trim(aName)//'_maxeta.dat'
        write(*,'(A)')" Writing Maximum elevations to file: '"//trim(newfile)//"'"

        call write_scal_dat_file(aName,'Maximum_WSE','maxeta',max_array)

        return
      end subroutine CMS_MaxWSE_dat
      
#ifdef _WIN32      
!*************************************************************
      subroutine CMS_MaxWSE_h5 (filename)
        use XMDF
        use out_lib, only: OPEN_CREATE_DATASET
        use comvarbl, only: reftime
        use diag_lib, only: diag_print_error
        
        implicit none
        character(len=*), intent(in) :: filename
        
        integer :: fid, gid, ierr, a_Number, a_MaxPathLength, iloc
        integer :: nfid, ndid
        integer :: i, j, a_NumTimes, nCells, complete
        real(4) :: Num
        real(4), allocatable :: Values(:), max_array(:)
        real(8), allocatable :: Times(:)
        character(100) :: newfile, a_Path
        character(100), allocatable :: Paths(:)
        logical :: exists

        call XF_OPEN_FILE(trim(filename),readwrite,fid,ierr) !Open XMDF file
        if(ierr<0)then
          write(*,'(A)') 'Error opening file'
        endif          
        call XF_GET_SCALAR_DATASETS_INFO (fid, a_Number, a_MaxPathLength, ierr)
        allocate(Paths(a_Number))
        do i=1,a_Number
          call XF_GET_SCALAR_DATASET_PATHS (fid, i, a_MaxPathLength, Paths, ierr)
          a_Path = Paths(i)
          iloc = index(a_Path,'Water_Elevation')
          if (iloc <= 0) then
            cycle
          else
            a_Path(iloc+15:100) = ' '
            exit
          endif
        enddo
        deallocate(Paths)
        
        if (i > a_Number) then
          write (*,'(A)') 'No Water Elevation Dataset found'
          stop
        endif 
        
        call OPEN_CREATE_DATASET(fid,trim(a_Path),gid,1,'m',ierr)
        
        !Parse through the datasets in the file and get max value
        call XF_READ_REFTIME(gid, Reftime, ierr)
        call XF_GET_DATASET_NUM_TIMES(gid, a_NumTimes, ierr)
        allocate(Times(a_NumTimes))
        times = 0.0 ;
        
        call XF_GET_DATASET_TIMES(gid, a_NumTimes, Times, ierr)
        call XF_GET_DATASET_NUMVALS(gid, nCells, ierr)
        allocate(max_array(nCells),Values(nCells))
        max_array = -1000.0 ; values = 0.0
        complete = 0
        write(*,*) ' '
        write(*,'(A)') ' Parsing existing Water Elevations'
        do i=1,a_NumTimes
          call XF_READ_SCALAR_VALUES_TIMESTEP(gid, i, nCells, Values, ierr)
          do j=1,nCells
            max_array(j) = max(Values(j),max_array(j))
          enddo
          num = real(i)/real(a_NumTimes)*100
          if (num .ge. complete) then           
            write(*,'(2x,I3,x,A)') complete, 'Percent Complete'
            complete = complete + 10
          endif
        enddo
        call XF_CLOSE_GROUP(gid,ierr)  !Close old dataset    
        call XF_CLOSE_FILE(fid,ierr)   !Close old file
        
        !Create new file for the Max Water Elevations based on original filename
        iloc = index(filename,'_wse')
        newfile = filename(1:iloc-1)//'_maxwse.h5'
        iloc = index(a_Path,'Water_Elevation')
        a_Path = a_Path(1:iloc-1)//'Max_Water_Elevation'

        inquire(FILE=trim(newfile), EXIST=exists)
        if (exists) then
          open(100,file=newfile,iostat=ierr)
          close(100,status='delete',iostat=ierr)
        endif            
        call XF_CREATE_FILE(trim(newfile),readwrite,nfid,ierr)
        call OPEN_CREATE_DATASET(nfid,trim(a_Path),ndid,1,'m',ierr)
        if (ierr.ge.0) then
          write(*,*) ' '
          write(*,'(A,A)')" Writing Maximum elevations to file: '"//trim(newfile)//"'"
        else
          call diag_print_error('Error creating dataset')
        endif
        call XF_WRITE_SCALAR_TIMESTEP(ndid,0.d0,nCells,max_array,ierr)
        
        call XF_CLOSE_GROUP(ndid,ierr) !Close new dataset
        call XF_CLOSE_FILE(nfid,ierr)  !Close new file 
        
        return
      end subroutine CMS_MaxWSE_h5
      
!**************************************************************************
      subroutine CMS_date2reftime
        use XMDF
        use comvarbl, only: iyr, imo, iday, ihr, imin, isec, reftime
        implicit none
        integer :: ierr
        
        
        write(*,'(A)') 'Enter date (yyyy-mm-dd hh:mm:ss)'
        call card_datetime(5,iyr,imo,iday,ihr,imin,isec) 
        call XF_CALENDAR_TO_JULIAN(0,iyr,imo,iday,ihr,imin,isec,reftime,ierr)
        write(*,*)''
        write(*,'(A,F0.6)') 'RT_JULIAN date is: ',reftime
        write(*,'(A)') 'Press <enter> key to continue.'
        read(*,*)
        write(*,*)''
        
        return
      end subroutine CMS_date2reftime
    
      subroutine CMS_reftime2date
        use XMDF
        use comvarbl, only: iyr, imo, iday, ihr, imin, isec, reftime
        implicit none
        integer :: ierr, era = 0

100     format('Calendar date is: ',i4.4,'-',i2.2,'-',i2.2, 1x ,i2.2,':',i2.2,':',i2.2)
        
        write(*,*) 'Enter reference time:'
        read (*,*) reftime 
        call XF_JULIAN_TO_CALENDAR(era,iyr,imo,iday,ihr,imin,isec,reftime,ierr)
        write(*,*)''
        write(*,100) iyr, imo, iday, ihr, imin, isec
        write(*,'(A)') 'Press <enter> key to continue.'
        read(*,*)
        write(*,*)''
        
        return
      end subroutine CMS_reftime2date
#endif
      
    end subroutine CMS_Tools_Dialog
!******************************************************    
    
!******************************************************
    logical function findCard(aFile,aCard,aValue)    
! Find a given card in the cardfile and return the rest of the selected line in the variable, aValue.
! written by Mitch Brown, USACE-CHL 06/25/2021
!
! The correct interface is given below
        !interface  
        ! function findCard(aFile,aCard,aValue)
        !    character(len=*),intent(in)    :: aFile
        !    character(len=*),intent(in)    :: aCard
        !    character(len=100),intent(out) :: aValue
        !    logical :: findCard
        !  end function
        !end interface
!******************************************************    
    implicit none
    character(len=*), intent(in)    :: aFile
    character(len=*), intent(in)    :: aCard
    character(len=100), intent(out) :: aValue
    
    character(len=80)  :: testCard
    character(len=100) :: aLine
    logical            :: foundCard = .false.
    integer            :: iloc
    
    open(77,file=aFile,status='unknown')
    do 
      read(77,*,end=100) testCard
      if (trim(testCard) == trim(aCard)) then
        backspace(77)
        read(77,'(A100)') aLine
        iloc = index(aLine,' ')
        aLine = adjustL(aLine(iloc:))
        read(aLine, '(A100)') aValue
        findCard = .true.
        return
        exit
      else
        findCard = .false.
      endif
    enddo

100 findCard = .false.

    return
    end function findcard
!******************************************************    
   
    
!******************************************************
    subroutine fileparts(astr,apath,aname,aext)    
! Determines the parts of a file name
! written by Alex Sanchez, USACE-CHL
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: apath,aname,aext
    !Internal
    integer :: i,k,nn    
    
    nn = len_trim(astr)    
    !Determine path
    apath = ''   
    do k=nn,1,-1
      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
        apath = astr(1:k)     
        exit
      endif
    enddo
    
    !Determine name and extension
    aname = astr(max(k+1,1):nn)
    aext = ''    
    do i=nn,k+1,-1
      if(astr(i:i)=='.')then
        aext = astr(i+1:nn)    
        aname = trim(astr(max(k+1,1):i-1))
        return
      endif
    enddo          

    return
    end subroutine fileparts
    
!!******************************************************
!    subroutine fileparts(astr,apath,aname,aext)    
!! Determines the parts of a file name
!! written by Alex Sanchez, USACE-CHL
!!******************************************************    
!    implicit none
!    !Input/Output
!    character(len=*),intent(in) :: astr
!    character(len=*),intent(inout) :: apath
!    character(len=*),intent(inout),optional :: aname,aext
!    !Internal
!    integer :: i,k,nn    
!    
!    nn = len_trim(astr)    
!    !Determmin path
!    apath = ''   
!    do k=nn,1,-1
!      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
!        apath = astr(1:k)     
!        exit
!      endif
!    enddo    
!    
!    if(.not.present(aname)) return
!    
!    !Determine name and extension
!    aname = astr(max(k+1,1):nn)
!    if(present(aext)) aext = ''    
!    do i=nn,k+1,-1
!      if(astr(i:i)=='.')then
!        if(present(aext)) aext = astr(i+1:nn)    
!        aname = trim(astr(max(k+1,1):i-1))
!        return
!      endif
!    enddo          
!
!    return
!    end subroutine fileparts
    
!******************************************************
    subroutine filepath(astr,apath)  
! Determines the path of a file
! written by Alex Sanchez, USACE-CHL
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: apath
    !Internal
    integer :: k,nn    
    
    nn = len_trim(astr)
    apath = ''   
    do k=nn,1,-1
      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
        apath = astr(1:k)     
        return
      endif
    enddo       

    return
    end subroutine filepath
    
!******************************************************
    subroutine fileext(astr,aext)  
! Determines the extension of a file
! written by Alex Sanchez, USACE-CHL    
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: aext
    !Internal
    integer :: i,k,nn
        
    !Determine extension
    nn = len_trim(astr)
    k = index(astr,'\')
    if(k==0) k = index(astr,'/')
    aext = ' '
    do i=nn,max(1,k+1),-1
      if(astr(i:i)=='.')then
        aext = astr(i+1:nn)        
        exit
      endif
    enddo    
    
    return
    end subroutine fileext
    
!******************************************************
    subroutine remove_underscores(astr)    
!******************************************************    
    implicit none
    character(len=*) :: astr
    integer :: i,nn
    
    nn = len_trim(astr)
    do i=1,nn
      if(astr(i:i)=='_')then
        astr(i:i) = ' '
      endif
    enddo    
    
    return
    end subroutine remove_underscores

!********************************************************
    subroutine uppercase(str)
!********************************************************
    implicit none
    character(len=*), intent(inout) :: str
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'a') .and. lle(str(i:i),'z'))then
        str(i:i) = achar(iachar(str(i:i)) - del)
      endif
    enddo

    return
    end subroutine uppercase
    
!********************************************************
    function toUpper(str) result(aString)
!********************************************************
    implicit none
    character(len=*), intent(in)  :: str
    character(len=len(str))       :: aString
    integer, parameter :: offset = iachar('a') - iachar('A')
    integer :: i

    aString = str
    do i=1, len(aString)
      select case(aString(i:i))
      case("a":"z")
        aString(i:i) = achar(iachar(aString(i:i)) - offset)
      end select
    enddo

    return
    end function toUpper
    
!********************************************************
    function toLower(str) result(aString)
!********************************************************
    implicit none
    character(len=*), intent(in)  :: str
    character(len=len(str))       :: aString
    integer, parameter :: offset = iachar('a') - iachar('A')
    integer :: i

    aString = str
    do i=1,len(aString)
      select case(aString(i:i))
      case("A":"Z")
        aString(i:i) = achar(iachar(aString(i:i)) + offset)
      end select
    enddo
    
    return
    end function toLower

!********************************************************
    subroutine lowercase(str)
!********************************************************
    implicit none
    character(len=*), intent(inout) :: str
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'A') .and. lle(str(i:i),'Z'))then
        str(i:i) = achar(iachar(str(i:i)) + del)
      endif
    enddo

    return
    end subroutine lowercase

!*************************************************************
    function trimspace(str)
! Return index of first space
!*************************************************************    
    integer :: nn,trimspace
    character(len=*),intent(in):: str
    
    nn = len(str)      
    do while(str(nn:nn)==' ' .and. nn>0)
      nn = nn - 1
    enddo
    trimspace = nn      
    
    return
    end function trimspace

!*****************************************************************
    subroutine removequotes(string)
!*****************************************************************
    integer :: i,j,js,je,nn
    character(len=*),intent(inout) :: string
    
    nn=len_trim(string)
    js=1; je=nn
    do i=1,nn
      if(string(i:i)=='"')then
        do j=i+1,nn
          if(string(j:j)=='"')then
            string = string(i+1:j-1)
            return
          endif
        enddo
      endif
    enddo

    return
    end subroutine removequotes
    
!*****************************************************************
    subroutine countquotes(string,count)
!*****************************************************************
    integer :: i,j,js,je,nn
    character(len=*),intent(in) :: string
    integer, intent(out) :: count
    
    count=0
    nn=len_trim(string)
    js=1; je=nn
    do i=1,nn
      if(string(i:i)=='"') then
        count=count+1  
      endif
    enddo

    return
    end subroutine countquotes

!!***********************************************************************    
!    subroutine copy_file_to_temp(file,tempfile)
!!***********************************************************************    
!    !use ifport, only: system
!    use DFLIB, only: systemqq
!    use diag_def
!    use diag_lib
!    implicit none
!    !Input
!    character(len=*),intent(in) :: file
!    !Output
!    character(len=*),intent(out) :: tempfile
!    !Internal
!    integer :: ierr
!    character(len=200) :: thepath
!    character(len=200) :: thename
!    character(len=10) :: theext
!    logical :: foundfile
!    
!    !Name of temporary file
!    call fileparts(file,thepath,thename,theext)
!    tempfile = trim(thepath)//trim(thename)//'_temporary.h5'
!    !call filepath(file,thepath)
!    !tempfile = trim(thepath)//'temp.h5'
!    
!    ! Check if temp file exists, if so, delete it
!    inquire(file=tempfile,exist=foundfile)
!    if(foundfile)then
!      open(100,file=tempfile,iostat=ierr)
!      close(100,status='delete',iostat=ierr)
!    endif
!
!    ! Make a copy of ICFILE to temp.h5 and read from that.
!    !res = system('copy '//trim(file)//' '//trim(tempfile)//' > trash.txt')
!    ierr = systemqq('copy '//trim(file)//' '//trim(tempfile)//' > trash.txt')
!    inquire(file='trash.txt',exist=foundfile)
!    if(foundfile)then  !remove scratch file
!      open(100,file='trash.txt',iostat=ierr)
!      close(100,status='delete',iostat=ierr)
!    endif
!    
!    return
!    end subroutine copy_file_to_temp
    
!************************************************************
    subroutine delete_file(file)
!************************************************************    
    implicit none
    integer :: ierr
    character(len=*),intent(in) :: file
    logical :: foundfile
    
    inquire(file=file,exist=foundfile)
    if(foundfile)then
      open(777,file=file,iostat=ierr)
      close(777,status='delete',iostat=ierr)
      !if(ierr/=0)then
        !res = systemqq('del '//trim(tfile))  
      !endif
    endif    
    
    return
    end subroutine delete_file
    
!-----------------------------------------------------------------------
!   converts a real*4 value with n significant digits to integers
!-----------------------------------------------------------------------
      subroutine split_real_to_integers (x, n, first, second) 

      integer*4      :: i, idot, n, lstring, first, second, loc
      real*4         :: x
      logical        :: significant_number
      character (len=20) :: string, part1, part2

      if (n < 2) n = 2
      if (n > 7) n = 7

      if (n == 2) then
          write (string, '(f0.2)') x
      else if (n == 3) then
          write (string, '(f0.3)') x
      else if (n == 4) then
          write (string, '(f0.4)') x
      else if (n == 5) then
          write (string, '(f0.5)') x
      else if (n == 6) then
          write (string, '(f0.6)') x
      else 
          write (string, '(f0.7)') x
      end if

      lstring = len_trim(string)
!
!     remove trailing zeros 
!
      idot = index (string(1:lstring),'.')
 
      do i = lstring, idot+2, -1
         if (string(i:i) == '0') then
             string(i:i) = ' '
             significant_number = .false.
         else
             significant_number = .true.
         end if
         if (significant_number) exit 
      end do

      loc=index(string,'.')
      part1 = string(:loc-1)
      part2 = string(loc+1:)

      read(part1,*) first
      read(part2,*) second
      
      return      
      end subroutine split_real_to_integers    
    
!************************************************************
    character(len=20) function Int2Str(k)
!   "Convert an integer to string."
!************************************************************
    implicit none
    integer, intent(in) :: k
    
    character(len=20) :: str
    
    write (str, '(i20)') k
    Int2Str = adjustl(str)
    
    end function Int2Str

!************************************************************
    subroutine check_TidalConstituent(name)
    !Drop '(' and change Lambda(2) to LDA2
!************************************************************
    character(len=*), intent(inout) :: name
    ! Locals
    integer k, iloc
    character(len=10) first,num,temp
    
    !Check names for parentheses (i.e., M(2) instead of M2) and if 'Lambda'
    iloc = 0
    iloc = index(name,'(')
    temp = ''
    if (iloc .ne. 0) then  !remove the ()
      if (name .eq. 'Lambda(2)') then
        name = 'Lda2'
      else
        temp=name
        first=temp(1:iloc-1)
        num=temp(iloc+1:iloc+1)
        name = trim(first)//trim(num)
      endif
    endif
!     name(k)=toUpper(trim(name(k)))
    
    return
    end subroutine
    
!************************************************************
Subroutine get_current_dir(dir_name)
!************************************************************
    USE ISO_C_BINDING

    Implicit NONE
    Character(LEN=:), ALLOCATABLE, Intent(INOUT) :: dir_name
    Type(C_PTR)             :: c_str_ptr
    Character(132), Pointer :: f_str_ptr
    Integer :: ilen, inull

    Interface
      Function c_getcwd(buf, size) BIND(C, NAME="getcwd")
        IMPORT :: C_PTR, C_SIZE_T
        Type(C_PTR),       VALUE :: buf
        Integer(C_SIZE_T), VALUE :: size
        Type(C_PTR)              :: c_getcwd
      End Function c_getcwd
    End Interface

    Interface
      Subroutine c_free(ptr) BIND(C,name="free")
        IMPORT :: C_PTR
        Implicit NONE
        Type(C_PTR), VALUE :: ptr
      End Subroutine c_free
    End Interface

    c_str_ptr = c_getcwd(C_NULL_PTR, 132_C_SIZE_T)

    Call C_F_POINTER(c_str_ptr, f_str_ptr)

    ilen  = LEN_TRIM(f_str_ptr)
    inull = INDEX(f_str_ptr, C_NULL_CHAR)
    If (inull /= 0) ilen = inull-1
    ilen = MAX(1, MIN(ilen,132))

    dir_name = TRIM(ADJUSTL(f_str_ptr(1:ilen)))

    NULLIFY(f_str_ptr)

    If (C_ASSOCIATED(c_str_ptr)) Then
      Call c_free(c_str_ptr)
    End If

  End Subroutine get_current_dir
