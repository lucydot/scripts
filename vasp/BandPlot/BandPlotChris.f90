
 Implicit NONE
 INTEGER i,j,k,l,m
 REAL,allocatable::Bands(:,:)
 INTEGER Nbands,Npts,NEigen,NIgnore
 Character(len=1)Junk
 LOGICAL :: file_exists

 INQUIRE(FILE="EIGENVAL", EXIST=file_exists)
 IF(.NOT. file_exists)GOTO 101

 OPEN(unit=1,file='EIGENVAL',Status='Unknown')

 READ(1,*)Junk
 READ(1,*)Junk
 READ(1,*)Junk
 READ(1,*)Junk
 READ(1,*)Junk
 READ(1,*)Junk,NEigen,Nbands
 ALLOCATE(Bands(NEigen,Nbands))
 DO i=1,NEigen
 READ(1,*)Junk
  DO j=1,NBands
   READ(1,*)l,Bands(i,j)
  ENDDO !j=1,NBands
 ENDDO !i=1,NEigen
 Close(1)
 WRITE(*,*)"Please enter the range over which you wish to plot the bands.(X Y)"
 READ(*,*) k,l
 WRITE(*,*)"Please enter the number of k-points to ignore"
 READ(*,*) NIgnore

 OPEN(unit=2,File="Bands.dat",Status='Unknown')
   DO j = k,l
    DO i=NIgnore+1,NEigen
    WRITE(2,*)i,(Bands(i,j))
   ENDDO
   WRITE(2,*) ""
 ENDDO
 GOTO 102

101 WRITE(*,*)"File EIGENVAL not found, aborting"
102 WRITE(*,*)"Program finished."

 END

