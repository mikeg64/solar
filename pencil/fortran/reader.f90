
PROGRAM reader


  IMPLICIT NONE

  ! file location settings
  character (len=*), parameter :: mag_field_dat = 'mag_field.dat'

  CHARACTER (LEN = 4) :: title = 'Mr.'
  CHARACTER (LEN = 15) :: name = 'Neil Smith'
  INTEGER :: age = 32
  CHARACTER (LEN = 35) :: address = '  2 Oak Lane, Manchester M23 4QD'
  INTEGER :: tel = 4455678
  integer, parameter :: bnx=64, bny=64 ! data in pencil shape
  real, dimension (:,:), allocatable :: Bz

  INTEGER :: rec_len		! record length to be inquired

  PRINT *, title, name, age, address, tel

  ! Assign the inquired record length to rec_len
  INQUIRE (IOLENGTH = rec_len) title, name, age, address, tel
  PRINT *, "The inquired record length = ", rec_len

  ! Open file for unformmated output using the inquired record length
  OPEN(UNIT = 1, FILE = 'temp.dat', RECL = rec_len, FORM = 'UNFORMATTED')
  WRITE(1) title, name, age, address, tel
  CLOSE(UNIT = 1)

  ! Open file for unformmated input
  OPEN(UNIT = 1, FILE = 'temp.dat', FORM = 'UNFORMATTED')
  READ(1) title, name, age, address, tel
  CLOSE(UNIT = 1)

  ! Print the list again (should be same as the first print)"
  PRINT *, title, name, age, address, tel
  inquire (iolength=rec_len) 1.0d0
  allocate (Bz(bnx,bny))
    ! Open file for unformmated input
  OPEN(UNIT = 1, FILE = mag_field_dat, FORM = 'UNFORMATTED', recl=rec_len*bnx*bny, access='direct')
  read (unit=1, rec=1) Bz
  print *, Bz
  CLOSE(UNIT = 1)

  Bz=1

  OPEN(UNIT = 1, FILE = 'new_mag_field.dat', FORM = 'UNFORMATTED', recl=rec_len*bnx*bny, access='direct')
  write (unit=1, rec=1) Bz
  print *, Bz
  CLOSE(UNIT = 1)

END PROGRAM
