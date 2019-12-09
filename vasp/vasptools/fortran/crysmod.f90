module crysmod
    ! This module provides utilities in order to manipulate a crystal object.
    ! Crystal object can be created from POSCAR/CONTCAR VASP files.
    ! In all functions and subroutines, angle are given or returned in degrees. If
    ! needed they are converted in radian localy.
    !

    implicit none

    type Atom
        ! atom type object
        character(len=3)           :: atname ! atom name
        real(kind=8), dimension(3) :: x      ! cartesian coordinates
        real(kind=8), dimension(3) :: r      ! reduce coordinates
    end type atom

    type Crystal
        ! crystal type object
        real(kind=8), dimension(3) :: veca, vecb, vecc
        ! lattice parameters :
        real(kind=8) :: a = -1.d0
        real(kind=8) :: b = -1.d0
        real(kind=8) :: c = -1.d0
        real(kind=8) :: alpha = -1.d0
        real(kind=8) :: beta  = -1.d0
        real(kind=8) :: gam   = -1.d0
        ! crystal name
        character(len=80) :: cname = ""
        ! atoms type
        integer                                     :: nType
        integer, dimension(:), allocatable          :: nAtomPerType
        character(len=3), dimension(:), allocatable :: atomType
        ! crystal atoms
        integer                               :: Natoms
        type(Atom), dimension(:), allocatable :: atoms
        ! transformation matrixes
        real(kind=8), dimension(3,3) :: Mcart2red, Mred2cart
    end type Crystal

    real(kind=8), parameter, private :: pi = 3.14159265359d0

    real(kind=8), dimension(81), parameter, private :: vtrans = (/ 0.d0, 1.d0, 1.d0, &
                                                                   0.d0, 1.d0, 0.d0, &
                                                                   0.d0, 1.d0,-1.d0, &
                                                                   0.d0, 0.d0, 1.d0, &
                                                                   0.d0, 0.d0, 0.d0, &
                                                                   0.d0, 0.d0,-1.d0, &
                                                                   0.d0,-1.d0, 1.d0, &
                                                                   0.d0,-1.d0, 0.d0, &
                                                                   0.d0,-1.d0,-1.d0, &
                                                                   1.d0, 1.d0, 1.d0, &
                                                                   1.d0, 1.d0, 0.d0, &
                                                                   1.d0, 1.d0,-1.d0, &
                                                                   1.d0, 0.d0, 1.d0, &
                                                                   1.d0, 0.d0, 0.d0, &
                                                                   1.d0, 0.d0,-1.d0, &
                                                                   1.d0,-1.d0, 1.d0, &
                                                                   1.d0,-1.d0, 0.d0, &
                                                                   1.d0,-1.d0,-1.d0, &
                                                                  -1.d0, 1.d0, 1.d0, &
                                                                  -1.d0, 1.d0, 0.d0, &
                                                                  -1.d0, 1.d0,-1.d0, &
                                                                  -1.d0, 0.d0, 1.d0, &
                                                                  -1.d0, 0.d0, 0.d0, &
                                                                  -1.d0, 0.d0,-1.d0, &
                                                                  -1.d0,-1.d0, 1.d0, &
                                                                  -1.d0,-1.d0, 0.d0, &
                                                                  -1.d0,-1.d0,-1.d0 /)
    ! translations
    real(kind=8), dimension(3,27), parameter :: trans = reshape(vtrans, (/3, 27/))

    contains

    subroutine cleanCrystal(cell)
        ! Clean crystal type :
        !
        !    * deallocate array if allocated
        !    * set default value to all attributes
        !

        implicit none

        type(Crystal), intent(inout) :: cell

        ! lattice vectors and parameters
        cell % veca(:) =  0.d0 ; cell % vecb(:) =  0.d0 ; cell % vecc(:) =  0.d0
        cell % a       = -1.d0 ; cell % b       = -1.d0 ; cell % c       = -1.d0
        cell % alpha   = -1.d0 ; cell % beta    = -1.d0 ; cell % gam     = -1.d0

        ! crystal name
        cell % cname = ""

        ! atom types
        cell % nType = 0
        if (allocated(cell % nAtomPerType)) deallocate(cell % nAtomPerType)
        if (allocated(cell % atomType)) deallocate(cell % atomType)

        ! atoms
        cell % Natoms = 0
        if (allocated(cell % atoms)) deallocate(cell % atoms)

        ! transformation matrix
        cell % Mcart2red(:,:) = 0.d0
        cell % Mred2cart(:,:) = 0.d0

    end subroutine cleanCrystal

    subroutine computeLatticeVectors(cell, verbose) 
        ! Compute lattice vectors from lattice parameters. Lattice parameters must
        ! be known. By default, lattice vector veca is along x, vecb is in the xOy
        ! plane and z is choosen such as [veca, vecb, vecc] is positive.

        implicit none

        ! arguments 
        type(Crystal), intent(inout)  :: cell    ! Crystal structure 
        logical, intent(in), optional :: verbose ! manage verbosity

        ! local :
        integer      :: i
        real(kind=8) :: alpha_rad, beta_rad, gamma_rad
        real(kind=8) :: cotgamma
        real(kind=8) :: coef2, coef3

        ! check
        if (cell % a < 0.) stop "lattice parameter a unknown"
        if (cell % b < 0.) stop "lattice parameter b unknown"
        if (cell % c < 0.) stop "lattice parameter c unknown"
        if (cell % alpha < 0.) stop "lattice parameter alpha unknown"
        if (cell % beta  < 0.) stop "lattice parameter beta unknown"
        if (cell % gam   < 0.) stop "lattice parameter gam unknown"

        if (present(verbose) .and. verbose .eqv. .true.) then
            write(*, *)
            write(*, "('# lattice parameters')")
            write(*, "('a      = ', ES10.3)") cell % a 
            write(*, "('b      = ', ES10.3)") cell % b
            write(*, "('c      = ', ES10.3)") cell % c
            write(*, "('alpha  = ', F10.2)") cell % alpha
            write(*, "('beta   = ', F10.2)") cell % beta
            write(*, "('gamma  = ', F10.2)") cell % gam
            write(*, *)
        end if

        ! radian conversion
        alpha_rad = cell % alpha * pi / 180.d0
        beta_rad  = cell % beta  * pi / 180.d0
        gamma_rad = cell % gam   * pi / 180.d0
        
        cotgamma = cos(gamma_rad) / sin(gamma_rad)

        ! matrix coefficients
        coef2 = (cos(alpha_rad) - cos(gamma_rad) * cos(beta_rad)) / sin(gamma_rad)
        coef3 = (sin(beta_rad))**2 - coef2**2             ! = 1 - cos^2( beta ) - coef2**2
        if( coef3 > 0. ) then
            coef3 = sqrt(coef3)
        else
            stop "subroutine computeLatticeVectors : bad matrix coefficient"
        end if

        ! vecteur unitaire de la maille dans le repere xyz
        cell % veca(1) = 1.d0
        cell % veca(2) = 0.d0
        cell % veca(3) = 0.d0

        cell % vecb(1) = cos(gamma_rad) 
        cell % vecb(2) = sin(gamma_rad)
        cell % vecb(3) = 0.d0

        cell % vecc(1) = cos(beta_rad)
        cell % vecc(2) = coef2
        cell % vecc(3) = coef3

        if (present(verbose) .and. verbose .eqv. .true.) then
            write(*,"('a ','(',F10.4,',',F10.4,',',F10.4,')')") (cell % veca(i), i = 1, 3)
            write(*,"('b ','(',F10.4,',',F10.4,',',F10.4,')')") (cell % vecb(i), i = 1, 3)
            write(*,"('c ','(',F10.4,',',F10.4,',',F10.4,')')") (cell % vecc(i), i = 1, 3)
        end if

    end subroutine computeLatticeVectors

    subroutine computeLatticeParameters(cell, verbose) 
        ! Compute lattice parameters from lattice vectors. Lattice vectors must be
        ! known.

        implicit none

        ! arguments 
        type(Crystal), intent(inout)  :: cell    ! crystal structure
        logical, intent(in), optional :: verbose ! manage verbosity

        ! local variables
        real(kind=8) :: scal

        cell % a = sqrt(sum((cell % veca(:))**2))
        cell % b = sqrt(sum((cell % vecb(:))**2))
        cell % c = sqrt(sum((cell % vecc(:))**2))

        ! compute alpha in degree
        scal = dot_product(cell % vecb(:) / cell % b, cell % vecc(:) / cell % c)
        if (abs(scal) < 1.) then
            cell % alpha = acos(scal) * 180.d0 / pi
        else
            write(*,"('alpha, dot product > 1., scal = ',G20.10)") scal
            stop
        end if

        ! compute beta in degree
        scal = dot_product(cell % veca(:) / cell % a, cell % vecc(:) / cell % c)
        if (abs(scal) < 1.) then
            cell % beta = acos(scal) * 180.d0 / pi
        else
            write(*,"('beta, dot product > 1., scal = ',G20.10)") scal
            stop
        end if

        ! compute gamma in degree
        scal = dot_product(cell % veca(:) / cell % a, cell % vecb(:) / cell % b)
        if (abs(scal) < 1.) then
            cell % gam = acos(scal) * 180.d0 / pi
        else
            write(*,"('gamma, dot product > 1., scal = ',G20.10)") scal
            stop
        end if

        ! print
        if (present(verbose) .and. verbose .eqv. .true.) then
            write(*, *)
            write(*, "('# lattice parameters')")
            write(*, "('a      = ', ES10.3)") cell % a 
            write(*, "('b      = ', ES10.3)") cell % b
            write(*, "('c      = ', ES10.3)") cell % c
            write(*, "('alpha  = ', F10.2)") cell % alpha
            write(*, "('beta   = ', F10.2)") cell % beta
            write(*, "('gamma  = ', F10.2)") cell % gam
            write(*, *)
        end if

    end subroutine computeLatticeParameters

    subroutine fillMatrix(cell, verbose)
        ! fill transformation matrix from lattice vectors.
        !
        !     * Mred2cart : reduce    -> cartesian
        !     * Mcart2red : cartesian -> reduce
        !
        ! Lattice vectors must be knwon and must be consistent with coordinates.

        implicit none

        ! arguments
        type(Crystal), intent(inout)  :: cell     ! crystal structure
        logical, intent(in), optional :: verbose  ! manage verbosity

        ! local variables
        integer      :: i, j
        real(kind=8), dimension(3,3) :: a
        !real(kind=8) :: alpha_rad, beta_rad, gamma_rad
        !real(kind=8) :: cotgamma
        !real(kind=8) :: coef2, coef3

        ! set matrix red2cart from lattice vectors
        cell % Mred2cart(:,1) = cell % veca(:)
        cell % Mred2cart(:,2) = cell % vecb(:)
        cell % Mred2cart(:,3) = cell % vecc(:)

        a(:,:) = cell % Mred2cart(:,:)
        call inverse(a, cell % Mcart2red, 3)

        ! check
        !if (cell % a < 0.) stop "lattice parameter a unknown"
        !if (cell % b < 0.) stop "lattice parameter b unknown"
        !if (cell % c < 0.) stop "lattice parameter c unknown"
        !if (cell % alpha < 0.) stop "lattice parameter alpha unknown"
        !if (cell % beta  < 0.) stop "lattice parameter beta unknown"
        !if (cell % gam   < 0.) stop "lattice parameter gam unknown"

        ! radian conversion
        !alpha_rad = cell % alpha * pi / 180.d0
        !beta_rad  = cell % beta  * pi / 180.d0
        !gamma_rad = cell % gam   * pi / 180.d0

        !cotgamma  = cos(gamma_rad) / sin(gamma_rad)

        ! matrix coefficients
        !coef2 = (cos(alpha_rad) - cos(gamma_rad) * cos(beta_rad)) / sin(gamma_rad)
        !coef3 = (sin(beta_rad))**2 - coef2**2             ! = 1 - cos^2( beta ) - coef2**2
        !if( coef3 > 0. ) then
            !coef3 = sqrt(coef3)
        !else
            !stop "subroutine fillMatrix : bad matrix coefficient"
        !end if
        
        ! transformation matrix reduce coordinates -> cartesian coordinates
        !cell % Mred2cart(:,:) = 0.d0
        !cell % Mred2cart(1,1) = cell % a                      ! bv 1

        !cell % Mred2cart(1,2) = cell % b * cos(gamma_rad)     ! bv 2
        !cell % Mred2cart(2,2) = cell % b * sin(gamma_rad)     ! bv 3

        !cell % Mred2cart(1,3) = cell % c * cos(beta_rad)      ! bv 4
        !cell % Mred2cart(2,3) = cell % c * coef2              ! bv 5
        !cell % Mred2cart(3,3) = cell % c * coef3              ! bv 6

        ! transformation matrix cartesian coordinates -> reduce coordinates
        !cell % Mcart2red(:,:) = 0.d0
        !cell % Mcart2red(1,1) = 1.d0 / cell % a
        !cell % Mcart2red(2,1) = - cotgamma / cell % a
        !cell % Mcart2red(3,1) = 1.d0 / cell % a * (coef2 / coef3 * cotgamma - cos(beta_rad) / coef3)

        !cell % Mcart2red(2,2) = 1.d0 / (cell % b * sin(gamma_rad))
        !cell % Mcart2red(3,2) = - coef2 / (cell % b * coef3 * sin(gamma_rad))

        !cell % Mcart2red(3,3) = 1.d0 / (cell % c * coef3)

        ! REMARK
        ! relations between Mcart2red et Mred2cart
        ! 
        ! t(Mred2cart^-1) = Mcart2red 
        ! Mred2cart * t(Mcart2red) = Id

        if (present(verbose) .and. verbose .eqv. .true.) then
          write(*,"('# reduce -> cartesian coordinates')")
          do i = 1, 3
            write(*,"(3F10.4)") (cell % Mred2cart(i,j), j = 1, 3)
          end do

          write(*,*)
          write(*,"('# cartesian -> reduce coordinates')")
          do i = 1, 3
            write(*,"(3F10.4)") (cell % Mcart2red(i,j), j = 1, 3)
          end do
          
          write(*,*)
        end if

    end subroutine fillMatrix

    function red2cart(cell, r)
        ! convert cartesian coordinate into reduce coordinate.
        implicit none
        type(Crystal), intent(in)              :: cell     ! crystal structure
        real(kind=8), dimension(3), intent(in) :: r        ! reduce coordinate
        real(kind=8), dimension(3)             :: red2cart ! cartesian  coordinate
        red2cart = matmul(cell % Mred2cart(:,:), r(:))
    end function red2cart

    function cart2red(cell, x)
        ! convert cartesian coordinate into reduce coordinate.
        implicit none
        type(Crystal), intent(in)              :: cell     ! crystal structure
        real(kind=8), dimension(3), intent(in) :: x        ! cartesian coordinate
        real(kind=8), dimension(3)             :: cart2red ! reduce coordinate
        cart2red = matmul(cell % Mcart2red(:,:), x(:))
    end function cart2red

    subroutine computeXYZCoord(cell)
        ! compute cartesian coordinates from reduce coordinates of all atoms into
        ! the crystal structure.
        implicit none
        type(Crystal), intent(inout) :: cell ! Crystal structure
        integer :: iat
        do iat = 1, cell % Natoms
            cell % atoms(iat) % x(:) = red2cart(cell, cell % atoms(iat) % r(:))
        end do
    end subroutine computeXYZCoord

    subroutine computeRedCoord(cell)
        ! compute reduce coordinates from cartesian coordinates of all atoms into
        ! the crystal structure.
        implicit none
        type(Crystal), intent(inout) :: cell ! Crystal structure
        integer :: iat
        do iat = 1, cell % Natoms
            cell % atoms(iat) % r(:) = cart2red(cell, cell % atoms(iat) % x(:))
        end do
    end subroutine computeRedCoord

    real(kind=8) function dist_i(cell, iat, jat)
        ! Compute the minimum distance between atoms iat and jat considering
        ! periodic boundary conditions. This function needs that reduce
        ! coordinates of atoms iat and jat are known.

        implicit none

        ! arguments
        type(Crystal), intent(in) :: cell ! Crystal structure
        integer, intent(in)       :: iat  ! atom number
        integer, intent(in)       :: jat  ! atom number

        ! local
        real(kind=8), dimension(3) :: dr, dx

        dr(:) = image(cell % atoms(iat) % r(:), cell % atoms(jat) % r(:))
        dx(:) = red2cart(cell, dr(:))

        dist_i = sqrt(sum(dx(:)**2))
        
    end function dist_i

   
    real(kind=8) function dist_r(cell, ri, rj)
        ! Compute the minimum distance between atoms at coordinate ri and rj
        ! considering periodic boundary conditions. ri and rj are reduce
        ! coordinates.

        implicit none

        ! arguments
        type(Crystal), intent(in)              :: cell ! Crystal structure
        real(kind=8), dimension(3), intent(in) :: ri   ! reduce coordinates of atom i
        real(kind=8), dimension(3), intent(in) :: rj   ! reduce coordinates of atom j

        ! local
        real(kind=8), dimension(3) :: dr, dx

        dr(:) = image(ri(:), rj(:))
        dx(:) = red2cart(cell, dr(:))

        dist_r = sqrt(sum(dx(:)**2))
        
    end function dist_r

    real(kind=8) function dist_x(cell, xi, xj)
        ! Compute the minimum distance between atoms at coordinate xi and xj
        ! considering periodic boundary conditions. xi and xj are cartesian
        ! coordinates.

        implicit none

        ! arguments
        type(Crystal), intent(in)              :: cell ! Crystal structure
        real(kind=8), dimension(3), intent(in) :: xi   ! cartesian coordinates of atom i
        real(kind=8), dimension(3), intent(in) :: xj   ! cartesian coordinates of atom j

        ! local
        real(kind=8), dimension(3) :: dr, dx

        dr(:) = image(cart2red(cell, xi(:)), cart2red(cell, xj(:)))
        dx(:) = red2cart(cell, dr(:))

        dist_x = sqrt(sum(dx(:)**2))
        
    end function dist_x

    function image(ri, rj)
        ! Returns the shortest vector in reduce coordinate considering the
        ! periodic boundary conditions. This function needs that reduce
        ! coordinate had been calculated. Use computeRedCoord if needed.

        implicit none

        ! arguments
        real(kind=8), dimension(3), intent(in) :: ri ! reduce coordinates of atom i
        real(kind=8), dimension(3), intent(in) :: rj ! reduce coordinates of atom j

        ! function
        real(kind=8), dimension(3) :: image

        ! local
        integer :: i

        image(:) = ri(:) - rj(:)
        do i = 1, 3, 1
            if (image(i) > 0.5) then
                image(i) = image(i) - 1.d0
            else if (image(i) < - 0.5) then
                image(i) = image(i) + 1.d0
            end if
        end do

    end function image

    real(kind=8) function computeVolume(cell)
        ! Compute the cell volume
    
        implicit none

        ! arguments :
        type(Crystal), intent(in) :: cell  ! Crystal structure

        ! locales :
        real(kind=8) :: noortho, cosa, cosb, cosg

        if (cell % a < 0. .or. cell % b < 0. .or. cell % c < 0. & 
           .or. cell % alpha < 0. .or. cell % beta < 0. .or. cell % gam < 0.) then
            stop "first, compute lattice paramters with subroutine computeLatticeParameters()"
        end if

        cosa = cos(cell % alpha * pi / 180.d0)
        cosb = cos(cell % beta  * pi / 180.d0)
        cosg = cos(cell % gam   * pi / 180.d0)

        noortho = 1.d0 - cosa**2 - cosb**2 - cosg**2 + 2.d0 * cosa * cosb * cosg

        if( noortho > 0.d0 ) then
          computeVolume = cell % a * cell % b * cell % c * sqrt( noortho )
        else
          stop "error computing the volume, check lattice vectors"
        end if

    end function computeVolume

    subroutine getComposition(cell)
        ! look for different atom type and number of each atom type.

        implicit none

        ! arguments
        type(Crystal), intent(inout) :: cell

        ! local
        integer, parameter                    :: ntypemax = 100
        character(len=3), dimension(ntypemax) :: tmpList
        integer                               :: iat, ntype, i
        logical                               :: newType

        ! look for type name and number of different type
        ntype = 1
        tmpList(1) = cell % atoms(1) % atname
        do iat = 2, cell % Natoms
            newType = .true.
            do i = 1, ntype
                if (cell % atoms(iat) % atname == tmpList(i)) then
                    newType = .false.
                    exit
                end if
            end do
            if (newtype) then
                ntype = ntype + 1
                if (ntype > ntypemax) stop "increase ntypemax in module crysmod"
                tmpList(ntype) = cell % atoms(iat) % atname
            end if
        end do

        ! allocate
        allocate(cell % nAtomPerType(nType), cell % atomType(nType))

        ! number of different type
        cell % nType = ntype

        ! names of types
        cell % atomType(:) = tmpList(1:ntype)

        ! look for the number of each type
        cell % nAtomPerType(:) = 0
        do iat = 1, cell % Natoms
            do i = 1, cell % nType
                if (cell % atoms(iat) % atname == cell % atomType(i)) then
                    cell % nAtomPerType(i) = cell % nAtomPerType(i) + 1
                end if
            end do
        end do

    end subroutine getComposition

    subroutine fromPOSCAR(cell, poscar, iunit, verbose)
        ! Build a crystal object from a POSCAR/CONTCAR VASP structure file.

        implicit none

        ! arguments
        type(Crystal), intent(out)             :: cell    ! built structure 
        integer, intent(in), optional          :: iunit   ! file unit for reading
        character(len=*), intent(in), optional :: poscar  ! file name from which the structure is read (default = "POSCAR")
        logical, intent(in), optional          :: verbose ! manage verbosity (default = .False.)

        ! local
        integer            :: i, j, n, io
        real(kind=8)       :: scal
        character(len=100) :: poscarname, title, line
        character(len=1)   :: c
        logical            :: existe, vasp4, cartesian, verb

        ! check POSCAR file
        if (.not. present(poscar)) then
            poscarname = "POSCAR"
        else
            poscarname = poscar
        end if

        inquire(file = trim(poscarname), exist = existe)
        if (.not. existe) then
            write(*,"('File ', a, ' does not exist')") trim(adjustl(poscarname))
            stop
        end if

        ! verbosity
        if (.not. present(verbose)) then
            verb = .false.
        else
            verb = verbose
        end if

        ! file unit
        if (present(iunit)) then
            io = iunit
        else
            io = 64
        end if

        ! open POSCAR file
        open(unit = io, file = trim(poscarname), action = "read")

        ! read header
        read(io, "(a)") title
        cell % cname = trim(adjustl(title))
        read(io, *) scal
        read(io, *) (cell % veca(i), i = 1, 3)
        read(io, *) (cell % vecb(i), i = 1, 3)
        read(io, *) (cell % vecc(i), i = 1, 3)
        
        cell % veca(:) = scal * cell % veca(:)
        cell % vecb(:) = scal * cell % vecb(:)
        cell % vecc(:) = scal * cell % vecc(:)

        call computeLatticeParameters(cell, verb)
        call fillMatrix(cell, verb)

        ! read line 6 : either number of ions or their type
        read(io, "(a)") line
        read(line, *) c
        if (.not. (ichar(c) >= 48 .and. ichar(c) <= 57)) then
            ! vasp 5 : ions type and number
            vasp4 = .false.
            ! number of element on the line
            line = adjustl(line)
            n = 1
            do i = 2, len_trim(line)
                if (line(i:i) == " " .and. line(i-1:i-1) /= " ") n = n + 1
            end do
            cell % nType = n
            allocate(cell % atomType(1:cell % nType), cell % nAtomPerType(1:cell % nType))
            read(line, *) (cell % atomType(i), i = 1, cell % nType)
            read(io, *) (cell % nAtomPerType(i), i = 1, cell % nType)
        else
            ! vasp 4 : ions number
            vasp4 = .true.
            ! number of element on the line
            line = adjustl(line)
            n = 1
            do i = 2, len_trim(line)
                if (line(i:i) == " " .and. line(i-1:i-1) /= " ") n = n + 1
            end do
            cell % nType = n
            allocate(cell % atomType(1:cell % nType), cell % nAtomPerType(1:cell % nType))
            read(line, *) (cell % nAtomPerType(i), i = 1, cell % nType)
            ! unknwon atom name => use Xi instead
            do i = 1, cell % nType
                write(cell % atomType(i), "('X', i2.2)") i
            end do
        end if

        ! compute the number of atom and allocate atom table
        cell % Natoms = 0
        do i = 1, cell % nType
            cell % Natoms = cell % Natoms + cell % nAtomPerType(i)
        end do
        allocate(cell % atoms(1:cell % Natoms))

        ! set atom name
        n = 0
        do i = 1, cell % nType
            do j = 1, cell % nAtomPerType(i)
                n = n + 1
                cell % atoms(n) % atname = cell % atomType(i)
            end do
        end do

        ! read coordinates
        read(io, "(a)") line
        read(line, *) c
        if (c == "s" .or. c == "S") then
            ! selective Dynamics
            read(io, "(a)") line
            read(line, *) c
        end if
        if (c == "c" .or. c == "C" .or. c == "k" .or. c == "K") then
            ! cartesian coordinates
            cartesian = .true.
        else if (c == "d" .or. c == "D") then
            ! reduce coordinates
            cartesian = .false.
        else
            stop "coordinates"
        end if

        ! read coordinates
        do i = 1, cell % Natoms
            if (cartesian) then
                read(io, *) (cell % atoms(i) % x(j), j = 1, 3)
                cell % atoms(i) % x(:) = scal * cell % atoms(i) % x(:)
                cell % atoms(i) % r(:) = 0.d0
            else
                read(io, *) (cell % atoms(i) % r(j), j = 1, 3)
                cell % atoms(i) % x(:) = 0.d0
            end if
        end do

        close(io)

        ! WARNING : case where a is not along x for cartesian coordinates ????

        ! compute reduce or cartesian coordinate
        if (cartesian) then
            call computeRedCoord(cell)
        else
            call computeXYZCoord(cell)
        end if

    end subroutine fromPOSCAR

    subroutine toPOSCAR(cell, filename, iunit)
        ! print cell into a VASP structure file POSCAR

        implicit none

        ! arguments
        type(Crystal), intent(inout)           :: cell      ! crystal structure
        integer, intent(in), optional          :: iunit     ! file logical unit (default : 64)
        character(len=*), intent(in), optional :: filename  ! file name (default : POSCAR)

        ! local
        integer :: io, i, iat, nat, j

        ! open output file on unit io
        if (present(iunit)) then
            io = iunit
        else 
            io = 10
        end if

        ! open output file
        if (present(filename)) then
            open(unit = io, file = trim(adjustl(filename)))
        else
            open(unit = io, file = "POSCAR")
        end if

        if (cell % cname == "") then
            write(io, "('crystal structure')")
        else
            write(io, "(a)") cell % cname
        end if
        write(io, "(' 1.0')")

        ! lattice vectors
        write(io, "(3F20.12)") (cell % veca(i), i = 1, 3)
        write(io, "(3F20.12)") (cell % vecb(i), i = 1, 3)
        write(io, "(3F20.12)") (cell % vecc(i), i = 1, 3)

        ! write composition
        call getComposition(cell)
        write(io, "(100a6)") (cell % atomType(i), i = 1, cell % nType)
        write(io, "(100i6)") (cell % nAtomPerType(i), i = 1, cell % nType)

        ! direct keyword
        write(io, "('direct')")

        ! sort atoms by name
        nat = 0
        do i = 1, cell % nType
            do iat = 1, cell % Natoms
                if (cell % atomType(i) == cell % atoms(iat) % atname) then
                    nat = nat + 1
                    write(io,"(3F20.12)") (cell % atoms(iat) % r(j), j = 1, 3)
                end if
            end do
        end do

        close(unit = io)

    end subroutine toPOSCAR

    subroutine fromCONFIG(cell, config, iunit, verbose)
        ! Build a crystal object from a CONFIG/REVCON DLPOLY structure file.

        implicit none

        ! arguments
        type(Crystal), intent(out)             :: cell    ! built structure 
        integer, intent(in), optional          :: iunit   ! file unit for reading
        character(len=*), intent(in), optional :: config  ! file name from which the structure is read (default = "CONFIG")
        logical, intent(in), optional          :: verbose ! manage verbosity (default = .False.)

        ! local
        integer            :: i, j, io, nat, levcfg, imcons
        character(len=100) :: configname, title
        logical            :: existe, verb

        ! CONFIG name
        if (.not. present(config)) then
            configname = "CONFIG"
        else
            configname = config
        end if

        ! check CONFIG file
        inquire(file = trim(configname), exist = existe)
        if (.not. existe) then
            write(*,"('File ', a, ' does not exist')") trim(adjustl(configname))
            stop
        end if

        ! verbosity
        if (.not. present(verbose)) then
            verb = .false.
        else
            verb = verbose
        end if

        ! file unit
        if (present(iunit)) then
            io = iunit
        else
            io = 64
        end if

        ! open CONFIG file
        open(unit = io, file = trim(configname), action = "read")

        ! read header
        read(io, "(a)") title
        cell % cname = trim(adjustl(title))
        read(io, *, iostat = j) levcfg, imcons, nat
        if (j /= 0) then
            write(*, "(/,'error maybe the number of atoms is missing, line 2 on')")
            write(*, "('CONFIG file')")
            stop
        end if
        read(io, *) (cell % veca(i), i = 1, 3)
        read(io, *) (cell % vecb(i), i = 1, 3)
        read(io, *) (cell % vecc(i), i = 1, 3)
        
        call computeLatticeParameters(cell, verb)
        call fillMatrix(cell, verb)

        ! allocate atoms array
        if (verb) write(*, "('Atom number = ', i6)") nat
        cell % Natoms = nat
        allocate(cell % atoms(1:nat))

        ! read coordinates
        do i = 1, nat
            read(io, *) cell % atoms(i) % atname
            read(io, *) (cell % atoms(i) % x(j), j = 1, 3)
            do j = 1, levcfg, 1
                read(io, *)
            end do
        end do
        close(io)

        ! compute reduce coordinates
        call computeRedCoord(cell)

    end subroutine fromCONFIG

    subroutine neighbor(cell, rcut, ns, nn, iatom, neigh, nneigh, ligand)
        ! look for neighbors of atoms in the list iatom. If ligand present, only
        ! atoms with name ligand are considered. If rcut < cell parameters / 2,
        ! the subroutine looks for multiple neighbors by considering all atoms
        ! in the 26 adjacent cells around the central the considered cell.
        !

        implicit none

        ! arguments
        type(Crystal), intent(in)               :: cell   ! crystal structure
        real(kind=8), intent(in)                :: rcut   ! cut off radius (in Angstrom)
        integer, intent(in)                     :: ns     ! number of selected atom for which neighbors list is built
        integer, intent(in)                     :: nn     ! maximum number of neighbors   
        integer, dimension(ns), intent(in)      :: iatom  ! list of atom index for which neighbors list is built
        integer, dimension(ns, nn), intent(out) :: neigh  ! neighbors list
        integer, dimension(ns), intent(out)     :: nneigh ! number of neighbors of each atom in iatom list
        character(len=3), intent(in), optional  :: ligand ! limit the neighbors list to atoms name ligand (defaut = "#" not any atom is excluded)

        ! local
        integer                    :: is, iat, it, jat
        real(kind=8)               :: d
        real(kind=8), dimension(3) :: rj, dr
        character(len=3)           :: lig
        
        if (present(ligand)) then
            lig = trim(adjustl(ligand))
        else
            lig = "#  "
        end if

        neigh(:,:) = 0
        nneigh(:) = 0

        if (rcut > cell % a / 2.d0 .or. &
            rcut > cell % b / 2.d0 .or. &
            rcut > cell % c / 2.d0) then
            ! small cell look for multiple neighbors due to translation
            write(*,"('cell parameters / 2 < rcut => look for multiple neighbors')")
            do is = 1, ns
                iat = iatom(is)
                do jat = 1, cell % Natoms
                    if (iat == jat) cycle
                    if (present(ligand) .and. trim(adjustl(cell % atoms(jat) % atname)) /= lig) cycle
                    rj(:) = cell % atoms(jat) % r(:)

                    do it = 1, 27
                        dr(:) = rj(:) + trans(:,it) - cell % atoms(iat) % r(:)
                        d = sqrt(sum((red2cart(cell, dr(:)))**2))
                        if (d <= rcut) then
                            nneigh(is) = nneigh(is) + 1
                            if (nneigh(is) > nn) stop "too much neighbors, increase nn"
                            neigh(is, nneigh(is)) = jat
                        end if
                    end do
                end do
            end do

        else
            ! simple case where two atoms can be neighbor only one time
            do is = 1, ns
                iat = iatom(is)
                do jat = 1, cell % Natoms
                    if (jat == iat) cycle
                    if (present(ligand) .and. trim(adjustl(cell % atoms(jat) % atname)) /= lig) cycle

                    d = dist_i(cell, iat, jat)
                    if (d <= rcut) then
                        nneigh(is) = nneigh(is) + 1
                        if (nneigh(is) > nn) stop "too much neighbors, increase nn"
                        neigh(is, nneigh(is)) = jat
                    end if
                end do
            end do
        end if

    end subroutine neighbor

    subroutine allneighbor(cell, rcut, nn, neigh, nneigh)
        ! Look for neighbors of all atoms in cell. If rcut > cell parameters/2,
        ! the subroutine looks for multiple neighbors by considering all atoms
        ! in the 26 adjacent cells around the considered cell.

        implicit none

        ! arguments
        type(Crystal), intent(in)                        :: cell   ! crystal structure
        real(kind=8), intent(in)                         :: rcut   ! cut off radius
        integer, intent(in)                              :: nn     ! maximum number of neighbors
        integer, dimension(cell%Natoms, nn), intent(out) :: neigh  ! neighbors list
        integer, dimension(cell%Natoms), intent(out)     :: nneigh ! number of neighbors of each atom

        ! local
        integer                    :: iat, jat, it
        real(kind=8)               :: d
        real(kind=8), dimension(3) :: dr
        
        neigh(:,:) = 0
        nneigh(:) = 0

        if (rcut > cell % a / 2.d0 .or. &
            rcut > cell % b / 2.d0 .or. &
            rcut > cell % c / 2.d0) then
            ! small cell look for multiple neighbors
            write(*,"('cell parameters / 2 < rcut => look for multiple neighbors')")
            do iat = 1, cell % Natoms - 1
                do jat = iat + 1, cell % Natoms
                    do it = 1, 27
                        dr(:) = cell % atoms(jat) % r(:) + trans(:,it) - cell % atoms(iat) % r(:)
                        d = sqrt(sum((red2cart(cell, dr(:)))**2))
                        if (d <= rcut) then
                            nneigh(iat) = nneigh(iat) + 1
                            nneigh(jat) = nneigh(jat) + 1
                            if (nneigh(iat) > nn .or. nneigh(jat) > nn) stop "too much neighbors, increase nn"
                            neigh(iat, nneigh(iat)) = jat
                            neigh(jat, nneigh(jat)) = iat
                        end if
                    end do
                end do
            end do

        else
            ! simple case where two atoms can be neighbor only one time
            do iat = 1, cell % Natoms - 1
                do jat = iat + 1, cell % Natoms
                    d = dist_i(cell, iat, jat)
                    if (d <= rcut) then
                        nneigh(iat) = nneigh(iat) + 1
                        nneigh(jat) = nneigh(jat) + 1
                        if (nneigh(iat) > nn .or. nneigh(jat) > nn) stop "too much neighbors, increase nn"
                        neigh(iat, nneigh(iat)) = jat
                        neigh(jat, nneigh(jat)) = iat
                    end if
                end do
            end do

        end if

    end subroutine allneighbor

    subroutine inverse(a,c,n)
    ! Inverse matrix
    !
    ! :Method: Based on Doolittle LU factorization for Ax=b
    !
    ! :author: Alex G. December 2009
    !
    ! :comments: the original matrix a(n,n) will be destroyed  during the calculation

    implicit none 

    ! arguments 
    integer, intent(in)                         :: n ! dimension
    real(kind=8), dimension(n,n), intent(inout) :: a ! array of coefficients for matrix A 
    real(kind=8), dimension(n,n), intent(out)   :: c ! inverse matrix of a

    ! local
    integer                      :: i, j, k
    real(kind=8)                 :: coeff
    real(kind=8), dimension(n,n) :: L, U
    real(kind=8), dimension(n)   :: b, d, x

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
      ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
      ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
    end subroutine inverse
end module crysmod
