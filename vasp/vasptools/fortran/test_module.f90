program testmod

    use crysmod

    implicit none

    type(Crystal) :: struct
    integer :: iat, jat

    call fromPOSCAR(struct, "POSCAR", verbose = .true.)

    do iat = 1, struct % Natoms - 1
        do jat = iat + 1, struct % Natoms
            write(*,"(a3,i2,2x,a3,i2,F8.3)") struct % atoms(iat) % atname, iat, &
                                             struct % atoms(jat) % atname, jat, &
                                             dist_i(struct, iat, jat)
        end do
    end do


end program testmod
