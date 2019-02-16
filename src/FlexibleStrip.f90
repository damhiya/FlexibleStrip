program FlexibleStrip
    use solverND
    implicit none

    real(8)     :: L, x_L, y_L, t_0, t_L

    integer(8)  :: n
    real(8)     :: dDt_0, dl_1, dl_2
    real(8)     :: et_L, ex_L, ey_L

    real(8)     :: Dt_0_min, Dt_0_max, l_1_min, l_1_max, l_2_min, l_2_max
    integer(8)  :: Dt_0_num, l_1_num, l_2_num

    real(8)     :: nDt_0, nl_1, nl_2
    real(8), dimension(:), allocatable  :: x, y

    real(8)     :: Delta_Dt_0, Delta_l_1, Delta_l_2
    real(8)     :: Dt_0, l_1, l_2
    integer(8)  :: i_Dt_0, i_l_1, i_l_2
    logical(8)  :: convergence

    real(8), dimension(:,:), allocatable :: result_list
    integer(8)  :: result_num

    integer(8)  :: i, j

    1 format (F20.10)
    2 format (I20)

    print '("# ", A, /)', "FlexibleStrip Solver"

    read 1, L
    read 1, x_L
    read 1, y_L
    read 1, t_0
    read 1, t_L

    read 2, n

    read 1, dDt_0
    read 1, dl_1
    read 1, dl_2

    read 1, et_L
    read 1, ex_L
    read 1, ey_L

    read 1, Dt_0_min
    read 1, Dt_0_max
    read 2, Dt_0_num

    read 1, l_1_min
    read 1, l_1_max
    read 2, l_1_num

    read 1, l_2_min
    read 1, l_2_max
    read 2, l_2_num

    print '("# ", A)', "Input"
    print '("# L        : ", F20.10)',    L
    print '("# x_L      : ", F20.10)',    x_L
    print '("# y_L      : ", F20.10)',    y_L
    print '("# t_0      : ", F20.10)',    t_0
    print '("# t_L      : ", F20.10)',    t_L

    print '("# n        : ", I20)',       n

    print '("# dDt_0    : ", F20.10)',    dDt_0
    print '("# dl_1     : ", F20.10)',    dl_1
    print '("# dl_2     : ", F20.10)',    dl_2

    print '("# et_L     : ", F20.10)',    et_L
    print '("# ex_L     : ", F20.10)',    ex_L
    print '("# ey_L     : ", F20.10)',    ey_L

    print '("# Dt_0_min : ", F20.10)',    Dt_0_min
    print '("# Dt_0_max : ", F20.10)',    Dt_0_max
    print '("# Dt_0_num : ", I20)',       Dt_0_num

    print '("# l_1_min  : ", F20.10)',    l_1_min
    print '("# l_1_max  : ", F20.10)',    l_1_max
    print '("# l_1_num  : ", I20)',       l_1_num

    print '("# l_2_min  : ", F20.10)',    l_2_min
    print '("# l_2_max  : ", F20.10)',    l_2_max
    print '("# l_2_num  : ", I20)',       l_2_num

    if (Dt_0_num == 1) then
        Delta_Dt_0 = 0.0
    else
        Delta_Dt_0 = (Dt_0_max - Dt_0_min) / dble(Dt_0_num - 1)
    end if

    if (l_1_num  == 1) then
        Delta_l_1  = 0.0
    else
        Delta_l_1  = (l_1_max  - l_1_min)  / dble(l_1_num  - 1)
    end if

    if (l_2_num  == 1) then
        Delta_l_2  = 0.0
    else
        Delta_l_2  = (l_2_max  - l_2_min)  / dble(l_2_num  - 1)
    end if

    Dt_0 = Dt_0_min
    l_1  = l_1_min
    l_2  = l_2_min

    allocate(x(0:n), y(0:n), result_list(3, Dt_0_num * l_1_num * l_2_num))

    result_num = 0

    do i_Dt_0 = 1, Dt_0_num
        do i_l_1 = 1, l_1_num
            Loop_l_2: do i_l_2 = 1, l_2_num
                call find_solution(L, x_L, y_L, t_0, t_L, &
                        n, &
                        Dt_0, l_1, l_2, &
                        dDt_0, dl_1, dl_2, &
                        et_L, ex_L, ey_L, &
                        nDt_0, nl_1, nl_2, &
                        x, y, &
                        convergence)

                if (.not. convergence) exit Loop_l_2

                do j = 1, result_num
                    if (abs(nDt_0 - result_list(1, j)) <= dDt_0 .and. &
                        abs(nl_1  - result_list(2, j)) <= dl_1  .and. &
                        abs(nl_2  - result_list(3, j)) <= dl_2) exit Loop_l_2
                end do

                result_num = result_num + 1

                result_list(1, result_num) = nDt_0
                result_list(2, result_num) = nl_1
                result_list(3, result_num) = nl_2

                print '(/,"# ", "Result", I20)', result_num
                print '("# Dt_0  :", F20.10)', nDt_0
                print '("# l_1   :", F20.10)', nl_1
                print '("# l_2   :", F20.10)', nl_2
                print '("# ", A)', "Curve :"
                do i = 0, n
                    print '(2F20.10)', x(i), y(i)
                end do

                l_2 = l_2 + Delta_l_2
            end do Loop_l_2

            l_1 = l_1 + Delta_l_1
        end do

        Dt_0 = Dt_0 + Delta_Dt_0
    end do

    deallocate(x, y, result_list)
    stop
end program