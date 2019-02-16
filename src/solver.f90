module solverND
    use util
    implicit none
    
    contains
    ! u(1) := t
    ! u(2) := t'
    ! f    := (l_1 sin t - l_2 cos t) / 2.0
    !       = (l_1 sin u(1) - l_2 cos u(1)) / 2.0
    !
    ! u'(1) = u(2)
    ! u'(2) = f(s, u(1))
    !
    ! u' = F(u)

    subroutine nsolve(L, t_0, Dt_0, l_1, l_2, n, nt_L, nsint, ncost)
        real(8),    intent(in)  :: L, t_0, Dt_0, l_1, l_2
        integer(8), intent(in)  :: n

        real(8),    intent(out) :: nt_L
        real(8), dimension(0:n), intent(out)    :: nsint, ncost

        real(8), dimension(2,0:n)   :: u
        real(8), dimension(2)       :: k_1, k_2, k_3, k_4
        real(8) :: h
        integer(8) :: i

        h      = L / dble(n)

        u(1,0) = t_0
        u(2,0) = Dt_0
        ncost(0) = cos(u(1,i))
        nsint(0) = sin(u(1,i))

        do i = 0, n-1
            k_1(1) = h * u(2,i)
            k_1(2) = h * (l_1 * sin( u(1,i) )               - l_2 * cos( u(1,i) )) / 2.0

            k_2(1) = h * (u(2,i) + k_1(2)/2.0)
            k_2(2) = h * (l_1 * sin( u(1,i) + k_1(1)/2.0 )  - l_2 * cos( u(1,i) + k_1(1)/2.0 )) / 2.0

            k_3(1) = h * (u(2,i) + k_2(2)/2.0)
            k_3(2) = h * (l_1 * sin( u(1,i) + k_2(1)/2.0 )  - l_2 * cos( u(1,i) + k_2(1)/2.0 )) / 2.0

            k_4(1) = h * (u(2,i) + k_2(2))
            k_4(2) = h * (l_1 * sin( u(1,i) + k_3(1) )      - l_2 * cos( u(1,i) + k_3(1) )) / 2.0

            u(1,i+1) = u(1,i) + ( k_1(1) + 2.0*k_2(1) + 2.0*k_3(1) + k_4(1) ) / 6.0
            u(2,i+1) = u(2,i) + ( k_1(2) + 2.0*k_2(2) + 2.0*k_3(2) + k_4(2) ) / 6.0

            ncost(i+1) = cos(u(1,i+1))
            nsint(i+1) = sin(u(1,i+1))
        end do

        nt_L = u(1,n)
    end subroutine nsolve

    subroutine correct_vars(L, x_L, y_L, t_0, t_L, n, &
                            Dt_0, l_1, l_2, &
                            dDt_0, dl_1, dl_2, &
                            uDt_0, ul_1, ul_2, &
                            nt_L, nx_L, ny_L, &
                            nsint, ncost)
        real(8),    intent(in)  :: L, x_L, y_L, t_0, t_L
        integer(8), intent(in)  :: n
        real(8),    intent(in)  :: Dt_0, l_1, l_2
        real(8),    intent(in)  :: dDt_0, dl_1, dl_2

        real(8), intent(out)    :: uDt_0, ul_1, ul_2
        real(8), intent(out)    :: nt_L, nx_L, ny_L

        real(8), dimension(0:n), intent(out)  :: nsint, ncost

        real(8), dimension(0:n) :: nsint_tmp1, nsint_tmp2, ncost_tmp1, ncost_tmp2

        real(8) :: nt_L_tmp1, nt_L_tmp2
        real(8) :: h
        real(8), dimension(3,3) :: J  ! Finite Difference Jacobian
        real(8), dimension(3)   :: Delta_func, Delta_vars

        h = L / dble(n)

        call nsolve(L, t_0, Dt_0, l_1, l_2, n, nt_L, nsint, ncost)
        nx_L = definite_integral(ncost, n, h)
        ny_L = definite_integral(nsint, n, h)

        call nsolve(L, t_0, Dt_0 + dDt_0/2.0, l_1, l_2, n, nt_L_tmp1, nsint_tmp1, ncost_tmp1)
        call nsolve(L, t_0, Dt_0 - dDt_0/2.0, l_1, l_2, n, nt_L_tmp2, nsint_tmp2, ncost_tmp2)
        J(1, 1) = ( nt_L_tmp1 - nt_L_tmp2 ) / dDt_0
        J(2, 1) = ( definite_integral(ncost_tmp1, n, h) - definite_integral(ncost_tmp2, n, h) ) / dDt_0
        J(3, 1) = ( definite_integral(nsint_tmp1, n, h) - definite_integral(nsint_tmp2, n, h) ) / dDt_0

        call nsolve(L, t_0, Dt_0, l_1 + dl_1/2.0, l_2, n, nt_L_tmp1, nsint_tmp1, ncost_tmp1)
        call nsolve(L, t_0, Dt_0, l_1 - dl_1/2.0, l_2, n, nt_L_tmp2, nsint_tmp2, ncost_tmp2)
        J(1, 2) = ( nt_L_tmp1 - nt_L_tmp2 ) / dl_1
        J(2, 2) = ( definite_integral(ncost_tmp1, n, h) - definite_integral(ncost_tmp2, n, h) ) / dl_1
        J(3, 2) = ( definite_integral(nsint_tmp1, n, h) - definite_integral(nsint_tmp2, n, h) ) / dl_1

        call nsolve(L, t_0, Dt_0, l_1, l_2 + dl_2/2.0, n, nt_L_tmp1, nsint_tmp1, ncost_tmp1)
        call nsolve(L, t_0, Dt_0, l_1, l_2 - dl_2/2.0, n, nt_L_tmp2, nsint_tmp2, ncost_tmp2)
        J(1, 3) = ( nt_L_tmp1 - nt_L_tmp2 ) / dl_2
        J(2, 3) = ( definite_integral(ncost_tmp1, n, h) - definite_integral(ncost_tmp2, n, h) ) / dl_2
        J(3, 3) = ( definite_integral(nsint_tmp1, n, h) - definite_integral(nsint_tmp2, n, h) ) / dl_2

        Delta_func(1) = nt_L - t_L
        Delta_func(2) = nx_L - x_L
        Delta_func(3) = ny_L - y_L

        Delta_vars = matmul(inv3(J), Delta_func)

        uDt_0   = Dt_0 - Delta_vars(1)
        ul_1    = l_1  - Delta_vars(2)
        ul_2    = l_2  - Delta_vars(3)
    end subroutine correct_vars

    subroutine find_solution(L, x_L, y_L, t_0, t_L, &
                                n, &
                                Dt_0, l_1, l_2, &
                                dDt_0, dl_1, dl_2, &
                                et_L, ex_L, ey_L, &
                                nDt_0, nl_1, nl_2, &
                                x, y, &
                                convergence)
        real(8),    intent(in)  :: L, x_L, y_L, t_0, t_L
        integer(8), intent(in)  :: n
        real(8),    intent(in)  :: Dt_0, l_1, l_2
        real(8),    intent(in)  :: dDt_0, dl_1, dl_2
        real(8),    intent(in)  :: et_L, ex_L, ey_L

        real(8), intent(out)    :: nDt_0, nl_1, nl_2
        real(8), dimension(0:n), intent(out)  :: x, y
        

        logical(8), intent(out) :: convergence

        real(8), dimension(0:n) :: nsint, ncost

        real(8)     :: h

        real(8)     :: unDt_0, unl_1, unl_2

        real(8)     :: nt_L, nx_L, ny_L

        integer(8)  :: i

        h = L / dble(n)
        convergence = .false.

        nDt_0 = Dt_0
        nl_1  = l_1
        nl_2  = l_2

        do i = 1, 30
            call correct_vars(L, x_L, y_L, t_0, t_L, n, &
                                nDt_0, nl_1, nl_2, &
                                dDt_0, dl_1, dl_2, &
                                unDt_0, unl_1, unl_2, &
                                nt_L, nx_L, ny_L, &
                                nsint, ncost)

            nDt_0   = unDt_0
            nl_1    = unl_1
            nl_2    = unl_2

            if (abs(nt_L - t_L) <= et_L .and. &
                abs(nx_L - x_L) <= ex_L .and. &
                abs(ny_L - y_L) <= ey_L) then
                convergence = .true.
                exit
            end if
        end do
        
        if (convergence) then
            x = integral(ncost, n, h)
            y = integral(nsint, n, h)
        end if
    end subroutine find_solution
end module
