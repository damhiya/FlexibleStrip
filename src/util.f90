module util
    implicit none
    
    contains
    pure function integral(func, n, h) result(Ifunc)
        real(8), dimension(0:n)               :: Ifunc

        real(8), dimension(0:n),    intent(in)  :: func
        integer(8),                 intent(in)  :: n
        real(8),                    intent(in)  :: h

        integer(8)  :: i

        Ifunc(0) = 0.0

        do i = 1, n
            Ifunc(i) = Ifunc(i-1) + (func(i-1) + func(i))/2.0
        end do

        Ifunc = Ifunc * h
    end function integral

    pure function definite_integral(func, n, h) result(Ifunc)
        real(8) :: Ifunc

        real(8), dimension(0:n),    intent(in)  :: func
        integer(8),                 intent(in)  :: n
        real(8),                    intent(in)  :: h

        integer(8)  :: i

        Ifunc = 0.0
        do i = 0, n
            Ifunc = Ifunc + func(i)
        end do
        
        Ifunc = Ifunc - (func(0) + func(n))/2.0

        Ifunc = Ifunc * h
    end function definite_integral

    function inv3(A) result(Ainv)
        real(8), dimension(3,3), intent(in) :: A

        real(8), dimension(3,3) :: Acof ! cofactor
        real(8), dimension(3,3) :: Ainv ! inverse
        real(8) :: det               ! determinant

        Acof(1,1) = + ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
        Acof(2,1) = - ( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
        Acof(3,1) = + ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
        
        Acof(1,2) = - ( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
        Acof(2,2) = + ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
        Acof(3,2) = - ( A(1,1)*A(3,2) - A(1,2)*A(3,1) )

        Acof(1,3) = + ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
        Acof(2,3) = - ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
        Acof(3,3) = + ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

        det =  A(1,1)*Acof(1,1) + A(1,2)*Acof(2,1) + A(1,3)*Acof(3,1)

        Ainv = Acof / det
    end function inv3
end module
