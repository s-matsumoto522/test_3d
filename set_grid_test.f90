module grid_test
    implicit none
    integer, parameter :: NXmin = -20, NXmax = 20     !x方向の計算領域の形状
    integer, parameter :: NYmin = -20, NYmax = 20     !y方向の計算領域の形状
    integer, parameter :: NZmin = -20, NZmax = 20     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 20.0d0, Ymax = 20.0d0, Zmax = 20.0d0  !各方向の計算領域の最大値
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ    !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2 !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
contains
!************************************
!   格子点を出力するサブルーチン    *
!************************************
    subroutine output_grid(X, Y, Z)
        double precision, intent(in) :: X(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision, intent(in) :: Y(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision, intent(in) :: Z(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        integer iX, iY, iZ
        open(11, file = 'set_grid_test.dat')
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    write(11, *) X(iX, iY, iZ), Y(iX, iY, iZ), Z(iX, iY, iZ)
                enddo
            enddo
        enddo
    end subroutine output_grid
!********************************
!   格子を設定するサブルーチン  *
!********************************
    subroutine set_grid(X, Y, Z)
        double precision, intent(out) :: X(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision, intent(out) :: Y(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision, intent(out) :: Z(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        integer iX, iY, iZ
        !---各方向の刻み幅を計算---
        dX = Xmax / dble(NXmax)
        dY = Ymax / dble(NYmax)
        dZ = Zmax / dble(NZmax)
        !---計算用数値の設定---
        ddX = 1.0d0 / dX
        ddY = 1.0d0 / dY
        ddZ = 1.0d0 / dZ
        ddX2 = 1.0d0 / (dX**2)
        ddY2 = 1.0d0 / (dY**2)
        ddZ2 = 1.0d0 / (dZ**2)
        !---各方向の計算領域の最小値を計算---
        Xmin = dX*dble(NXmin)
        Ymin = dY*dble(NYmin)
        Zmin = dZ*dble(NZmin)
        !---格子点の設定---
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    X(iX, iY, iZ) = dX*dble(iX)
                    Y(iX, iY, iZ) = dY*dble(iY)
                    Z(iX, iY, iZ) = dZ*dble(iZ)
                enddo
            enddo
        enddo
    end subroutine set_grid
end module grid_test

program main
    use grid_test
    implicit none
    double precision X(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
    double precision Y(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
    double precision Z(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
    call set_grid(X, Y, Z)
    call output_grid(X, Y, Z)
end program main

