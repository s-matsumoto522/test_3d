module timestep_test
    implicit none
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = -10, NXmax = 10     !x方向の計算領域の形状
    integer, parameter :: NYmin = -10, NYmax = 10     !y方向の計算領域の形状
    integer, parameter :: NZmin = -10, NZmax = 10     !z方向の計算領域の形状
    integer, parameter :: Nstep = 500             !総ステップ数
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi, Zmax = 2.0d0*pi  !各方向の計算領域の最大値
    double precision, parameter :: NU = 1.0d0       !動粘性係数
    double precision, parameter :: dt = 1.0d-3/16.0d0      !時間の刻み幅
    integer, save :: Ng                             !格子点数
    double precision, save :: ddt                   !時間の刻み幅の逆数
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ    !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2 !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
contains
!************************************
!   時間積分を計算するサブルーチン  *
!************************************
    subroutine cal_timeint(V, R1, R2)
        double precision, intent(out) :: V
        double precision R1, R2
        V = V + 0.5d0*dt*(3.0d0*R2 - R1)
    end subroutine cal_timeint
end module timestep_test

program main
    use timestep_test
    implicit none
    integer istep
    double precision V, error, R1, R2
    V = 0.0d0
    do istep = 1, Nstep
        R1 = exp(-dt*dble(istep-1))
        R2 = exp(-dt*dble(istep))
        call cal_timeint(V, R1, R2)
    enddo
    error = abs(V - (1.0d0 - exp(-dt*dble(Nstep))))
    write(*, *) V, (1.0d0 - exp(-dt*dble(istep)))
    open(11, file = 'chk_timestep_precision.dat', position='append')
    write(11, *) dt, error
    close(11)
end program main
