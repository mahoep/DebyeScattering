module angular_functions_mod
  implicit none
contains
  subroutine angularFunctions(theta, Nmx, pi, tau)
    use iso_fortran_env, wp => real128  ! quad  
    ! Input parameters
    real(16), intent(in) :: theta(:)
    integer, intent(in) :: Nmx
    
    ! Output arrays
    real(16), allocatable, intent(out) :: pi(:,:), tau(:,:)
    
    ! Local variables
    integer :: i, n
    real(16) :: en
    integer :: theta_size
    real(16), allocatable :: cos_theta(:)
    
    ! Get the size of the theta array
    theta_size = size(theta)
    
    ! Allocate the pi and tau arrays
    allocate(pi(0:Nmx, theta_size))
    allocate(tau(0:Nmx, theta_size))
    allocate(cos_theta(theta_size))
    
    ! Initialize the arrays with zeros
    pi = 0.0_wp
    tau = 0.0_wp
    
    ! Compute cos(theta) and store it in cos_theta for reuse
    cos_theta = cos(theta)
    
    ! Set initial values for pi and tau
    pi(0, :) = 0.0_wp
    pi(1, :) = 1.0_wp
    
    tau(0, :) = 0.0_wp
    tau(1, :) = cos_theta * pi(1, :) - 2.0_wp * pi(0, :)
    
    ! Loop to calculate the angular functions for n >= 2
    do n = 2, Nmx
        en = real(n)
        pi(n, :) = (2.0_wp*en-1.0_wp) / (en-1.0_wp) * cos_theta * pi(n-1, :) - n / (en-1.0_wp) * pi(n-2, :)
        tau(n, :) = en * cos_theta * pi(n, :) - (en+1.0_wp) * pi(n-1, :)
    end do
  
  end subroutine angularFunctions
end module angular_functions_mod