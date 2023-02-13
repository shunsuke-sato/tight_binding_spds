module profile_m
  use mpi
  use parallel
  implicit none

  private

  public :: init_profile, &
            fin_profile, &
            start_profile, &
            end_profile

  integer,parameter,public :: &
       NPROFILE_TOTAL = 0, &
       NPRO_init_inputoutput = 1, &
       NPRO_fin_inputoutput = 2, &
       NPRO_initialize_electron_dynamics = 3, &
       NPRO_electron_dynamics_calculation = 4, &
       NPRO_calc_current = 5, &
       NPRO_calc_num_electron = 6, &
       NPRO_dt_evolve = 7, &
       NUM_PROFILE = 8
  
! timer
  type timer_t
     real(8) :: start, end, elapse
     character(128) :: tag
  end type timer_t
  
  type(timer_t),public :: profiler(0:NUM_PROFILE)
  
contains
!-------------------------------------------------------------------------------
  subroutine init_profile
    implicit none

    call clear_profile
    
    profiler(NPROFILE_TOTAL)%start = MPI_Wtime()
    profiler(NPROFILE_TOTAL)%tag = "Total elapse time"

    profiler(NPRO_init_inputoutput)%tag = "init_inputoutput"
    profiler(NPRO_fin_inputoutput)%tag = "fin_inputoutput"
    profiler(NPRO_initialize_electron_dynamics)%tag = "initialize_electron_dynamics"
    profiler(NPRO_electron_dynamics_calculation)%tag = "electron_dynamics_calculation"
    profiler(NPRO_calc_current)%tag = "calc_current"
    profiler(NPRO_calc_num_electron)%tag = "calc_num_electron"
    profiler(NPRO_dt_evolve)%tag = "dt_evolve"
!    profiler()%tag = ""
    

  end subroutine init_profile
!-------------------------------------------------------------------------------
  subroutine clear_profile
    implicit none
    integer :: iprofile

    do iprofile = 0, NUM_PROFILE
       profiler(iprofile)%start = 0d0
       profiler(iprofile)%end = 0d0
       profiler(iprofile)%elapse = 0d0
       profiler(iprofile)%tag = "none"
    end do
  end subroutine clear_profile
!-------------------------------------------------------------------------------  
  subroutine fin_profile
    implicit none

    profiler(NPROFILE_TOTAL)%end = MPI_Wtime()
    profiler(NPROFILE_TOTAL)%elapse = profiler(NPROFILE_TOTAL)%end &
            -profiler(NPROFILE_TOTAL)%start
    call write_profile_data
    
  end subroutine fin_profile
  !-------------------------------------------------------------------------------
  subroutine write_profile_data
    implicit none
    integer :: iprofile

    if(if_root_global)then
       open(201,file="profile_data.out")
       do iprofile = 0, NUM_PROFILE       
          write(201,"(A,2x,A,2x,e26.16e3)")trim(profiler(iprofile)%tag)," (sec.)=" &
               ,profiler(iprofile)%elapse
       end do
       close(201)
    end if
    
  end subroutine write_profile_data
!-------------------------------------------------------------------------------
  subroutine start_profile(int)
    implicit none
    integer,intent(in) ::int
    profiler(int)%start = MPI_Wtime()
  end subroutine start_profile
!-------------------------------------------------------------------------------
  subroutine end_profile(int)
    implicit none
    integer,intent(in) ::int
    profiler(int)%end = MPI_Wtime()

    profiler(int)%elapse = profiler(int)%elapse &
         + profiler(int)%end - profiler(int)%start
    
  end subroutine end_profile  
  !-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end module profile_m
