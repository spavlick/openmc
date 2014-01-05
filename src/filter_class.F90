module filter_class

  use global,          only: meshes
  use mesh,            only: get_mesh_bin
  use mesh_header,     only: StructuredMesh
  use particle_header, only: Particle

  implicit none
  private

  ! General Filter Class, you must extend from this class when you create your
  ! own filter.
  type, abstract, public :: FilterClass
    integer :: n_bins = 0
  contains
    procedure(initialize_interface), deferred :: initialize
    procedure(finalize_interface), deferred :: finalize
    procedure(get_filter_bin_interface), deferred :: get_filter_bin
!   procedure :: get_next_bin => get_next_bin_fun
  end type FilterClass

  ! Define all deferred procedure interfaces
  abstract interface
    subroutine initialize_interface(this)
      import FilterClass
      class(FilterClass) :: this
    end subroutine initialize_interface
    subroutine finalize_interface(this)
      import FilterClass
      class(FilterClass) :: this
    end subroutine finalize_interface    
    function get_filter_bin_interface(this, p) result(bin)
      import FilterClass
      import Particle
      class(FilterClass) :: this
      type(Particle), intent(in) :: p
      integer :: bin
    end function get_filter_bin_interface 
  end interface

  ! Mesh filter class definition
  type, public, extends(FilterClass) :: MeshFilterClass
    integer, allocatable :: bins(:)
  contains
    procedure :: initialize => initialize_filter_mesh
    procedure :: finalize => finalize_filter_mesh
    procedure :: get_filter_bin => get_filter_bin_mesh
  end type MeshFilterClass

  ! Set up pointer type for this class
  type, public :: FilterClassPointer
    class(FilterClass), pointer :: p
  end type FilterClassPointer

contains

!===============================================================================
! INITIALIZE_FILTER_MESH initializes a mesh filter
!===============================================================================

  subroutine initialize_filter_mesh(this)

    class(MeshFilterClass) :: this ! mesh filter instance

    ! Allocate bins array
    if(.not.allocated(this % bins)) allocate(this % bins(this % n_bins))
    
  end subroutine initialize_filter_mesh

!===============================================================================
! FINALIZE_FILTER_MESH finalizes a mesh filter
!===============================================================================

  subroutine finalize_filter_mesh(this)

    class(MeshFilterClass) :: this ! mesh filter instance

    ! Free all memory
    if (allocated(this % bins)) deallocate(this % bins)

  end subroutine finalize_filter_mesh

!===============================================================================
! GET_FILTER_BIN_MESH finds the mesh bin where the particle is located and saves
! it to matching_bins
!===============================================================================

  function get_filter_bin_mesh(this, p) result(bin)

    class(MeshFilterClass) :: this ! mesh filter instance
    type(Particle), intent(in) :: p ! particle to tally
    integer :: bin ! integer bin location for matching bins array

    type(StructuredMesh), pointer :: m => null()

    ! Get the mesh
    m => meshes(this % bins(1))

    ! Get the mesh bin from the particles coordinates
    call get_mesh_bin(m, p % coord0 % xyz, bin)

    ! Set pointers to null
    m => null()

  end function get_filter_bin_mesh

end module filter_class
