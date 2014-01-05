module tally_class

  use constants,     only: N_FILTER_TYPES
  use filter_class,  only: FilterClassPointer
  use score_class,   only: ScoreClassPointer
  use tally_header,  only: TallyResult

  implicit none
  private

  ! General Tally Class, you must extend from this class when you create your
  ! own tally.
  type, abstract, public :: TallyClass

    ! Basic data
    integer :: id                   ! user-defined identifier
    character(len=52) :: label = "" ! user-defined label
    real(8) :: volume               ! volume of region

    ! Filter information
    integer :: n_filters ! number of filters
    type(FilterClassPointer), allocatable :: filters(:)

    ! The stride attribute is used for determining the index in the results
    ! array for a matching_bin combination. Because multiple dimensions are
    ! mapped onto one dimension in the results array, the stride attribute gives
    ! the stride for a given filter type within the results array
    integer, allocatable :: stride(:)

    ! This array provides a way to lookup what index in the filters array a
    ! certain filter is. For example, if find_filter(FILTER_CELL) > 0, then the
    ! value is the index in filters.
    integer :: find_filter(N_FILTER_TYPES) = 0

    ! Values to score, e.g. flux, absorption, etc.
    ! scat_order is the scattering order for each score.
    ! It is to be 0 if the scattering order is 0, or if the score is not a
    ! scattering response
    integer :: n_scores = 0
    type(ScoreClassPointer), allocatable :: scores(:)

    ! Results for each bin
    integer :: total_nuclide_bins = 1
    integer :: total_filter_bins
    integer :: total_score_bins
    type(TallyResult), allocatable :: results(:,:,:)

  end type TallyClass

  ! Define all deferred procedure interfaces
  abstract interface
  end interface

contains

end module tally_class
