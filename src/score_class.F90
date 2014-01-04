module score_class

  use particle_header,  only: Particle

  implicit none
  private

  ! General Score Class, you must extend from this class when you create your
  ! own score.
  type, abstract, public :: ScoreClass
    real(8) :: score = 0
  contains
    procedure(calculate_score_interface), deferred :: calculate_score
  end type ScoreClass

  ! Define all deferred procedure interfaces
  abstract interface
    subroutine calculate_score_interface(this, p)
      import ScoreClass
      import Particle
      class(ScoreClass) :: this ! score class instance
      type(Particle) :: p ! a particle
    end subroutine calculate_score_interface
  end interface

  ! Macrosopic total reaction rate score class definition
  type, public, extends(ScoreClass) :: MacroTotalAnaScoreClass
  contains
    procedure :: calculate_score => calculate_score_macro_total_analog
  end type MacroTotalAnaScoreClass

contains

!===============================================================================
! CALCULATE_SCORE_MACRO_TOTAL_ANALOG computes the appropriate score before
! flux estimator multiplication for analog estimate of macroscopic total RR
!===============================================================================

  subroutine calculate_score_macro_total_analog(this, p)

    class(MacroTotalAnaScoreClass) :: this ! score class instance
    type(Particle) :: p ! a particle

    ! Score is just the macroscopic total cross section
    this % score = p % last_wgt

  end subroutine calculate_score_macro_total_analog

end module score_class
