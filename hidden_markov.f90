module hidden_markov_mod
use kind_mod, only: dp
contains
subroutine hidden_markov_normal_sim(trans_prob,xmean,xsd,xx)
real(kind=dp), intent(in)  :: trans_prob(:,:)
real(kind=dp), intent(in)  :: xmean(:)
real(kind=dp), intent(in)  :: xsd(:)
real(kind=dp), intent(out) :: xx(:)
integer                    :: nstates
nstates = size(trans_prob,1)
end subroutine hidden_markov_normal_sim
end module hidden_markov_mod