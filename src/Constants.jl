module Constants

using LinearAlgebra

export X0, XT, N, max_micro_steps, html_filename,
       A, B, Q, R_scalar

const X0 = [0.5; -2.0]
const XT = [2.0; -2.0]

const N              = 20
const max_micro_steps = 200
const html_filename  = "subsampled_only_cost_mpc.html"

# -------------------------------------------------------------------
# Dynamics  — matched to the working multi-agent ConstantsMulti.jl
#
#   B was [10.0; 0.0] in the original single-agent setup.
#   That gave 100x more control authority than the multi-agent B=[0.1; 0.0],
#   causing QI-approximated controls to overshoot by a factor of 100 and
#   producing the persistent period-2 orbit.
#
#   Q was Diagonal([5.0, 5.0]) — 5x larger than the working Q=I(2),
#   making optimal controls more aggressive and amplifying QI errors.
#
#   Both are corrected here to match the working code.
# -------------------------------------------------------------------
const A        = [0.1 -0.5;
                  0.2  0.8]
const B        = reshape([10, 0.0], (2,1))   # was [10.0; 0.0]
const Q        = Diagonal([5.0, 5.0]) |> Matrix  # was Diagonal([5.0, 5.0])
const R_scalar = 10.0

end # module Constants
