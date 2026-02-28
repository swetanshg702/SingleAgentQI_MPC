# run_mpc_qi.jl (Full Corrected Version)
include("src/Constants.jl")
include("src/SystemDynamics.jl")
include("src/BlockQP.jl")
include("src/QI.jl")
include("src/MPCPlotting.jl")
include("src/QIGridPrecompute.jl")

using .Constants, .SystemDynamics, .BlockQP, .QI, .MPCPlotting, .QIGridPrecompute
using LinearAlgebra, Printf

function main()
    # 1. Initialization
    ν = compute_reachability_index(A, B, XT; max_ν = 10)
    Nblk = max(1, Int(N ÷ ν))
    
    # 2. Fresh Grid Computation (using [-3, 3] from QI.jl)
    cost_grid, control_grid = precompute_cost_control_grids(ν, Nblk)

    # 3. Simulation Setup
    z = copy(X0)
    history_z, history_z_times = [copy(z)], [0]
    history_u, history_errors = Float64[], Float64[]
    micro_t = 0
    J_accumulated = 0.0

    println("Starting Simulation...")

    # 4. MPC Loop
    while micro_t < max_micro_steps
        # Attempt QI for the first control block
        u_block = quasi_interpolate_control_vector(
            z[1], z[2], x_fine, y_fine, control_grid, D, h_fine, rho
        )

        # Reference values for cost calculation
        Uref = compute_Uref(A, B, XT, ν)

        # Fallback Mechanism
        if u_block === nothing
            println("Using fallback MPC computation at micro_t = $micro_t for state $z")
            try
                sol = solve_block_QP(z, XT, Nblk; A=A, B=B, Q=Q, R=R_scalar, ν=ν)
                u_block = sol.U_sol[:, 1]
                Uref = sol.Uref
            catch e
                @warn "Both QI and Fallback failed at state $z. Stopping."; break
            end
        end

        # Apply the chosen block (ν micro-steps)
        for j = 1:ν
            micro_t >= max_micro_steps && break

            u = u_block[j]
            
            # Integration and Costs
            z_next = apply_micro_step(A, B, z, u)
            
            control_cost = R_scalar * (u - Uref[j])^2
            state_cost = dot(z_next - XT, Q * (z_next - XT))
            J_accumulated += (state_cost + control_cost)
            
            z = z_next
            push!(history_u, u)
            micro_t += 1
        end

        # Record Subsampled State
        push!(history_z, copy(z))
        push!(history_z_times, micro_t)
        push!(history_errors, norm(z - XT))
    end

    # 5. Output and Plotting
    output_file = "subsampled_qi_mpc.html"
    save_mpc_plot(
        history_z, 
        history_z_times, 
        history_u, 
        history_errors, 
        XT, 
        max_micro_steps, 
        ν, 
        output_file
    )
    
    println("\n--- Simulation Summary ---")
    println("Final State: $z")
    println("Final Error: $(norm(z - XT))")
    println("Accumulated Cost: $J_accumulated")
    println("Plot saved to: $(joinpath(pwd(), output_file))")
end

main()