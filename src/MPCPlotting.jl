module MPCPlotting

using PlotlyJS
using PlotlyJS: attr

export save_mpc_plot

function save_mpc_plot(
    history_z,
    history_z_times,
    history_u,
    history_errors,
    XT,
    max_micro_steps,
    ν,
    html_filename
)
    n_z = length(history_z)
    n_u = length(history_u)
    n_e = length(history_errors)

    t_z = history_z_times[1:n_z]
    t_u = collect(1:n_u)

    p_sub = [s[1] for s in history_z]
    v_sub = [s[2] for s in history_z]
    # Ensure no space between . and (
    log_err = log10.(history_errors[1:n_e] .+ 1e-16)

    tr_p_sub = scatter(
        x = t_z, y = p_sub,
        mode = "markers+lines",
        name = "position (subsampled)",
        marker = attr(size=8, symbol="diamond"),
        yaxis = "y1"
    )

    tr_p_target = scatter(
        x = [0, max_micro_steps], y = [XT[1], XT[1]],
        mode = "lines", name = "p_target",
        line = attr(dash="dash", width=2), yaxis = "y1"
    )

    tr_v_sub = scatter(
        x = t_z, y = v_sub,
        mode = "markers+lines",
        name = "velocity (subsampled)",
        marker = attr(size=8, symbol="diamond"),
        yaxis = "y2"
    )

    tr_v_target = scatter(
        x = [0, max_micro_steps], y = [XT[2], XT[2]],
        mode = "lines", name = "v_target",
        line = attr(dash="dash", width=2), yaxis = "y2"
    )

    tr_u = scatter(
        x = t_u, y = history_u,
        mode = "lines", name = "u (micro)",
        line = attr(width=1.5), yaxis = "y3"
    )

    tr_logerr = scatter(
        x = t_z[1:n_e], y = log_err,
        mode = "lines+markers", name = "log10(error)",
        line = attr(width=2), marker = attr(size=6),
        yaxis = "y4"
    )

    layout = Layout(
        title = "QI-MPC: Subsampled States (nu=$ν)",
        xaxis = attr(title="micro-steps", range=[0, max_micro_steps]),
        yaxis = attr(domain=[0.76, 1.00], title="position"),
        yaxis2 = attr(domain=[0.52, 0.74], title="velocity"),
        yaxis3 = attr(domain=[0.26, 0.50], title="control u"),
        yaxis4 = attr(domain=[0.00, 0.24], title="log10 error"),
        margin = attr(l=70, r=140, t=80, b=60)
    )

    fig = Plot([tr_p_sub, tr_p_target, tr_v_sub, tr_v_target, tr_u, tr_logerr], layout)
    savefig(fig, html_filename)
end

end