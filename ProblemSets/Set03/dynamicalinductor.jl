using DynamicalSystems, Plots

@inline @inbounds function eom_inductor(u, p, t)
    k = p[1]; B = p[2]
    du1 = u[2]
    du2 = -k*u[2] - u[1]^3 +  B*cos(t)
    return SVector{2}(du1, du2)
end


@inline @inbounds function inductor_jac(u, p, t)
    k, B = p
    J = @SMatrix [0  1;
    -3*u[1]^2  -k]
    return J
end

function poinc(ttotal,tstep,B)
    ds = ContinuousDynamicalSystem(eom_inductor,zeros(2),[0.1,B],inductor_jac,t0 = 0.0)
    a = trajectory(ds, ttotal, dt = tstep,Ttr=100.0) # every period T = 2π/ω
    sc = scatter(a[:,1],a[:,2], marker = (:auto, 1, 0.0, :blue,0,0.2,:black,:solid),dpi=800,size=(600,400))
    #sc = scatter(a[:,1],a[:,2], line = (:line, :solid, nothing, 0.6, 1, :red))
    return sc
end



poinc(10.0,10.0,9.80)
#poinc(0.001,9.80)
