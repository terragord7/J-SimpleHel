function simulate(acfun,finp,Tsim,dt,x0,u0)
    # set up simulation
    # vector of times
    time=collect(0:dt:Tsim)
    # initialize state time history
    state=zeros(NSTATES,length(time))
    # set initial condition
    state[:,1]=x0
    # simulate open-loop system
    for i=1:length(time)-1
        # integrate
        sol=rk4(acfun,finp,time[i],dt,state[:,i],u0)
        state[:,i+1]=sol
    end
    # plot simulation
    #gr(size=(1000,1000))
    #plt3=plot(time,state[4,:],linewidth=2,xaxis="Time [s]",yaxis="Pitch rate, q [rad/s]",label="A",reuse = false)
    #display(plt3)
    #
    return state, time
end
