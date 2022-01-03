# DESCRIPTION
# Trims the aircraft dynamics at an arbitrary flight condition.
#
# INPUT
# - acfun: name of system dynamics function
# - x0: initial guess for trim state vector
# - u0: initial guess for trim control vector
# - targ_des: vector of trim targets
#
# OUTPUT
# - x0trim: trim state vector
# - u0trim: trim control vector 
# - itrim: number of trim iterations
#
# ------------------------------------------------------------------------------

function trimmer(acfun,x0,u0,targ_des)
    #
    #x0trim=zeros(NSTATES)
    #u0trim=zeros(NCTRLS)
    #x0trim[:]=x0[:]
    #u0trim[:]=u0[:]
    x0trim=x0
    u0trim=u0
    #
    it=0
    conv=0
    trim_tol=5e-5
    itmax=80
    #
    @printf("trimming the aircraft\n")
    while ((it<itmax) && (conv==0))
        it=it+1
        # azimuthal average of xdot
        xdot0, y0=acfun(x0trim,u0trim)
        targvec=[xdot0;y0]
        targvec=targvec[TRIMTARG]
        targ_err=targvec-targ_des
        #
        @printf("%d          %f\n",it,norm(targ_err))
        #
        conv=1
        for k=1:NSTATES
            if (abs(targ_err[k])>TOL[k])
                conv=0
            end
        end
        #
        if (conv==0)
            #@show x0trim
            #@show u0trim
            A, B, C, D =linearize(acfun,x0trim,u0trim)
            #@show x0trim
            #@show u0trim
            #return
            # jacobian matrix
            Jac=[A B
                 C D]
            Jac2 = Jac[TRIMTARG,TRIMVARS]
            trimvec=[x0trim;u0trim]
            trimvec[TRIMVARS]=trimvec[TRIMVARS]-0.5*pinv(Jac2)*targ_err
            #trimvec[TRIMVARS]=trimvec[TRIMVARS]-Jac2\targ_err
            x0trim=trimvec[1:NSTATES]
            u0trim=trimvec[NSTATES+1:NSTATES+NCTRLS]
        end
    end
    #
    if (conv==0)
        @printf("warning: trim not achieved\n\n")
        itrim=0
    else
        @printf("successful trim\n\n")
        itrim=1
    end
    #
    return x0trim, u0trim, itrim
end
