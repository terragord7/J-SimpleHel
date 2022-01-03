# DESCRIPTION
# Simple model of a helicopter using static aerodynamic model for rotor,
# tail rotor, horzontal and vertical tail, and fuselage.
#
# INPUT
# - state: state vector
# - ctrls: stick control inputs
#
# OUTPUT
# - dstate: state derivative
# - FM: total forces and moments
#
# ------------------------------------------------------------------------------

function SimpleHel!(state,sctrls)

      # unpack states
      u=state[1]
      v=state[2]
      w=state[3]
      p=state[4]
      q=state[5]
      r=state[6]
      phi=state[7]
      theta=state[8]
      psi=state[9]
      x=state[10]
      y=state[11]
      z=state[12]
      # convert to stick inputs to swashplate and TR collective inputs
      ctrls = control_mixing(sctrls)
      # unpack controls
      theta1c=ctrls[1]*pi/180
      theta1s=ctrls[2]*pi/180
      theta0=ctrls[3]*pi/180
      theta0tr=ctrls[4]*pi/180

      # ----------------------------- MAIN ROTOR -------------------------------

      # hub velocities
      Vh=[u;v;w]+cross([p;q;r],rMR)
      uh=Vh[1]
      vh=Vh[2]
      wh=Vh[3]
      # advance ratios
      Vinplane=sqrt(uh^2+vh^2)
      mu=Vinplane/(OMEGA*RADIUS)
      muz=wh/(OMEGA*RADIUS)
      mu2=mu*mu
      # get wind axis transformation
      psiw=atan2(vh,uh)
      cpsiw=cos(psiw)
      spsiw=sin(psiw)
      Th2w=[ cpsiw spsiw
            -spsiw cpsiw]
      Tw2h=Th2w'
      # non-dimensional angular rates in wind axes
      temp=Th2w*[p; q]
      phw=temp[1]/OMEGA
      qhw=temp[2]/OMEGA
      # cyclic pitch translated to wind axes
      temp=Th2w*[theta1s;theta1c]
      theta1sw=temp[1]
      theta1cw=temp[2]
      # CT and lambda0 iteration, solves for inflow and thrust coefficient
      # first guess
      lambda0=0.05
      # set tolerance greater than error
      delta_lambda=1
      iter=0
      function inflow_MR()
            while (abs(delta_lambda)>1e-10 && iter<200)
                  CT=0.5*A0*SOL*((1/3+0.5*mu2)*theta0+0.5*mu*(theta1sw+0.5*phw)+
                     0.5*(muz-lambda0)+0.25*(1+mu2)*TWIST)
                  lamtot2=mu2+(lambda0-muz)^2
                  delta_lambda=-(2.0*lambda0*sqrt(lamtot2)-CT)*lamtot2/(2.0*lamtot2^1.5+0.25*
                               A0*SOL*lamtot2-CT*(muz-lambda0))
                  lambda0=lambda0+0.5*delta_lambda
            iter=iter+1
            end
            return CT
      end
      CT=inflow_MR()
      # thrust
      T=CT*RHO*OMEGA^2*pi*RADIUS^4
      # uniform inflow assumption in this model
      lambda1cw=0.0
      lambda1sw=0.0
      # quasi-Steady flapping in wind axes (Padfield p. 107)
      Abt = [-8/3*mu*(1+0.5*mu2)    -2*mu*(1+0.5*mu2)        -(1+2*mu^2)      0.0
             -GAMMA/6*mu*(1+0.5*mu2) -2*GAMMA*mu/15*(1+mu2/3) -2*GAMMA*mu2/9.0 (1-0.5*mu2*mu2)]
      Abl = [-2*mu*(1+0.5*mu2)          (1+0.5*mu^2)                    0.0
             -2*GAMMA*mu/9.0*(1-0.5*mu2)  GAMMA/9.0*mu^2+0.5*GAMMA/9.0*mu2  -(1.0-0.5*mu2) ]
      Abo = [-(1+0.5*mu^2)                          16/GAMMA*(1+0.5*mu2)
             16/GAMMA*(1.0-0.5*mu2)+GAMMA/9.0*mu^2  (1.0-0.5*mu2)]
      #
      temp=Abt*[theta0; TWIST; theta1sw; theta1cw]+Abl*
               [muz-lambda0; lambda1sw; lambda1cw]+Abo*[phw; qhw]
      beta1cw=temp[1]
      beta1sw=temp[2]
      # convert back to hub system
      temp=Tw2h*[beta1sw; beta1cw]
      beta1s=temp[1]
      beta1c=temp[2]
      # torque
      delta=DELTA0+DELTA2*CT*CT
      CQ=-(muz-lambda0)*CT+0.125*delta*SOL*(1+(7/3)*mu^2)
      Q=CQ*RHO*OMEGA^2*pi*RADIUS^5
      # translate rotor forces and moments to CG
      Xr=T*beta1c
      Yr=-T*beta1s
      Zr=-T
      MvecMR=[0.0; 0.0; Q]+cross(rMR, [Xr; Yr; Zr])
      Lr=MvecMR[1]
      Mr=MvecMR[2]
      Nr=MvecMR[3]

      # ------------------------------ FUSELAGE --------------------------------

      # velocity vector
      Vf=[u; v; w]
      # drag force vector
      Ff=-0.5*RHO*sqrt(Vf'*Vf)*Vf*CDSF
      # transform to body forces
      Xf=Ff[1]
      Yf=Ff[2]
      Zf=Ff[3]
      MvecF=cross(rFS,[Xf;Yf;Zf])
      Lf=MvecF[1]
      Mf=MvecF[2]
      Nf=MvecF[3]

      # ----------------------- HORIZONTAL STABILIZER --------------------------

      # local angle of attack
      Vht=[u; v; w]+cross([p; q; r], rHT)
      uht=Vht[1]
      wht=Vht[3]
      alphaht=atan2(wht,uht)
      # lift force
      CLht=AHT*alphaht
      LIFTht=0.5*RHO*(Vht'*Vht)*SHT*CLht
      # transform to body forces at CG
      Xht=0.0
      Yht=0.0
      Zht=-LIFTht
      MvecHT=cross(rHT,[Xht;Yht;Zht])
      Lht=MvecHT[1]
      Mht=MvecHT[2]
      Nht=MvecHT[3]

      # ----------------------- VERTICAL STABILIZER --------------------------

      # local angle of attack
      Vvt=[u; v; w]+cross([p; q; r], rVT)
      uvt=Vvt[1]
      vvt=Vvt[2]
      alphavt=atan2(vvt,uvt)
      # lift force
      CLvt=AVT*alphavt
      LIFTvt=0.5*RHO*(Vvt'*Vvt)*SVT*CLvt
      # transform to body forces at CG
      Xvt=0.0
      Yvt=LIFTvt
      Zvt=0.0
      MvecVT=cross(rVT,[Xvt;Yvt;Zvt])
      Lvt=MvecVT[1]
      Mvt=MvecVT[2]
      Nvt=MvecVT[3]

      # ----------------------------- TAIL ROTOR -------------------------------

      # local velocities at tail rotor hub in body CS
      Vtr=[u;v;w]+cross([p;q;r],rTR)
      # transfrom to TR coordinate system
      Ttrg=[1.0  0.0 0.0
            0.0  0.0 1.0
            0.0 -1.0 0.0]
      Vtr_tr=Ttrg*Vtr
      utr=Vtr_tr[1]
      vtr=Vtr_tr[2]
      wtr=Vtr_tr[3]
      mutr=sqrt(utr^2+vtr^2)/(OMEGTR*RTR)
      muztr=wtr/(OMEGTR*RTR)
      # CT and lambda0 iteration for tail rotor
      # first guess
      lambda0tr=0.05
      #set tolerance greater than error
      delta_lambda=1
      iter=0
      function inflow_TR()
            while(abs(delta_lambda)>1e-10 && iter<200)
                  CTtr=0.5*A0TR*SOLTR*((1/3+0.5*mutr^2)*theta0tr+0.5*(muztr-lambda0tr)+
                       0.25*(1+mutr^2)*TWISTTR)
                  lamtot2=mutr^2+(lambda0tr-muztr)^2
                  delta_lambda=-(2.0*lambda0tr*sqrt(lamtot2)-CTtr)*lamtot2/(2.0*lamtot2^1.5+
                               0.25*A0TR*SOLTR*lamtot2-CTtr*(muztr-lambda0tr))
                  lambda0tr=lambda0tr+0.5*delta_lambda
                  iter=iter+1
            end
            return CTtr
      end
      CTtr=inflow_TR()
      # tail rotor thrust
      Ttr=CTtr*RHO*OMEGTR^2*pi*RTR^4
      # force vector in body coordinates
      Ftr=Ttrg'*[0.0;0.0;-Ttr]
      Xtr=Ftr[1]
      Ytr=Ftr[2]
      Ztr=Ftr[3]
      # resolve force and moment about CG
      MvecTR=cross(rTR,Ftr)
      Ltr=MvecTR[1]
      Mtr=MvecTR[2]
      Ntr=MvecTR[3]

      # ------------------------ EQUATIONS OF MOTION ---------------------------

      # sum forces and moments
      X=Xr+Xf+Xtr+Xht+Xvt
      Y=Yr+Yf+Ytr+Yht+Yvt
      Z=Zr+Zf+Ztr+Zht+Zvt
      L=Lr+Lf+Ltr+Lht+Lvt
      M=Mr+Mf+Mtr+Mht+Mvt
      N=Nr+Nf+Ntr+Nht+Nvt
      # total forces and moments
      FM=zeros(6,1)
      FM=[X;Y;Z;L;M;N]

      # EQUATIONS OF MOTION
      # trigonometric functions
      cphi=cos(phi)
      sphi=sin(phi)
      cthe=cos(theta)
      sthe=sin(theta)
      cpsi=cos(psi)
      spsi=sin(psi)
      # declare state derivative vector
      dstate=zeros(12,1)
      # equations of motion
      # translational dynamics
      dstate[1]=X/mass-g*sthe-q*w+r*v
      dstate[2]=Y/mass+g*cthe*sphi-r*u+p*w
      dstate[3]=Z/mass+g*cthe*cphi-p*v+q*u
      # rotational dynamics
      gam=Ix*Iz-Ixz^2
      dstate[4]=(Iz*L+Ixz*N+Ixz*(Ix-Iy+Iz)*p*q-(Iz^2-Iy*Iz+Ixz^2)*q*r)/gam
      dstate[5]=(M+(Iz-Ix)*p*r-Ixz*(p^2-r^2))/Iy
      dstate[6]=(Ix*N+Ixz*L-Ixz*(Ix-Iy+Iz)*q*r+(Ix^2-Ix*Iy+Ixz^2)*p*q)/gam
      # rotational kinematics
      dstate[7]=p+q*sphi*sthe/cthe+r*cphi*sthe/cthe
      dstate[8]=q*cphi-r*sphi
      dstate[9]=q*sphi/cthe+r*cphi/cthe
      # position
      dstate[10]=u*cthe*cpsi+v*(sphi*sthe*cpsi-cphi*spsi)+w*(cphi*sthe*cpsi+sphi*spsi)
      dstate[11]=u*cthe*spsi+v*(sphi*sthe*spsi+cphi*cpsi)+w*(cphi*sthe*spsi-sphi*cpsi)
      dstate[12]=-u*sthe+v*sphi*cthe+w*cphi*cthe
      #
      return dstate, FM
end
