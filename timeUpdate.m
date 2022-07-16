function track=timeUpdate(track,model,algorithm)

ny=model.measurement.n;
F=model.state.F;
Q=model.state.Q;
lambda=algorithm.VB.lambda;

track.x=F*track.x;
track.P=F*track.P*F'+Q;
track.v=max(lambda*track.v,2*ny+1);%avoids v falling below the minimum required for IW density
track.V=lambda*track.V;
track.a=lambda*track.a;
