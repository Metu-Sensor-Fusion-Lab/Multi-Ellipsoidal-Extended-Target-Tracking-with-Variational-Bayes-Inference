L=2; %number of ellipsoidal extent components

T=0.1;
model.state.T=T;

ny=2;
model.measurement.n=ny;%measurement dimension
Zny=zeros(ny);
Iny=eye(ny);
model.measurement.R=0.001*Iny;%measurement noise covariance
model.measurement.s=1/4;%s pararameter
model.measurement.gamma=30;%Average Number of Measurements of Extended Targets

tempF=eye(ny);%F for non-center components
F=[Iny T*Iny; Zny Iny];
sigmaQ=1;
tempQ=eye(ny);%Q for non-center components
Q=sigmaQ^2*[T^3/3*Iny T^2/2*Iny; T^2/2*Iny  T*Iny];
model.measurement.H=cell(1,L);

model.state.F=F;
model.state.Q=Q;
for ell=1:L-1
    model.state.F=blkdiag(model.state.F,tempF);
    model.state.Q=blkdiag(model.state.Q,tempQ);
end
model.measurement.H=repmat([Iny repmat(Zny,1,L)],1,1,L);
for l = 2:L
    model.measurement.H(1,2*l+1,l) = 1;
    model.measurement.H(2,2*l+2,l) = 1;
end

P0prior=blkdiag(10000*Iny,10000*Iny); %Prior covariance
vprior=10;%Prior v value
Vprior=10*Iny;%Prior V value
aprior=1;%Prior a value
Pe=100^2*eye(2); %Prior covariance of  mixture component positions from the center

algorithm.Prior.a=aprior*ones(1,L);
algorithm.Prior.x=zeros(4,1); %prior target state is zero
algorithm.Prior.P=P0prior;
algorithm.Prior.v=vprior*ones(1,L);
algorithm.Prior.V=repmat(Vprior,1,1,L);
for ell = 2:L
    algorithm.Prior.x = [algorithm.Prior.x; zeros(ny,1)];
    algorithm.Prior.P = blkdiag(algorithm.Prior.P,Pe);
end

%Algorithm Parameters
algorithm.VB.lambda=0.9; %time update forgetting factor
algorithm.VB.maxNiter=30; %Maximum number of iterations in VB
algorithm.VB.convergenceThreshold=1e-3; %convergence threshold in VB
algorithm.VB.L=L;%number of ellipsoid components
algorithm.VB.s=1/4;%s value for VB
