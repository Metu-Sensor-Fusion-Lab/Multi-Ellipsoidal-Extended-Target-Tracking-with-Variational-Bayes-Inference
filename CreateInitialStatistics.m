function initstats=CreateInitialStatistics(Y,algorithm)
%Using Multi Ellipse VB

ny=size(Y,1);
Nmeas=size(Y,2);

s=algorithm.VB.s;
L=algorithm.VB.L;

meanY=mean(Y,2);
if Nmeas<=ny
    covY=1e-2*eye(ny);
else
    covY=cov(Y')/s;
end

mu=zeros(ny,L);%initial means
sigma=repmat(covY,1,1,L);%initial covariances
piVB=ones(1,L)/L;%equal initial probabilities

[V, D]=eig(covY);
[maxEig, maxind]=max(sqrt(diag(D)));
maxVec=V(:,maxind);
h=2*maxEig/L;%distance between initial means
if mod(L,2)
    mu(:,1)=meanY;
    ind=1;
    for i=1:(L-1)/2
        ind=ind+1;
        mu(:,ind)=meanY+i*h*maxVec;
        ind=ind+1;
        mu(:,ind)=meanY-i*h*maxVec;
    end
else
    ind=0;
    for i=1:L/2
        ind=ind+1;
        mu(:,ind)=meanY+(2*i-1)*h*maxVec/2;
        ind=ind+1;
        mu(:,ind)=meanY-(2*i-1)*h*maxVec/2;
    end
end
initstats.a=algorithm.Prior.a;
initstats.x = [mu(:,1); zeros(ny,1)];
initstats.P = zeros(4);
for ell = 2:L
    initstats.x = [initstats.x; mu(:,ell)-mu(:,1)];
    initstats.P = blkdiag(initstats.P,zeros(2));
end
initstats.v=algorithm.Prior.v;
for ell=1:L
    initstats.V(:,:,ell)=(algorithm.Prior.v(ell)-2*ny-2)*(sigma(:,:,ell)/L^2);
end




