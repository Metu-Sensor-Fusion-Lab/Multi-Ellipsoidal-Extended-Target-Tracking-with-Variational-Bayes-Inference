function [statsout, logpyy] = MultiEllipseVB(Y,stats,initstats,model,algorithm)
% Implementation of the measurement update of
% "Multi-Ellipsoidal Extended Target Tracking with Variational Bayes Inference"
% By Barkın Tuncer, Umut Orguner, Emre Özkan

% INPUTS
% Y: Measurements (2 x Ny double) (Ny is the number of measurements)

% stats: Predicted Statistics (structure)
% stats.x: Predicted kinematic mean (2*L+1 x 1 double)
% stats.P: Predicted kinematic covariance (2*L+1 x 2*L+1 double)
% stats.v: Predicted extent dof (1xL double)
% stats.V: Predicted extent scale matrices (ny x ny x L double)
% stats.a: Predicted Dirichlet parameters (1xL double)

% initstats: Initial Statistics (structure) (Optional. Set [] if not necessary. When set [] the VB iterations are started based on predicted statistics (i.e., stats))
% initstats.x: Initial kinematic mean (2*L+1 x 1 double)
% initstats.P: Initial kinematic covariance (2*L+1 x 2*L+1 double)
% initstats.v: Initial extent dof (1xL double)
% initstats.V: Initial extent scale matrices (ny x ny x L double)
% initstats.a: Initial Dirichlet parameters (1xL double)

% model: Measurement Model Parameters (structure)
% model.N: Measurement dimension ny (1 x 1 integer)
% model.R: Measurement noise Covariance (ny x ny double)
% model.H: Measurement matrices for extent components (ny x 2*L+1 x L double)
% model.s: Measurement scale parameter (1 x 1 double)

% algorithm: Algorithm Parameters (structure)
% algorithm.L=L: Number of extent ellipsoids (1 x 1 integer)
% algorithm.convergenceThreshold: Convergence threshold in VB (1 x 1 double)
% algorithm.maxNiter: Maximum number of iterations in VB (1 x 1 integer)

% OUTPUTS
% statsout: Updated Statistics (structure)
% statsout.x: Updated kinematic mean (2*L+1 x 1 double)
% statsout.P: Updated kinematic covariance (2*L+1 x 2*L+1 double)
% statsout.v: Updated extent dof (1xL double)
% statsout.V: Updated extent scale matrices (ny x ny x L double)
% statsout.a: Updated Dirichlet parameters (1xL double)

% logpyy: Predictive likelihood (1 x 1 double)

ny=model.n;
s=model.s;
R=model.R;
H=model.H;

L=algorithm.L;
maxNiter=algorithm.maxNiter;
convergenceThreshold=algorithm.convergenceThreshold;

statsout=stats;

%Predicted quantities
x = stats.x;
P = stats.P;
v = stats.v;
V = stats.V;
a = stats.a;

nx = length(x);

if isempty(initstats)
    %initialize VB iteration
    xkk_ = stats.x;
    %Pkk_= stats.P;
    Pkk_= zeros(nx);
    vkk_ = stats.v;
    Vkk_ = stats.V;
    akk_ = stats.a;
else
    xkk_ = initstats.x;
    Pkk_= initstats.P;
    vkk_ = initstats.v;
    Vkk_ = initstats.V;
    akk_ = initstats.a;
end
zkk_ = Y;
Nmeas = size(Y,2);
Pzkk_ = repmat(zeros(ny),1,1,Nmeas);

%Allocate
xkk = xkk_;
Pkk = Pkk_;
vkk = vkk_;
Vkk = Vkk_;
akk = akk_;
zkk = zkk_;
Pzkk = Pzkk_;

loggammakk = zeros(Nmeas,L);
invsX=zeros(ny,ny,L);
invsXHxkk_=zeros(ny,L);
Hxkk_=zeros(ny,L);
HPkk_H=zeros(ny,ny,L);
Wjl=zeros(ny,ny,Nmeas,L);

invR=inv(R);%R^-1
invP=inv(P);% inv(P_k|k-1)

iter=0;
convergenceStat=inf;

plotFlag=0; % For debugging purposes, make this 1, then the extent ellipsoids will be drawn at each VB iteration.
if plotFlag %Plot Extent Ellipsoids with Measurements (For debugging purposes)
    plotstats.x=xkk;
    plotstats.P=Pkk;
    plotstats.v=vkk;
    plotstats.V=Vkk;
    plotstats.a=akk;

    figure(547)
    plot(Y(1,:),Y(2,:),'.')
    hold on
    plot(zkk(1,:),zkk(2,:),'r.')

    [GM, ~, ~]=convertMultiEllipseTracktoGM(plotstats);
    for ell=1:L
        drawEllipse(GM.x(:,ell),GM.P(:,:,ell),1,gca,1.5,[0, 0.4470, 0.7410]);  
    end
    title(['Iteration No:', num2str(iter)])
    pause(0.5)
    cla
end

while iter<maxNiter &&  convergenceStat>convergenceThreshold
    iter=iter+1;
    oldxkk=xkk;
    for ell=1:L%Do these calculations only once in each iteration. Not as many as the number of measurements in the double for loop below
        invsX(:,:,ell)=(vkk_(ell) - ny - 1)*inv(s*Vkk_(:,:,ell)); %Find E[inv(sXk)] only once
        invsXHxkk_(:,ell) = invsX(:,:,ell) * H(:,:,ell) * xkk_;
        Hxkk_(:,ell)=H(:,:,ell)*xkk_;
        HPkk_H(:,:,ell)=H(:,:,ell) * Pkk_ * H(:,:,ell)';
        temploggammatildekk=psi(akk_(ell)) - psi(sum(akk_))- 0.5 * log(det(Vkk_(:,:,ell))) + 0.5*ny*log(2) + 0.5 *sum(psi((vkk_(ell) - ny - [1:ny])/2));
        for j=1:Nmeas
            Wjl(:,:,j,ell) = (zkk_(:,j) - Hxkk_(:,ell)) * (zkk_(:,j) - Hxkk_(:,ell))'+ Pzkk_(:,:,j)+ HPkk_H(:,:,ell);%22
            %loggammakk(j,ell) =   temploggammatildekk - 0.5*trace(invsX(:,:,ell)* Wjl(:,:,j,ell));%23
            loggammakk(j,ell) =   temploggammatildekk - 0.5*(invsX(1,:,ell)*Wjl(:,1,j,ell)+invsX(2,:,ell)*Wjl(:,2,j,ell));%23
        end
    end
    for j=1:Nmeas
        loggammakk(j,:)=loggammakk(j,:)-logSum(loggammakk(j,:));%Normalize
    end
    gammakk = exp(loggammakk);%26

    vkk = v + sum(gammakk,1);%34
    akk = a + sum(gammakk,1);%36
    Vkk=V; %Set Vkk to predicted V initially
    for j = 1 : Nmeas
        tmpCov = zeros(ny);
        for ell = 1 : L
            tmpCov = tmpCov + gammakk(j,ell) * invsX(:,:,ell);%30
            Vkk(:,:,ell)=Vkk(:,:,ell)+(gammakk(j,ell)/s)*Wjl(:,:,j,ell);%35
        end
        Pzkk(:,:,j) = inv(invR + tmpCov);%30
        zkk(:,j) = Pzkk(:,:,j) * (invR * Y(:,j) + invsXHxkk_*gammakk(j,:)');%31
    end

    invPkk=invP;% Set posterior information matrix to prior information matrix initially
    invPkkxkk=invP*x;% Set posterior information vector to prior information vector initially
    for ell=1:L
        invPkkxkk=invPkkxkk+H(:,:,ell)'*invsX(:,:,ell)*zkk_*gammakk(:,ell);%42
        invPkk=invPkk+sum(gammakk(:,ell))*H(:,:,ell)'*invsX(:,:,ell)*H(:,:,ell);%43
    end
    Pkk=inv(invPkk);
    xkk=Pkk*invPkkxkk;

    convergenceStat=max(max(abs(xkk-oldxkk)));

    xkk_ = xkk;
    Pkk_ = Pkk;
    Vkk_ = Vkk;
    vkk_ = vkk;
    zkk_ = zkk;
    Pzkk_ = Pzkk;
    akk_ = akk;


    if plotFlag %Plot Extent Ellipsoids with Measurements (For debugging purposes)
        plotstats.x=xkk;
        plotstats.P=Pkk;
        plotstats.v=vkk;
        plotstats.V=Vkk;
        plotstats.a=akk;

        figure(547)
        plot(Y(1,:),Y(2,:),'.')
        hold on
        plot(zkk(1,:),zkk(2,:),'r.')

        [GM, ~, ~]=convertMultiEllipseTracktoGM(plotstats);
        for ell=1:L
            drawEllipse(GM.x(:,ell),GM.P(:,:,ell),1,gca,1.5,[0, 0.4470, 0.7410]);  
        end        
        title(['Iteration No:', num2str(iter)])
        pause(0.5)
        cla
    end
end


statsout.x=xkk;
statsout.P=Pkk;
statsout.v=vkk;
statsout.V=Vkk;
statsout.a=akk;

%Calculate the predictive likelihood
logpyy=sum(logmvnpdf((Y-zkk),zeros(2,1),R));
for j=1:Nmeas
    logpyy=logpyy+0.5*log(det(2*Pzkk(:,:,j)/s))-0.5*(trace(invR*Pzkk(:,:,j))-ny);
end
for ell=1:L
    logpyy=logpyy-0.5*(vkk(ell)-ny-1)*log(det(Vkk(:,:,ell)))+0.5*(v(ell)-ny-1)*log(det(V(:,:,ell)));
    logpyy=logpyy+logGammad(0.5*(vkk(ell)-ny-1),ny)-logGammad(0.5*(v(ell)-ny-1),ny);
end
logpyy=logpyy-0.5*(trace(Pkk*invP)-nx)+0.5*log(det(2*pi*Pkk))+logmvnpdf(xkk,x,P);
logpyy=logpyy-sum(sum(exp(loggammakk).*loggammakk));
suma=sum(a);
sumakk=sum(akk);
logpyy=logpyy+gammaln(suma)-sum(gammaln(a))-gammaln(sumakk)+sum(gammaln(akk));

function outlogGammad=logGammad(x,d)
D=1:d;
outlogGammad=0.25*d*(d-1)*log(pi)+sum(gammaln(x-0.5*(D-1)));

function out=logmvnpdf(y,x,P)
Ny=size(y,2);
difference=y-repmat(x,1,Ny);
out=-0.5*log(det(2*pi*P))-0.5*sum(difference'.*(difference'/P),2);

function [logOfSum] = logSum(logA)

N=size(logA,1);

maxLogA = max(logA);
logOfSum = maxLogA+log(sum(exp(logA-repmat(maxLogA,N,1))));
