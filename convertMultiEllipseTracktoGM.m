function [GM, meanTarget, covTarget]=convertMultiEllipseTracktoGM(stats)

ny=size(stats.V,1);
L=length(stats.x)/ny-1;
for i=1:L
    if i==1
        GM.x(:,1)=stats.x(1:ny);
    else
        GM.x(:,i)=stats.x(1:ny)+stats.x(ny*i+1:ny*i+ny);
    end
    GM.P(:,:,i)=stats.V(:,:,i)/(stats.v(i)-2*ny-2);
end

GM.pi=stats.a/sum(stats.a);

meanTarget=GM.x*GM.pi';
covTarget=zeros(ny);
for i=1:L
    covTarget=covTarget+GM.pi(i)*(GM.P(:,:,i)+(GM.x(:,i)-meanTarget)*(GM.x(:,i)-meanTarget)');
end