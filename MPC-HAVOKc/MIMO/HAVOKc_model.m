function [Av, Bv, Pup, Pdown,S,Sr] =HAVOKc_model(x,u,nd,r)
%x=data
%u=inputs 
%nd=embedding length 
%r=truncation value for the hankel matrix 

%Av=system matrix in v space 
%Bv input matrix in v space 
%Pup projection matrix from v space to hankel 
%Pdown projection matrix from hankel space to v space 
%S all singular values of the Hankel matrix of the measured data

%% Starting by making the hankel 


%measured values/states: 
index=1;
n=size(x,1); %Determine the amount of states 
for i=0:nd
    H(index:index+n-1,:)=x(:,nd-i+1:end-i);
    index=index+n;
end
%input
index=1;
nu=size(u,1); %Determine the amount of states 

for i=0:nd
    Hu(index:index+nu-1,:)=u(:,nd-i+1:end-i);
    index=index+nu;
end


%% Taking the SVD of the Hankel matrix for measurement and truncation it  
[U,S,V]=svd(H); 

Ur=U(:,1:r);
Sr=S(1:r,1:r);
Vr=V(:,1:r);

%% Creating the projection matrix  Pup and Pdown, based on the truncated SVD
Pdown=pinv(Ur*Sr); 
Pup=Ur*Sr; 
%% Projection Hu to V space 
Vu=Pdown*Hu; 

%% Creating V_(r,k) and V_(r,k+1) 
VrCon=Vr(1:end-1,:);
VrPred=Vr(2:end,:);

%% Making the G matrix: 

G=[transpose(VrCon);Vu(:,1:end-1)]; 


%% Approximation Av and Bv

dummy=transpose(VrPred)*pinv(G); 
Av=dummy(:,1:r); 

Bv=dummy(:,r+1:end); 







end 
