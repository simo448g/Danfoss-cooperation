function [F_eta, Fu, F0, b,F_r] =lifting_dynamics(Av,Bv,Pup,Pdown,nd,r, Hp,Hc,nu,DeltaZ0,eta_0,RefV2,ny)
%Output is the lifted matrixes, F_eta, Fu, F0,b, F_r 

%As input it is the dynamics system matrix Av and input Bv, which 
%Dynamics unfolded on the V subspace. Futhermore it is the projection 
%Matrix to and from the V subspace to the hankel space Pup and Pdown. 

%Futhermore the embedding delay nd is nedded, 
%r=SVD trunucation value
%Hp=prediction horizion 
%nu=number of inputs
%DeltaZ0 inital input sekvens in the hankel domain 
%eta_0 inital "state" value 
% Ref value for the system
%ny=antal målingerne


%% 

%Setting up Cc
n=size(Av,1); %Number of states
Cc=[];
for i=1:ny
Cc=[Cc; zeros(1,i-1),1,zeros(1,nd*ny-i+ny)]; %No idea if Cc is correct 
end 

%Cc=[ones(1,1),zeros(1,nd-1)];
%% 
An=[Av, zeros(n,ny);Cc*Pup*Av, eye(ny) ];
Bn=[Bv; Cc*Pup*Bv]; 






%% Starting by making F_eta 
%First the dignoal consisting of -1 i made 
RowAn=size(An,1);

F_eta=-eye(Hp*RowAn+RowAn); %Was -eye((Hp+1)*RowAn);

%Afterward An is added below the diagnaol

index=RowAn+1; 
for i=1:Hp
    F_eta(index:index+RowAn-1,index-RowAn:index-1)=An;
    index=index+RowAn; 
end 


%% Making Fu
Fu=zeros(size(F_eta,1),nu*Hc); %Was zeros(size(F_eta,1),nu*size(F_eta,2)-RowAn);
% Making the BnPdownI_i
Bnu=[];
for i=1:nd+1
    Bnu=[Bnu;Bn*Pdown*[zeros(i-1,nu);eye(nu);zeros((nd+1)*nu-i-(nu-1),nu)]];
end 

%First the "diagnoal" is made consisting of Bn*Pdown*Ik
indexRow=RowAn+1; %The first ones needs to be zero
indexCol=1;

for i=1:Hp %was size(Fu,1)-size(Bnu,1)-RowAn+1
    Fu(indexRow:indexRow+size(Bnu,1)-1,indexCol:indexCol+size(Bnu,2)-1)=Bnu;
    indexRow=indexRow+size(An,1); 
    %indexRow=indexRow+size(Bn*Pdown*I1,1); 
    indexCol=indexCol+size(Bnu,2); 
end 

Fu=Fu(1:size(F_eta,1),:);

% %Making the last part where less and less is used 
% for i=1:size(Bnu,1)-1
%     Fu(indexRow:indexRow+size(Bnu,1)-1-i,indexCol:indexCol+size(Bnu,2)-1)=Bnu(1:end-i,:);
%     indexRow=indexRow+nu; 
%     %indexRow=indexRow+size(Bn*Pdown*I1,1); 
%     indexCol=indexCol+size(Bnu,2); 
% end 
% 



%% Making F_0 

%First D, is design 
D=zeros((nd+1)*nu,(nd+1)*nu); 

%Making the off diagonal one: 
for i=2:(nd)*nu+1
    D(i,i-1)=1; 
end 


F0=zeros(Hp*RowAn+RowAn,1);
indexRow=RowAn+1; 
for i=1:Hp
    %F0(indexRow:indexRow+RowAn-1,:)=Bn*Pdown*D^i*DeltaZ0;
    F0(indexRow:indexRow+RowAn-1,1)=Bn*Pdown*D^i*DeltaZ0;
    %F0(indexRow,:)=Bn*Pdown*D^i*DeltaZ0;
    indexRow=indexRow+RowAn;
end 

%% Making b 
b=[-eta_0;zeros(size(F_eta,1)-size(eta_0,1),1)];

%% Making the refarence change F_r 
F_r=zeros(size(F_eta,1),1);
%Making sure Ref is the right dimension! 
% if size(Ref,1)>size(Ref,2)
%     Ref=Ref';
% end 

index=1; 
for i=2*RowAn+ny+1:RowAn:size(F_eta,1)%-RowAn%i=RowAn+r+1:RowAn:size(F_eta,1)%-RowAn
    F_r(i:i+ny-1,1)=RefV2(:,index);
    index=index+1;
end 



end
