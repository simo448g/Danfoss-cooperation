function [QR] =lifting_QR(Q,R,Hp)
%Lifting Q and R, and putte them into one matrix. 
%Can take both Q and R as matrix or as vectors 

%% Checing if R and Q is on matrix form
% if this is not the case, they are put on matrix form
if size(Q,1) == 1 
    Q=diag(Q); 
elseif size(Q,2) == 1 
    Q=diag(Q); 
end 

if size(R,1) == 1 
    R=diag(R); 
elseif size(R,2) == 1 
    R=diag(R); 
end 
%% Predefine the size of Qlift and Rlift
Qlift=zeros(size(Q,1)*Hp+size(Q,1),size(Q,1)*Hp+size(Q,1)); 

Rlift=zeros(size(R,1)*Hp,size(R,1)*Hp); 

%% Putting Q on the diagnoal of the lift matrix
Qsize=size(Q,1);

for i=1:Qsize:size(Qlift,1)
    Qlift(i:i+Qsize-1,i:i+Qsize-1)=Q;
end 

%% Putting R on the diagnoal of the lift matrix

Rsize=size(R,1);
for i=1:Rsize:size(Rlift,1)
    Rlift(i:i+Rsize-1,i:i+Rsize-1)=R;
end 

%% Collecting it all into one matrix: 

QR=[Qlift, zeros(size(Qlift,1),size(Rlift,1)); zeros(size(Rlift,1),size(Qlift,1)),Rlift]; 



end 