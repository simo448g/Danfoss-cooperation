function [F,f1,fLift,Fi] =lifting_input_constraints(Hp,Fis,f)

%First the matrix Fi (with i*nu-nu:i*nu being the colums corresponding to
%the time i) 

%Fi=zeros(Hc*size(Fis,1),Hc*size(Fis,2));
Fi=zeros(Hp*size(Fis,1),Hp*size(Fis,2));

indexCol=1;
for i=1:size(Fis,1):size(Fi,1)
    Fi(i:i+size(Fis,1)-1,indexCol:indexCol+size(Fis,2)-1)=Fis;
    indexCol=indexCol+size(Fis,2);
end 

%% Making F, by making the addition as described in the report  

F=zeros(size(Fi,1),size(Fi,2)); 


nu=size(Fis,2);
indexCol=1;
for i=1:size(Fis,2):size(F,2)% i=1:size(Fis,2):size(F,2)
    %Dummy used to sum the columens of Fi
    dummy=zeros(size(Fi,1),size(Fis,2));
    %summing the columens of Fi
    for ii=i:nu:size(Fi,2)-nu
        dummy=dummy+Fi(:,ii:ii+size(Fis,2)-1);
    end 
    %If it is the first time summing the answear is also f1 
    if i==1
        f1=dummy;
    end 

    F(:,i:i+size(Fis,2)-1)=dummy;
    indexCol=indexCol+size(Fis,2);
end 

%% Next it is desired to lift f

fLift=zeros(size(F,1),size(f,2)); 
for i=1:size(f,1):size(fLift,1)
    fLift(i:i+size(f,1)-1,:)=f;
end 

end 
