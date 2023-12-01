function [WLift,wLift] =lifting_slew_rate_constraints(Hc,W,w)
%In this function the slew rate constraints is lifted. 

%% First W is lifted: 
WLift=zeros(Hc*size(W,1),Hc*size(W,2)); 


indexCol=1;
for i=1:size(W,1):size(WLift,1)
    WLift(i:i+size(W,1)-1,indexCol:indexCol+size(W,2)-1)=W;
    indexCol=indexCol+size(W,2);
end 


%% Next w is lifted: 
wLift=zeros(size(WLift,1),size(w,2)); 
for i=1:size(w,1):size(wLift,1)
    wLift(i:i+size(w,1)-1,:)=w;
end 



end