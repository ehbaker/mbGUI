function yhat = modelnet(beta,x)
% beta(1)=rat
% x(1,1)     = numYears, the width of each block
% x(2,1)     = numSites, the length of each block
% x(:,block1)= ks*PDD     for bw*rat>=-snmelt
% x(:,block1)= 0          for bw*rat< -snmelt
% x(:,block2)= 0          for bw*rat>=-snmelt
% x(:,block2)= -bw        for bw*rat< -snmelt
% x(:,block3)= 0          for bw*rat>=-snmelt
% x(:,block3)= ki*PDD     for bw*rat< -snmelt
% x(:,block4)= 0          for bw*rat>=-snmelt
% x(:,block4)= (ki/ks)*bw for bw*rat< -snmelt
% x(:,block5)= bw
% x(:,block6)= zone_weights
w=x(1,1);
L=x(2,1);
bs=(x(:,2:w+1)+ beta(1)*x(:,w+2:2*w+1) + x(:,2*w+2:3*w+1) + beta(1)*x(:,3*w+2:4*w+1)).*x(:,5*w+2:6*w+1);
bw=(beta(1)*x(:,4*w+2:5*w+1)).*x(:,5*w+2:6*w+1);
bs2=[];
bw2=[];
for i=1:L
	if ~isnan(bs(i,:))
 		bs2=[bs2;bs(i,:)]; %#ok<AGROW>
 	end
 	if ~isnan(bw(i,:))
 		bw2=[bw2;bw(i,:)]; %#ok<AGROW>
 	end
end
[L2,w2]=size(bs2);
yhat= ones(1,L2)*bs2*ones(w2,1) + ones(1,L2)*bw2*ones(w2,1);



