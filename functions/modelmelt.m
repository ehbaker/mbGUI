function yhat = modelmelt(beta,x)
% beta(1)=ki
% x(:,1)= PDD
% x(:,2)= bw for bw< ks*PDD
% x(:,2)= 0  for bw>=ks*PDD
% x(:,3)= 1  for bw< ks*PDD
% x(:,3)= 0  for bw>=ks*PDD
% x(:,4)= ks 
ks =x(1,4);    
yhat = ( beta(1)*x(:,3) + ks.*(1-x(:,3)) ).*x(:,1)+(1-beta(1)/ks)*x(:,2); 