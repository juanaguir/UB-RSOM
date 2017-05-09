function [ alphac] = alphacalc(xhigh,xlow)
%This function calculates the alpha values to perform frequency band
%equalization. (methods, equation 2)
%input: 
% xhigh.  Hi frequency data
% xlow.   Lo9w frequency darta
 alphaval=[0:0.01:100];
 differ=zeros(size(alphaval));
 for(i=1:length(alphaval))
 differ(i)=sum(sum((xlow - alphaval(i)*xhigh).^2));
 end
  [dum,ind]=min(differ);
  alphac=alphaval(ind);
figure; plot(alphaval,differ);
end

