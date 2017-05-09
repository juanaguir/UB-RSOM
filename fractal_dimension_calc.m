function [cden] = fractal_dimension_cal(R,threshold)
%% this function makes a threshold to an image, binarize it and then calculates the fractal dimension
%% Input: R. Reconstructed image
%%        Threshold: threshold value
%% Output: cden, fractal number

index= find(R<max(max(max(R)))/threshold);
index2= find(R>=max(max(max(R)))/threshold);
R(index)=0;
R(index2)=1;

cden=hausDim(R);


end






