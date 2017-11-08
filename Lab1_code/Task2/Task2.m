addpath ..\LBP

load Task2_data
% Use LBP (lbp.m and getmapping.m) for one image
P = 8;
R = 1;
mapping = getmapping(P, 'u2'); % help getmapping

H1 = lbp(norm_face,R,P,mapping,'hist');
subplot(2,1,1),stem(H1);
title('The LBP histogram of cropped face')

H2 = lbp(norm_face);
subplot(2,1,2),stem(H2);

SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
I2 = lbp(norm_face,SP,0,'i');
imshow(I2, []);
%LBP code image using sampling points in SP and no mapping. Now H2 is equal to histogram of I2.
