clear all

addpath ..\Chehra_v0.1_MatlabFit\mex_functions
addpath ..\Chehra_v0.1_MatlabFit\
addpath ..\FaceCrop\

load Task1_data

% Write the code by yourself for the function
%% 1.	Perform face preprocessing using the provided facial image. Plot the results.
whos

% face tracking
fitting_model='../Chehra_v0.1_MatlabFit/models/Chehra_f1.0.mat';
fpt =  FaceTrack(example_img, fitting_model, 1);

% face registration
tform = cp2tform(fpt, model, 'lwm', 49);
reg_im = imtransform(example_img, tform);

% face crop
norm_face = FaceCrop(example_img, fpt, 1);

%imshow(norm_face, []);

hist(double(norm_face(:)), 255)
title('The gray histogram of cropped face')
