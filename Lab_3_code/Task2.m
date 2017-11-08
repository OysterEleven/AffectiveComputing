clear all
clc
load lab3_data

% Set reduced dimension
options_facial = [];
options_facial.ReducedDim = 20;
options_audio = [];
options_audio.ReducedDim = 15;

% Construct two subspace U1 and U2 for facial expression features and audio
% features

[U1, eigvalue_f] = PCA(training_data, options_facial);
[U2, eigvalue_a] = PCA(training_data_proso, options_audio);

% Extract low-dimensional features for training and testing features:
training_data1 = training_data * U1;
testing_data1 = testing_data * U1;
mean1 = mean(training_data1);
std1 = std(training_data1);

training_data2 = training_data_proso * U2;
testing_data2 = testing_data_proso * U2;
mean2 = mean(training_data2);
std2 = std(training_data2);

% Use Z-score normalization to normalize for facial expression
training_data1 = (training_data1 - mean1) ./ std1;
testing_data1 = (testing_data1 - mean1) ./ std1;
% Use Z-score normalization to normalize for audio
training_data2 = (training_data2 - mean2) ./ std2;
testing_data2 = (testing_data2 - mean2) ./ std2;

% Use ?anoncorr?to construct the Canonical Projective Vector (CPV)
[A, B, r] = canoncorr (training_data1, training_data2);
correlation_coefficient = r

% Construct Canonical Correlation Discriminant Features (CCDF)
CCA_training_data1=training_data1*A;
CCA_testing_data1=testing_data1*A;
CCA_training_data2=training_data2*B;
CCA_testing_data2=testing_data2*B;

% Concatenate multiple feature
training_CCDF = [CCA_training_data1 CCA_training_data2];
testing_CCDF = [CCA_testing_data1 CCA_testing_data2];

% linear kernel SVM for fused features
SVM = svmtrain(training_CCDF, training_class,'showplot',true);

% test linear SVM classifier with training and testing data 
training_result = svmclassify(SVM, training_CCDF,'showplot',true);
testing_result = svmclassify(SVM, testing_CCDF,'showplot',true);

% Accuracies of  linear kernel SVM for fused features
training_acc = sum(training_class == training_result)/50
testing_acc = sum(testing_class == testing_result)/50

% Confusion matrices 
training_conMat = confusionmat(training_class, training_result)
testing_conMat = confusionmat(testing_class, testing_result)
