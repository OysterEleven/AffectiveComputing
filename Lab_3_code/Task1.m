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

% Concatenate multiple feature
combined_trainingData = [training_data1 training_data2];
combined_testingData = [testing_data1 testing_data2];

% linear kernel SVM for fused features
SVM = svmtrain(combined_trainingData, training_class,'showplot',true);

% test linear SVM classifier with training and testing data 
training_result = svmclassify(SVM, combined_trainingData,'showplot',true);
testing_result = svmclassify(SVM, combined_testingData,'showplot',true);

% Accuracies of  linear kernel SVM for fused features
training_acc = sum(training_class == training_result)/50
testing_acc = sum(testing_class == testing_result)/50

% Confusion matrices
training_conMat = confusionmat(training_class, training_result)
testing_conMat = confusionmat(testing_class, testing_result)

%% Optional task (cross-validation)
data = [training_data ; testing_data];
data_proso = [training_data_proso ; testing_data_proso];
class = [training_class ; testing_class];
data_personID = [training_personID ; testing_personID];
training_acc_cv = [];
testing_acc_cv = [];
for i = 0:9
    cv_train_data = data;
    cv_test_data = data(i*10+1:(i+1)*10,:);
    cv_train_data(i*10+1:(i+1)*10,:) = [];
    
    cv_train_data_proso = data_proso;
    cv_test_data_proso = data_proso(i*10+1:(i+1)*10,:);
    cv_train_data_proso(i*10+1:(i+1)*10,:) = [];
    
    cv_train_class = class;
    cv_test_class = class(i*10+1:(i+1)*10,:);
    cv_train_class(i*10+1:(i+1)*10,:) = [];
    
    [V1, eigvalue_f] = PCA(cv_train_data, options_facial);
    [V2, eigvalue_a] = PCA(cv_train_data_proso, options_audio);

    % Extract low-dimensional features for training and testing features:
    cv_train_data1 = cv_train_data * V1;
    cv_test_data1 = cv_test_data * V1;
    mean_cv1 = mean(cv_train_data1);
    std_cv1 = std(cv_train_data1);

    cv_train_data2 = cv_train_data_proso * V2;
    cv_test_data2 = cv_test_data_proso * V2;
    mean_cv2 = mean(cv_train_data2);
    std_cv2 = std(cv_train_data2);
    % Use Z-score normalization to normalize for facial expression
    cv_train_data1 = (cv_train_data1 - mean_cv1) ./ std_cv1;
    cv_test_data1 = (cv_test_data1 - mean_cv1) ./ std_cv1;
    % Use Z-score normalization to normalize for audio
    cv_train_data2 = (cv_train_data2 - mean_cv2) ./ std_cv2;
    cv_test_data2 = (cv_test_data2 - mean_cv2) ./ std_cv2;

    combined_train = [cv_train_data1 cv_train_data2];
    combined_test = [cv_test_data1 cv_test_data2];

    SVM_cv = svmtrain(combined_train, cv_train_class,'showplot',true);
    training_result_cv = svmclassify(SVM_cv, combined_train,'showplot',true);
    testing_result_cv = svmclassify(SVM_cv, combined_test,'showplot',true);
    training_acc_cv(i+1) = sum(cv_train_class == training_result_cv)/90;
    testing_acc_cv(i+1) = sum(cv_test_class == testing_result_cv)/10;

end
% mean and std of CV fold performances
mean_train_cv = mean(training_acc_cv)
mean_test_cv = mean(testing_acc_cv)
std_train_cv = std(training_acc_cv)
std_test_cv = std(testing_acc_cv)