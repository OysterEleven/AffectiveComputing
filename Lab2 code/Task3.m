%% Task 3
load lab2_data.mat;
% Train an SVM with the prosody data using the ?raining_data_proso?features and a 3rd order polynomial kernel.
SVM_proso = svmtrain(training_data_proso, training_class,'kernel_function','polynomial','showplot',true);

% test the prosodic feature based SVM classifier with training and testing data 
training_proso_result = svmclassify(SVM_proso, training_data_proso,'showplot',true);
testing_proso_result = svmclassify(SVM_proso, testing_data_proso,'showplot',true);

% Accuracies of 3rd order polynomial kernel SVM with prosodic features
training_proso_acc = sum(training_class == training_proso_result)/50
testing_proso_acc = sum(testing_class == testing_proso_result)/50



% Train an SVM with the MFCC data using the ?raining_data_mfcc?features and a 3rd order polynomial kernel.
SVM_mfcc = svmtrain(training_data_mfcc, training_class,'kernel_function','polynomial','showplot',true);

% test the MFCC feature based SVM classifier with training and testing data 
training_mfcc_result = svmclassify(SVM_mfcc, training_data_mfcc,'showplot',true);
testing_mfcc_result = svmclassify(SVM_mfcc, testing_data_mfcc,'showplot',true);

% Accuracies of 3rd order polynomial kernel SVM with MFCC derived features
training_mfcc_acc = sum(training_class == training_mfcc_result)/50
testing_mfcc_acc = sum(testing_class == testing_mfcc_result)/50

% Confusion matrices of 3rd order polynomial kernel SVM with prosodic features (training and testing data)
training_proso_conMat = confusionmat(training_class, training_proso_result)
testing_proso_conMat = confusionmat(testing_class, testing_proso_result)

% Confusion matrices of 3rd order polynomial kernel SVM with MFCC derived features (training and testing data)
training_mfcc_conMat = confusionmat(training_class, training_mfcc_result)
testing_mfcc_conMat = confusionmat(testing_class, testing_mfcc_result)



%% Optional task (cross-validation)
data = [training_data_proso 
    testing_data_proso];
class = [training_class
    testing_class];
data_personID = [training_personID
    testing_personID];
training_acc = [];
testing_acc = [];
for i = 0:9
    cv_data = data;
    cv_data(i*10+1:(i+1)*10,:) = [];
    cv_class = class;
    cv_class(i*10+1:(i+1)*10,:) = [];
    SVM = svmtrain(cv_data, cv_class,'kernel_function','polynomial','showplot',true);
    training_result = svmclassify(SVM, cv_data,'showplot',true);
    testing_result = svmclassify(SVM, data(i*10+1:(i+1)*10,:),'showplot',true);
    training_acc(i+1) = sum(cv_class == training_result)/90
    testing_acc(i+1) = sum(class(i*10+1:(i+1)*10,:) == testing_result)/10

end
% mean and std of CV fold performances
mean_train_cv = mean(training_acc)
mean_test_cv = mean(testing_acc)
std_train_cv = std(training_acc)
std_test_cv = std(testing_acc)