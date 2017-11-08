load Task3_data
whos
% linear kernel SVM for LBP-TOP features
SVM_lbp = svmtrain(training_data, training_class,'showplot',true);

% test linear SVM classifier with training and testing data 
training_result = svmclassify(SVM_lbp, training_data,'showplot',true);
testing_result = svmclassify(SVM_lbp, testing_data,'showplot',true);

% Accuracies of  linear kernel SVM with LBP-TOP features
training_acc = sum(training_class == training_result)/50
testing_acc = sum(testing_class == testing_result)/50

% Confusion matrices of 3rd order polynomial kernel SVM with lindic features (training and testing data)
training_conMat = confusionmat(training_class, training_result)
testing_conMat = confusionmat(testing_class, testing_result)
