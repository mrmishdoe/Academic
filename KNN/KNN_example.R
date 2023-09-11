#K nearest neighbors example BIOS 694

heart_data <- read.csv("cleveland_data.csv", header = F)

names(heart_data) <- c("Age", "Sex", "Chest_pain_type", "Blood_pressure", "cholestoral", "Blood_sugar", 
                       "electrocardiographic", 
                       "Max_heart_rate", "Excercise_induced_angina", 
                       "Oldpeak", "Slope", "Major_Vessels", "thal", "Disease")

#Model matrix is all but last column of heart_data
x_data <- heart_data[, -dim(heart_data)[2]]
head(x_data)

#Normalize data so everything is in same scale (important for KNN because distance with will be dominated by features with large scale)

x_data_scale= scale(x_data)
head(x_data_scale)

#Create y labels, binary so disease indicated by 1, absence by 0
y_label <- heart_data[, dim(heart_data)[2]]
y_label_binary <- ifelse(y_label>0, 1, 0)

#Train/Test split

set.seed(2023)
library(caTools)
train_split <- 0.7
sample_indices <- sample.split(heart_data[,1], SplitRatio = 0.70)

x_train = as.matrix(x_data_scale[sample_indices,])
y_train = as.vector(y_label_binary[sample_indices])
x_test = as.matrix(x_data_scale[-sample_indices,])
y_test = as.vector(y_label_binary[-sample_indices])

#Train KNN

library(class)

test_pred <- knn(train = x_train, test = x_test,cl = y_train, k=5)

#Evaluate model via confusion matrix

actual <- y_test

confusion_matrix <- table(actual,test_pred)
confusion_matrix
# Can calculate accuracy via dividing sum of diagonal elements by length of y_test

accuracy <- sum(diag(confusion_matrix))/length(actual)
cat("Accuracy: ", accuracy*100)

# Different values of K

#K = 3
test_pred <- knn(
  train = x_train, 
  test = x_test,
  cl = y_train, 
  k=3
)

actual <- y_test

confusion_matrix <- table(actual,test_pred)
confusion_matrix
accuracy <- sum(diag(confusion_matrix))/length(actual)
cat("Accuracy: ", accuracy*100)

#K = 20
test_pred <- knn(
  train = x_train, 
  test = x_test,
  cl = y_train, 
  k=20
)

actual <- y_test

confusion_matrix <- table(actual,test_pred)
confusion_matrix
accuracy <- sum(diag(confusion_matrix))/length(actual)
cat("Accuracy: ", accuracy*100)

#Model evaluation with caret

library(caret)
heart_data$disease_binary <- ifelse (heart_data$Disease ==0 , 0, 1)
heart_data$disease_binary <- factor(heart_data$disease_binary)

heart_data <- heart_data[,-14]


trainIndex<-caret::createDataPartition(heart_data[,"disease_binary"], 
                                   times=1, 
                                   p = .75, 
                                   list = FALSE)
train <- heart_data[trainIndex, ]
test <- heart_data[-trainIndex, ]

#Preprocessing (scaling)

preProcValues <- preProcess(train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, train)
testTransformed <- predict(preProcValues, test)

#Model tuning
knnModel <- train(
  disease_binary ~ ., 
  data = trainTransformed, 
  method = "knn", 
  trControl = trainControl(method = "cv"), 
  tuneGrid = data.frame(k = c(3,5,7))
)
knnModel

#Training best performing model
best_model<- knn3(
  disease_binary ~ .,
  data = trainTransformed,
  k = knnModel$bestTune$k
)
best_model

#Model evaluation

predictions <- predict(best_model, testTransformed, type = "class")
# Calculate confusion matrix
cm <- confusionMatrix(predictions, testTransformed$disease_binary)
cm
