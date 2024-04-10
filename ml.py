
######################################################


import pandas as pd

landmark_gene = pd.read_csv("/home/albin/Desktop/New_thesis/excel/goid_6357.csv")

cleaned_landmark = landmark_gene.drop(columns=['GOID', 'Genes', 'Term'])
cleaned_landmark = cleaned_landmark.drop_duplicates(subset='EntrezgeneID', keep='first')
cleaned_landmark.set_index('EntrezgeneID', inplace= True, drop= True)

###################################################

target_gene = pd.read_csv("/home/albin/Desktop/New_thesis/excel/target_6357.csv",engine='python')
cleaned_target= target_gene.drop(columns=['Genes'])
cleaned_target.set_index('EntrezgeneID', inplace= True, drop= True)

###########################################################

training_data = cleaned_landmark.transpose()
testing_data  = cleaned_target.transpose()

X =training_data.values
y=testing_data.values

#####################################################################

from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.datasets import make_regression
from sklearn.preprocessing import StandardScaler

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2, random_state=42)

# Standardize features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# Initialize the Linear Regression model
model = LinearRegression()

# Train the model
model.fit(X_train_scaled, y_train)

# Validate the model
predictions = model.predict(X_test_scaled)
mse = mean_squared_error(y_test, predictions)
print(f'Mean Squared Error: {mse}')


############################################
#Predict on test data
predictions = model.predict(X_test_scaled[:5])
print("Predicted values are: ", predictions)
print("Real values are: ", y_test[:5])
##############################################

############################################################

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import KFold

# Assuming X_train_scaled and y_train are already defined and preprocessed

# Define the number of output neurons
output_neurons = y_train.shape[1]

# Define the number of splits for cross-validation
n_splits = 5
kf = KFold(n_splits=n_splits)

# Initialize variables to store the best score and the best model
best_score = float('inf')
best_model = None

# Perform cross-validation
for train_index, val_index in kf.split(X_train_scaled):
    # Split data into training and validation sets
    X_train_fold, X_val_fold = X_train_scaled[train_index], X_train_scaled[val_index]
    y_train_fold, y_val_fold = y_train[train_index], y_train[val_index]

    # Build a neural network model with hyperparameters for each fold
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(128, activation='relu', input_shape=(X_train_scaled.shape[1],)),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(64, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(output_neurons)
    ])

    # Compile the model
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.001, beta_1=0.5)
    model.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['mae'])

    # Fit data to model
    history = model.fit(X_train_fold, y_train_fold,
                        batch_size=64,
                        epochs=50,
                        verbose=1,
                        validation_data=(X_val_fold, y_val_fold))

    # Evaluate the model on the validation set
    score = model.evaluate(X_val_fold, y_val_fold, verbose=0)

    # If the score is better than the best score, update the best score and best model
    if score[1] < best_score:
        best_score = score[1]
        best_model = model

# Use the best model for prediction
predictions = best_model.predict(X_test_scaled)


# Print the best validation mean absolute error
print(f"Best validation MAE: {best_score}")



from matplotlib import pyplot as plt
#plot the training and validation accuracy and loss at each epoch
loss = best_model.history['loss']
val_loss = best_model.history['val_loss']
epochs = range(1, len(loss) + 1)
plt.plot(epochs, loss, 'y', label='Training loss')
plt.plot(epochs, val_loss, 'r', label='Validation loss')
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.show()





############################################
#Predict on test data
predictions =best_model.predict(X_test_scaled[:5])
print("Predicted values are: ", predictions)
print("Real values are: ", y_test[:5])
##############################################

print(best_model.metrics_names)
#Comparison with other models..
#Neural network - from the current code
mse_neural, mae_neural = best_model.evaluate(X_test_scaled, y_test)
print('Mean squared error from neural net: ', mse_neural)
print('Mean absolute error from neural net: ', mae_neural)

########################################################

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error


model = RandomForestRegressor(n_estimators = 30, random_state=30)
model.fit(X_train_scaled, y_train)

y_pred_RF = model.predict(X_test_scaled)

mse_RF = mean_squared_error(y_test, y_pred_RF)
mae_RF = mean_absolute_error(y_test, y_pred_RF)
print('Mean squared error using Random Forest: ', mse_RF)
print('Mean absolute error Using Random Forest: ', mae_RF)

param_grid = {
    'n_estimators': [100, 200, 300, 400, 500],
    'max_depth': [None, 10, 20, 30, 40, 50],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['auto', 'sqrt', 'log2']
}

# Initialize the grid search model
grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='neg_mean_squared_error', n_jobs=-1)

# Fit the grid search to the data
grid_search.fit(X_train_scaled, y_train)

# Print the best parameters
print("Best parameters found: ", grid_search.best_params_)

y_pred_RF = grid_search.predict(X_test_scaled)

mse_RF = mean_squared_error(y_test, y_pred_RF)
mae_RF = mean_absolute_error(y_test, y_pred_RF)
print('Mean squared error using Random Forest: ', mse_RF)
print('Mean absolute error Using Random Forest: ', mae_RF)

#Feature ranking...
import pandas as pd
feature_list = list(X.columns)
feature_imp = pd.Series(model.feature_importances_, index=feature_list).sort_values(ascending=False)
print(feature_imp)
