#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 08:41:55 2024

@author: Manasa Kesapragada
"""

import pandas as pd
import numpy as np
import math
from keras import backend as K
from keras.models import Sequential
from keras.models import load_model
from keras.layers import LSTM, Dense, Dropout, Bidirectional
from keras.optimizers.legacy import Adam
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import os
import random
import tensorflow as tf
import pdb
from keras.regularizers import l2

base = "/Users/alexandranava/Desktop/Spores/"
experiment1 = "M4576_s7" # 135
experiment2 = "M4581_s7" # 107 
df_test = pd.read_csv('/Users/kesaprm/FY23/Summer23/ARO_Project/Train1_Model_Data.csv') #10mMl Set1
df_train = pd.read_csv('/Users/kesaprm/FY23/Summer23/ARO_Project/Test2_Model_Data.csv') #10mMl Set2



columns_to_drop = ['X_POSITION', 'Y_POSITION',  'ELLIPSE_MINOR','ELLIPSE_MAJOR']
df_train = df_train.drop(columns=columns_to_drop)
df_test = df_test.drop(columns=columns_to_drop)
#df['drug'] = df['drug'].fillna(0)



# Reset particle IDs to start from 1
df_train['new_spore_id'], unique_ids = pd.factorize(df_train['SPORE_ID'])
df_train['new_spore_id'] = df_train['new_spore_id'] + 1

df_test['new_spore_id'], unique_ids = pd.factorize(df_test['SPORE_ID'])
df_test['new_spore_id'] = df_test['new_spore_id'] + 1

df_train = df_train.drop(['Unnamed: 0'], axis=1)
df_test = df_test.drop(['Unnamed: 0'], axis=1)


#plt.hist(df['ecc'])
#plt.hist(df['size'])
plt.hist(df_train['PERIMETER'])

input_cols = ['INTENSITY','AREA','GERMINANT_EXPOSURE','PERIMETER', 'CIRCULARITY' , 'GERMINATION']
output_cols = 'GERMINATION'


def create_dataset(df, cells, lookback=3, in_cols= input_cols, out_cols=output_cols, drug_conc=[1]):
    trainX, trainY = [], [] #lists of training and testing inputs/outputs
    for drug in drug_conc:
        for track in range(cells[0], cells[1]):
            cell = df.loc[(df["new_spore_id"] == track)] #all rows of data pertaining to this cell
            cell = cell[in_cols] #reduce it to our columns of interest
            for i in range(len(cell)-lookback-1):
                trainX.append(cell[i:i+lookback])
            cell = cell[out_cols]
            #pdb.set_trace()
            for i in range(len(cell)-lookback-1):
                trainY.append(cell[i+lookback+1:i+lookback+2])

    trainX = np.array(list(map(lambda x: x.to_numpy(), trainX)))
    trainY = np.array(list(map(lambda x: x.to_numpy(), trainY)))
    return np.array(trainX), np.array(trainY)

# train with 80% val with 20%

trainX, trainY = create_dataset(df_train,cells=(1,66))
valX, valY = create_dataset(df_train, cells=(66,76))
testX, testY = create_dataset(df_test, cells=(1,54))


#trains numerous models using a list of numbers to initialize the models
def typical_model(trainX,trainY,valX,valY,testX,testY,numbers):
    models = [] #list of models
    predictions = [] #list of prediction vectors

    for i in numbers:
        #build the model
        print('Training model number {}'.format(i))
        model = Sequential()
        model.add(LSTM(80, input_shape=(trainX.shape[1], trainX.shape[2]))) #(Lookback of 1)
        model.add(Dense(1, activation='sigmoid'))
        model.add(Dropout(0.01))
        model.compile(loss='binary_crossentropy', optimizer=Adam(), metrics=['accuracy'])
        
        early_stopping = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3, restore_best_weights=True)
        
        history = model.fit(trainX, trainY, validation_data=(valX, valY), epochs=30, batch_size=1, verbose=1, callbacks=[early_stopping])
        models.append(model)
        model.summary()

        model.save("initial_{}.h5".format(i))
        print("Saved model {} to disk".format(i))
        
        # plot training and validation loss
        loss = history.history['loss']
        val_loss = history.history['val_loss']
        epochs = range(1, len(loss) + 1)
        
        plt.figure()
        plt.plot(epochs, loss, 'bo', label='Training loss')
        plt.plot(epochs, val_loss, 'b', label='Validation loss')
        plt.title('Training and validation loss for model {}'.format(i))
        plt.legend()
        plt.show()
        
        # plot training and validation acc
        loss = history.history['accuracy']
        val_loss = history.history['val_accuracy']
        epochs = range(1, len(loss) + 1)
        
        plt.figure()
        plt.plot(epochs, loss, 'bo', label='Training accuracy')
        plt.plot(epochs, val_loss, 'b', label='Validation accuracy')
        plt.title('Training and validation accuracy for model {}'.format(i))
        plt.legend()
        plt.show()

        #Predict on training, validation, and test sets
        trainPredict = model.predict(trainX)
        valPredict = model.predict(valX)
        testPredict = model.predict(testX)

        # Threshold the predictions to get binary outputs
        trainPredict_binary = (trainPredict > 0.5).astype(int)
        valPredict_binary = (valPredict > 0.5).astype(int)
        testPredict_binary = (testPredict > 0.5).astype(int)

        # Calculate RMSEs
        trainScore = math.sqrt(mean_squared_error(trainY, trainPredict_binary))
        print('Training RMSE: {}'.format(trainScore))
        valScore = math.sqrt(mean_squared_error(valY, valPredict_binary))
        print('Validation RMSE: {}'.format(valScore))
        testScore = math.sqrt(mean_squared_error(testY, testPredict_binary))
        print('Testing RMSE: {}'.format(testScore))

        # add predictions to list of prediction vectors
        preds = np.concatenate([trainPredict_binary, valPredict_binary, testPredict_binary], axis=None)
        predictions.append(preds)


    #find the average of all prediction vectors
    mean_pred = np.mean(predictions,axis=0)

    #find the prediction vector that is closest to the mean
    closest = 0 #index of the prediction vector that is closest to the mean
    dist = 100 #the distance of that closest vector to the mean vector
    for i in range(len(predictions)):
        thisdist = math.sqrt(mean_squared_error(predictions[i], mean_pred))

        print('Model number: {}, distance from mean: {}'.format(i,thisdist))

        if thisdist < dist:
            dist = thisdist
            closest = i

    #return the "most average" model
    print('Returning model {}, whose distance is {}'.format(closest,dist))
    return models, closest


models, closest = typical_model(trainX,trainY,valX,valY,testX,testY,range(1))
model = models[closest]
#model.save("alex_data_model_lk2.h5")
#model.save("alex_data_model_lk2.keras")
print("Saved model to disk")


#Predict on training, validation, and test sets
trainPredict = model.predict(trainX)
trainPredict =(trainPredict > 0.5).astype(int)
valPredict = model.predict(valX)
valPredict =(valPredict > 0.5).astype(int)
testPredict = model.predict(testX)
testPredict =(testPredict > 0.5).astype(int)


# #Calculate RMSEs
# trainScore = math.sqrt(mean_squared_error(trainY, trainPredict))
# print('Training RMSE: {}'.format(trainScore))
# valScore = math.sqrt(mean_squared_error(valY, valPredict))
# print('Validation RMSE: {}'.format(valScore))
# testScore = math.sqrt(mean_squared_error(testY, testPredict))
# print('Testing RMSE: {}'.format(testScore))

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# Calculate accuracy, precision, recall, and F1-score
train_accuracy = accuracy_score(trainY, trainPredict)
val_accuracy = accuracy_score(valY, valPredict)
test_accuracy = accuracy_score(testY, testPredict)
print('Training Accuracy: {}'.format(train_accuracy))
print('Validation Accuracy: {}'.format(val_accuracy))
print('Testing Accuracy: {}'.format(test_accuracy))

train_precision = precision_score(trainY, trainPredict)
val_precision = precision_score(valY, valPredict)
test_precision = precision_score(testY, testPredict)
print('Training Precision: {}'.format(train_precision))
print('Validation Precision: {}'.format(val_precision))
print('Testing Precision: {}'.format(test_precision))

train_recall = recall_score(trainY, trainPredict)
val_recall = recall_score(valY, valPredict)
test_recall = recall_score(testY, testPredict)
print('Training Recall: {}'.format(train_recall))
print('Validation Recall: {}'.format(val_recall))
print('Testing Recall: {}'.format(test_recall))

train_f1 = f1_score(trainY, trainPredict)
val_f1 = f1_score(valY, valPredict)
test_f1 = f1_score(testY, testPredict)
print('Training F1-Score: {}'.format(train_f1))
print('Validation F1-Score: {}'.format(val_f1))
print('Testing F1-Score: {}'.format(test_f1))


from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve

# Plot ROC curve
def plot_roc_curve(y_true, y_pred, title):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    auc_score = roc_auc_score(y_true, y_pred)
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (area = {auc_score:0.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.show()

# Training set
train_auc = roc_auc_score(trainY, trainPredict)
print('Training AUC: {}'.format(train_auc))
plot_roc_curve(trainY, trainPredict, 'ROC curve - Training set')

# Validation set
val_auc = roc_auc_score(valY, valPredict)
print('Validation AUC: {}'.format(val_auc))
plot_roc_curve(valY, valPredict, 'ROC curve - Validation set')

# Test set
test_auc = roc_auc_score(testY, testPredict)
print('Testing AUC: {}'.format(test_auc))
plot_roc_curve(testY, testPredict, 'ROC curve - Testing set')

# Confusion matrices
train_cm = confusion_matrix(trainY, trainPredict)
val_cm = confusion_matrix(valY, valPredict)
test_cm = confusion_matrix(testY, testPredict)
print('Training Confusion Matrix:')
print(train_cm)
print('Validation Confusion Matrix:')
print(val_cm)
print('Testing Confusion Matrix:')
print(test_cm)




# print('Writing results for model {}'.format(1))
# maxtrack = int(max(df_train['new_spore_id']))
# for track in range(1, maxtrack+1):
#     cell = df_train.loc[(df_train['new_spore_id']==track)]
#     if len(cell)==0:
#         continue
#     maxslice = max(df_train.loc[(df_train['new_spore_id']==track), 'FRAME'])+1
#     minslice = min(df_train.loc[(df_train['new_spore_id']==track), 'FRAME'])
#     for sl in range(int(minslice),int(maxslice+1)):
#         x = cell.loc[(cell['FRAME']>sl-1) & (cell['FRAME']<=sl)]
#         x = x[input_cols].to_numpy()
#         if x.size > 0:
#           x=x.reshape(1, 1, 6)
#           df_train.loc[(df_train['new_spore_id']==track) & (df_train['FRAME']==sl), 'pred_germ{}'.format(1)] = model.predict(x)
 
# df_train['pred_error{}'.format(1)] = df_train['pred_germ{}'.format(1)] - df_train['GERMINATION']


# plt.errorbar(df_train['FRAME'].unique(),df_train.groupby(['FRAME']).mean()['pred_germ1'],
#             yerr = df_train.groupby(['FRAME']).sem()['pred_germ1'],color = 'b', label = 'Predicted values')
# plt.errorbar(df_train['FRAME'].unique(),df_train.groupby(['FRAME']).mean()['GERMINATION'],
#             yerr = df_train.groupby(['FRAME']).sem()['GERMINATION'] , color = 'orange', label = 'Original values')
# #plt.xticks([21,28],fontsize=20)
# #plt.xticks(np.arange(16, 29, 1.0),fontsize=20)
# #plt.xlim(16,28)
# plt.yticks(fontsize=20)
# plt.xlabel('Timestep',fontsize=14)
# plt.ylabel('Intensity',fontsize=14)
# plt.legend()


print('Writing results for model {}'.format(1))
maxtrack = int(max(df_test['new_spore_id']))
for track in range(1, maxtrack+1):
    cell = df_test.loc[(df_test['new_spore_id']==track)]
    if len(cell)==0:
        continue
    maxslice = max(df_test.loc[(df_test['new_spore_id']==track), 'FRAME'])+1
    minslice = min(df_test.loc[(df_test['new_spore_id']==track), 'FRAME'])
    for sl in range(int(minslice),int(maxslice+1)):
        x = cell.loc[(cell['FRAME']>sl-1) & (cell['FRAME']<=sl)]
        x = x[input_cols].to_numpy()
        if x.size > 0:
          x=x.reshape(1, 1, 6)
          df_test.loc[(df_test['new_spore_id']==track) & (df_test['FRAME']==sl), 'pred_germ{}'.format(1)] = model.predict(x)
 
df_test['pred_error{}'.format(1)] = df_test['pred_germ{}'.format(1)] - df_test['GERMINATION']


df_test['pred_germ1'] =(df_test['pred_germ1'] > 0.5).astype(int)

# plt.errorbar(df_test['FRAME'].unique(),df_test.groupby(['FRAME']).mean()['pred_germ1'],
#             yerr = df_test.groupby(['FRAME']).sem()['pred_germ1'],color = 'b', label = 'Predicted values')
# plt.errorbar(df_test['FRAME'].unique(),df_test.groupby(['FRAME']).mean()['GERMINATION'],
#             yerr = df_test.groupby(['FRAME']).sem()['GERMINATION'] , color = 'orange', label = 'Original values')
# #plt.xticks([21,28],fontsize=20)
# #plt.xticks(np.arange(16, 29, 1.0),fontsize=20)
# #plt.xlim(16,28)
# plt.yticks(fontsize=20)
# plt.xlabel('Timestep',fontsize=14)
# plt.ylabel('Mean spores germinated',fontsize=14)
# #plt.title('Polarity Reversal Directedness with Transfer Learning',fontsize=16)
# plt.legend()


plt.errorbar(df_test['FRAME'].unique(),df_test.groupby(['FRAME']).mean()['pred_germ1'],
            color = 'b', label = 'Predicted values')
plt.errorbar(df_test['FRAME'].unique(),df_test.groupby(['FRAME']).mean()['GERMINATION'],
             color = 'orange', label = 'Original values')
#plt.xticks([21,28],fontsize=20)
#plt.xticks(np.arange(16, 29, 1.0),fontsize=20)
#plt.xlim(16,28)
plt.yticks(fontsize=20)
plt.xlabel('Timestep',fontsize=14)
plt.ylabel('Mean spores germinated',fontsize=14)
#plt.title('Polarity Reversal Directedness with Transfer Learning',fontsize=16)
plt.legend()


import matplotlib.pyplot as plt

# Get unique spore IDs
unique_spore_ids = df_test['new_spore_id'].unique()

# Loop through each spore ID and create separate plots
for spore_id in unique_spore_ids:
    spore_data = df_test[df_test['new_spore_id'] == spore_id]
    
    # Create a new figure for each spore ID
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    
    # Plot GERMINATION values
    ax[0].errorbar(
        spore_data['FRAME'].unique(),
        spore_data.groupby(['FRAME']).mean()['GERMINATION'],
        yerr=spore_data.groupby(['FRAME']).sem()['GERMINATION'],
        color='orange',
        label='Original values'
    )
    ax[0].set_title(f'Spore ID {spore_id} - GERMINATION')
    ax[0].set_xlabel('FRAME')
    ax[0].set_ylabel('GERMINATION')
    ax[0].legend()
    
    # Plot pred_germ1 values
    ax[1].errorbar(
        spore_data['FRAME'].unique(),
        spore_data.groupby(['FRAME']).mean()['pred_germ1'],
        yerr=spore_data.groupby(['FRAME']).sem()['pred_germ1'],
        color='blue',
        label='Predicted values'
    )
    ax[1].set_title(f'Spore ID {spore_id} - Predicted GERMINATION')
    ax[1].set_xlabel('FRAME')
    ax[1].set_ylabel('pred_germ1')
    ax[1].legend()
    
    # Adjust layout
    plt.tight_layout()
    plt.show()


# Get unique spore IDs
unique_spore_ids = df_test['SPORE_ID'].unique()

# Define the number of rows and columns for the grid
n_rows = 18
n_cols = 3

# Create a figure with a grid of subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 30))

for idx, spore_id in enumerate(unique_spore_ids):
    row = idx // n_cols
    col = idx % n_cols
    spore_data = df_test[df_test['SPORE_ID'] == spore_id]
    
    # Plot GERMINATION values
    axes[row, col].errorbar(
        spore_data['FRAME'].unique(),
        spore_data.groupby(['FRAME']).mean()['GERMINATION'],
        yerr=spore_data.groupby(['FRAME']).sem()['GERMINATION'],
        color='orange',
        label='Original values'
    )
    
    # Plot pred_germ1 values
    axes[row, col].errorbar(
        spore_data['FRAME'].unique(),
        spore_data.groupby(['FRAME']).mean()['pred_germ1'],
        yerr=spore_data.groupby(['FRAME']).sem()['pred_germ1'],
        color='blue',
        label='Predicted values'
    )
    
    # Set title and labels
    axes[row, col].set_title(f'Spore ID {spore_id}')
    axes[row, col].set_xlabel('FRAME')
    axes[row, col].set_ylabel('Values')
    axes[row, col].legend()

# Adjust layout
plt.tight_layout()
plt.show()


# Number of subplots per figure
n_rows, n_cols = 3, 3

# Number of figures needed
n_figures = len(unique_spore_ids) // (n_rows * n_cols) + (len(unique_spore_ids) % (n_rows * n_cols) != 0)

for fig_num in range(n_figures):
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 15))
    axes = axes.flatten()  # Flatten the 2D array of axes to make indexing easier
    
    for i in range(n_rows * n_cols):
        spore_idx = fig_num * (n_rows * n_cols) + i
        if spore_idx >= len(unique_spore_ids):
            break

        spore_id = unique_spore_ids[spore_idx]
        spore_data = df_test[df_test['SPORE_ID'] == spore_id]
        
        # Plot GERMINATION values
        axes[i].errorbar(
            spore_data['FRAME'].unique(),
            spore_data.groupby(['FRAME']).mean()['GERMINATION'],
            yerr=spore_data.groupby(['FRAME']).sem()['GERMINATION'],
            color='orange',
            label='Original values'
        )
        
        # Plot pred_germ1 values
        axes[i].errorbar(
            spore_data['FRAME'].unique(),
            spore_data.groupby(['FRAME']).mean()['pred_germ1'],
            yerr=spore_data.groupby(['FRAME']).sem()['pred_germ1'],
            color='blue',
            label='Predicted values'
        )
        
        # Set title and labels
        axes[i].set_title(f'Spore ID {spore_id}')
        axes[i].set_xlabel('FRAME')
        axes[i].set_ylabel('Values')
        axes[i].legend()

    # Adjust layout
    plt.tight_layout()
    plt.show()

###########################################################################################################
#### Interpretation### And Sindy
# Prepare the data
X_train = trainX.reshape((trainX.shape[0], -1))
y_train = trainPredict.reshape((trainPredict.shape[0], -1))
X_val = valX.reshape((valX.shape[0], -1))
y_val = valPredict.reshape((valPredict.shape[0], -1))
X_test = testX.reshape((testX.shape[0], -1))
y_test = testPredict.reshape((testPredict.shape[0], -1))


# Combine training and validation data for SINDy
X = np.concatenate((X_train, X_val), axis=0)
y = np.concatenate((y_train, y_val), axis=0)

import pysindy as ps

library = ps.PolynomialLibrary(degree=2)
optimizer = ps.STLSQ(threshold=0.1)
sindy_model = ps.SINDy(feature_library=library, optimizer=optimizer)


# Fit the SINDy model
sindy_model.fit(X, t=5)  # Assuming time step t=1
sindy_model.print()


# Predict on test data
test_predictions = sindy_model.predict(X_test)

# Plotting
plt.figure()
plt.plot(y_test, label='True values')
plt.plot(test_predictions)
plt.xlabel('Time')
plt.ylabel('Value')
plt.legend()
plt.show()