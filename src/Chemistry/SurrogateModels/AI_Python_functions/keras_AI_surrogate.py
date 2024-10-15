import tensorflow as tf
import numpy as np
import os

def initiate_model(model_file_path):
    model = tf.keras.models.load_model(model_file_path)
    return model
    
def training_step(model, x, y, x_val, y_val, batch_size, epochs):
    history = model.fit(x, y, 
                        epochs=epochs,
                        batch_size=batch_size,
                        validation_data=(x_val, y_val))
    print(history, flush=True)
    return history["val_loss"]

def prediction_step(model, x, batch_size):
    prediction = model.predict(x, batch_size)
    return np.array(prediction, dtype=np.float64)
    #with tf.device('/cpu:0'):
    #    prediction = model.predict(x, batch_size)
    #    return np.array(prediction, dtype=np.float64)

def get_weights(model):
    weights = model.get_weights()
    return weights

def training_step(model, x, y, batch_size, epochs, output_file_path):
    history = model.fit(x, y, 
                        epochs=epochs,
                        batch_size=batch_size)
    if output_file_path:
        if output_file_path[-6:] != ".keras":
            output_file_path += ".keras"
        model.save(output_file_path)
    return history
    