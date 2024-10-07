import tensorflow as tf
import numpy as np
import os

def initiate_model(model_file_path):
    model = tf.keras.models.load_model(model_file_path)
    return model
    
def training_step(model, x, y, x_val, y_val, batch_size, epochs):
    epochs = 2000 # This is a constant parameter during all experiments
    history = model.fit(x, y, 
                        epochs=epochs,
                        batch_size=batch_size,
                        validation_data=(x_val, y_val))
    print(history, flush=True)
    return history["val_loss"]

def prediction_step(model, x, batch_size):
    prediction = model.predict(x, batch_size)
    print("Prediction from Python", flush=True)
    print(prediction, flush=True)
    return np.array(prediction, dtype=np.float64)

def get_weights(model):
    weights = model.get_weights()
    # Convert all weight and bias arrays to float64 (double precision)
    weights_as_double = [w.astype(np.float64) for w in weights]
    return weights_as_double