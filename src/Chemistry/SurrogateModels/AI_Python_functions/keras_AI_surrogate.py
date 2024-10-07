import tensorflow as tf
import os

def initiate_model(model_file_path):
    print(model_file_path, flush=True)
    #model = tf.keras.load_model(model_file_path)
    #return model
    return

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
    return prediction