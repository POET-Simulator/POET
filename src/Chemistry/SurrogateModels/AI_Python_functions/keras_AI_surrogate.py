import tensorflow as tf
import numpy as np
from sklearn.cluster import KMeans
import os

os.environ["TF_XLA_FLAGS"] = "--tf_xla_cpu_global_jit"
os.environ["XLA_FLAGS"] = "--xla_gpu_cuda_data_dir=" + cuda_dir

def k_means(data, k=2, tol=1e-6):
    kmeans = KMeans(n_clusters=k, tol=tol)
    labels = kmeans.fit_predict(data)
    return labels

def initiate_model(model_file_path):
    print("AI: Model loaded from: " + model_file_path, flush=True)
    model = tf.keras.models.load_model(model_file_path)
    return model

def prediction_step(model, model_reactive, x, cluster_labels, batch_size)
    # Catch input size mismatches
    model_input_shape = model.input_shape[1:]
    if model_input_shape != x.shape[1:]:
        print(f"Input data size {x.shape[1:]} does not match model input size {model_input_shape}",
              flush=True)        

    # Predict separately if clustering is used
    if cluster_labels:
        cluster_labels = np.asarray(cluster_labels, dtype=bool)
        # Combine results
        prediction = np.zeros_like(x)
        prediction[cluster_labels] = model_reactive.predict(x[cluster_labels], batch_size)
        prediction[~cluster_labels] = model.predict(x[~cluster_labels], batch_size)
    else:
        prediction = model.predict(x, batch_size)
    return np.array(prediction, dtype=np.float64)


def get_weights(model):
    weights = model.get_weights()
    return weights

def training_step(model, x, y, batch_size, epochs, 
                  output_file_path):
    history = model.fit(x, y, 
                        epochs=epochs,
                        batch_size=batch_size)
    if output_file_path:
        model.save(output_file_path)
    return history
