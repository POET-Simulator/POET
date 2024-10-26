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

def prediction_step(model, x, batch_size, cluster_labels):
    # TODO PREDICT ACCORDING TO CLUSTER

    prediction = model.predict(x, batch_size)
    return np.array(prediction, dtype=np.float64)


def get_weights(model):
    weights = model.get_weights()
    return weights

def training_step(model, x, y, batch_size, epochs, 
                  train_cluster, output_file_path):
    # Check clustering of input data
    # and only train for the cluster where nothing is happening
    
    # TODO TRAIN ACCORDING TO CLUSTER
    print("Training cluster: " + str(train_cluster), flush=True)
    print("Training data size is: " + str(len(x)), flush=True)

    history = model.fit(x, y, 
                        epochs=epochs,
                        batch_size=batch_size)
    if output_file_path:
        if output_file_path[-6:] != ".keras":
            output_file_path += ".keras"
        model.save(output_file_path)
    return history
