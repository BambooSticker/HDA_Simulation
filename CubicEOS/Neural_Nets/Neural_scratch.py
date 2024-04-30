import numpy as np

# Define the activation function and its derivative
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def sigmoid_derivative(x):
    return x * (1 - x)

# Define the structure of the neural network
inputLayer_neurons = 2  # number of features in data set
hiddenLayer_neurons = 3  # number of hidden layers neurons
output_neurons = 1  # number of neurons at output layer

# Randomly initialize the weight and biases
weights_input_hidden = np.random.uniform(size=(inputLayer_neurons, hiddenLayer_neurons))
bias_input_hidden = np.random.uniform(size=(1, hiddenLayer_neurons))
weights_hidden_output = np.random.uniform(size=(hiddenLayer_neurons, output_neurons))
bias_hidden_output = np.random.uniform(size=(1, output_neurons))

# Define the training data
X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])  # training inputs
y = np.array([[0], [1], [1], [0]])  # training outputs

epochs = 12000  # number of training iterations
learning_rate = 0.1  # learning rate

for _ in range(epochs):

    # Forward Propagation
    hiddenLayer_linearTransform = np.dot(X, weights_input_hidden) + bias_input_hidden
    hiddenLayer_activations = sigmoid(hiddenLayer_linearTransform)

    outputLayer_linearTransform = np.dot(hiddenLayer_activations, weights_hidden_output) + bias_hidden_output
    output = sigmoid(outputLayer_linearTransform)

    # Backward Propagation
    error = y - output
    d_output = error * sigmoid_derivative(output)

    error_hidden_layer = d_output.dot(weights_hidden_output.T)
    d_hiddenLayer = error_hidden_layer * sigmoid_derivative(hiddenLayer_activations)

    # Updating Weights and Biases
    weights_hidden_output += hiddenLayer_activations.T.dot(d_output) * learning_rate
    bias_hidden_output += np.sum(d_output, axis=0, keepdims=True) * learning_rate
    weights_input_hidden += X.T.dot(d_hiddenLayer) * learning_rate
    bias_input_hidden += np.sum(d_hiddenLayer, axis=0, keepdims=True) * learning_rate

print(output)