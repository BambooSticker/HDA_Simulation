import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pylab as plt
from pyomo.environ import *
from idaes.core.util.model_statistics import degrees_of_freedom
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from omlt import OmltBlock, OffsetScaling
from omlt.neuralnet import FullSpaceNNFormulation, NetworkDefinition
from idaes.core.surrogate.sampling.scaling import OffsetScaler
from omlt.io import load_keras_sequential
from time import time
from IPython.display import display 
from idaes.core import FlowsheetBlock
from tensorflow.keras.models import load_model

def load_surrogate(b): 
    model = load_model('./surrogate_model/F102_Surrogate.h5')
    train_data = pd.read_csv("./data/train_outcome.csv",index_col = 0)
    
    input_data = train_data.iloc[:, :3]
    output_data = train_data.iloc[:, 3:5]
    
    input_labels = input_data.columns
    output_labels = output_data.columns
    
    x = input_data
    y = output_data
    
    input_scaler = OffsetScaler.create_normalizing_scaler(x)
    output_scaler = OffsetScaler.create_normalizing_scaler(y)
    
    offset_inputs = input_scaler.offset_series().to_numpy()
    factor_inputs = input_scaler.factor_series().to_numpy()
    offset_outputs = output_scaler.offset_series().to_numpy()
    factor_outputs = output_scaler.factor_series().to_numpy()
    
    m.F101.VLE = OmltBlock()
    
    #scale_y = (-0.25, 0.125) #(mean,stdev) of the output
    scaler = OffsetScaling(offset_inputs=offset_inputs,
                        factor_inputs=factor_inputs,
                        offset_outputs=offset_outputs,
                        factor_outputs=factor_outputs)
    
    scaled_input_bounds = {0: (0, 1),
                           1: (0, 1),
                           2: (0, 1)
                          }
    
    net = load_keras_sequential(model, scaler, scaled_input_bounds)
    formulation = FullSpaceNNFormulation(net)
    m.F101.VLE.build_formulation(formulation)
    
    @b.Constraint()
    def connect_Temperature(b):
        return b.outlet.T == b.VLE.inputs[0]
    @b.Constraint()
    def connect_Pressure(b):
        return b.outlet.P == b.VLE.inputs[1]
    @b.Constraint()
    def connect_Fraction_B(b):
        return b.inlet.z['benzene'] == b.VLE.inputs[2]
    
    @b.Constraint()
    def prediction_Liq_B(b):
        return b.outlet.x['benzene'] == b.VLE.outputs[0]
    @b.Constraint()
    def prediction_Vap_B(b):
        return b.outlet.y['benzene'] == b.VLE.outputs[1]