import pytest 
import sys
sys.path.append('../')
from DNN.rcnn_sat import b_d_net_readout_0, b_d_net_readout_1, b_d_net_readout_2, b_d_net_readout_3, b_d_net_readout_4, b_d_net_readout_5, b_d_net_readout_final
from utils import non_trainable
import tensorflow as tf
import numpy as np
@pytest.mark.parametrize("input_shape", [(128, 128,3)])
@pytest.mark.parametrize("num_readouts", [7])
@pytest.mark.parametrize("num_classes", [2])
@pytest.mark.parametrize("num_images_all",[10])
def test_readout(input_shape, num_readouts, num_classes,num_images_all):
    all_models = [b_d_net_readout_0, b_d_net_readout_1, b_d_net_readout_2, b_d_net_readout_3, b_d_net_readout_4, b_d_net_readout_5, b_d_net_readout_final]
    input_layer = tf.keras.layers.Input(input_shape)
    predictions_all_readouts = np.zeros([num_readouts,num_images_all,num_classes])
    predictions_all_readouts[:] = np.nan

    for m,readout_model in enumerate(all_models):
        tf.keras.backend.clear_session() 
        model = readout_model(input_layer,classes_scenes=2)
        #define trainable layers (Readout)
        layer = non_trainable(model)
        if 'Readout' in layer.name:
            assert layer.trainable == True
        else:
            assert layer.trainable == False   


             
        
        
    
    
    

