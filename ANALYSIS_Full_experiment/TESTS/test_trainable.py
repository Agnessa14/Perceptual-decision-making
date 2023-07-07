import pytest 
import sys
sys.path.append('../')
from DNN.rcnn_sat import b_d_net_readout_0, b_d_net_readout_1, b_d_net_readout_2, b_d_net_readout_3, b_d_net_readout_4, b_d_net_readout_5, b_d_net_readout_final
from utils import non_trainable, train_readout, assertion_block
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

    for readout_model in all_models:
        tf.keras.backend.clear_session() 
        model = readout_model(input_layer,classes_scenes=2)
       
        non_trainable(model)
        
                
            
@pytest.mark.parametrize("x_train_shape", [(1000,128, 128,3)])
@pytest.mark.parametrize("y_train_shape", [(1000,1)])
@pytest.mark.parametrize("x_val_shape", [(10,128, 128,3)])
@pytest.mark.parametrize("y_val_shape", [(10,1)])
@pytest.mark.parametrize("classes_scenes", [(2)])
def test_train(x_train_shape,y_train_shape,x_val_shape,y_val_shape, classes_scenes):
    all_models = [b_d_net_readout_0, b_d_net_readout_1, b_d_net_readout_2, b_d_net_readout_3, b_d_net_readout_4, b_d_net_readout_5, b_d_net_readout_final]
    input_shape = (x_train_shape[1],x_train_shape[2],x_train_shape[3])
    x_train = tf.random.normal(x_train_shape)
    y_train = tf.random.uniform(y_train_shape, maxval=classes_scenes, dtype = tf.dtypes.int32)
    x_val = tf.random.normal(x_val_shape)
    y_val = tf.random.uniform(y_val_shape, maxval = classes_scenes, dtype=tf.dtypes.int32)
    batch_size = y_val_shape[0]
    input_layer = tf.keras.layers.Input(input_shape)
    for readout_model in all_models:
        tf.keras.backend.clear_session() 
        
        model = readout_model(input_layer,classes_scenes=classes_scenes)
        non_trainable(model)
        
        train_readout(x_train, y_train, x_val, y_val,batch_size, model)
    

if __name__ == "__main__":
    test_train((1000,128, 128,3),(1000,1),(10,128, 128,3),(10,1), 2)             
        
        
    
    
    

