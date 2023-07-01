'''
Keras implementation of B networks
'''

import tensorflow as tf


def b_layer(x, filters, kernel, layer_num, pooling=True):
    '''Base layer for B models
    '''
    if pooling:
        x = tf.keras.layers.MaxPool2D(
            pool_size=(2, 2),
            name='MaxPool_Layer_{}'.format(layer_num))(x)

    x = tf.keras.layers.Conv2D(
        filters, kernel, padding='same', use_bias=False,
        kernel_initializer='glorot_uniform',
        kernel_regularizer=tf.keras.regularizers.l2(1e-6),
        name='Conv_Layer_{}'.format(layer_num))(x)

    data_format = tf.keras.backend.image_data_format()
    norm_axis = -1 if data_format == 'channels_last' else -3
    x = tf.keras.layers.BatchNormalization(
        norm_axis,
        name='BatchNorm_Layer_{}'.format(layer_num))(x)

    x = tf.keras.layers.Activation(
        'relu', name='ReLU_Layer_{}'.format(layer_num))(x)

    return x


def readout(x, classes_scenes):
    '''Readout layer
    '''
    x = tf.keras.layers.GlobalAvgPool2D(name='GlobalAvgPool')(x)
    # x = tf.keras.layers.Dense(
    #     classes, kernel_initializer='glorot_uniform',
    #     kernel_regularizer=tf.keras.regularizers.l2(1e-6),
    #     name='ReadoutDense')(x)
    x = tf.keras.layers.Dense(
    classes_scenes, kernel_initializer='glorot_uniform',
    kernel_regularizer=tf.keras.regularizers.l2(1e-6),
    bias_initializer='zeros',
    name='ReadoutDense_Scenes')(x)
    x = tf.keras.layers.Activation('softmax', name='Softmax')(x)
    return x


def intermediate_readout(x, classes_scenes, layer_num):
    '''Readout layer
    '''
    x = tf.keras.layers.GlobalAvgPool2D(name='GlobalAvgPool')(x)
    x = tf.keras.layers.Dense(
    classes_scenes, kernel_initializer='glorot_uniform',
    kernel_regularizer=tf.keras.regularizers.l2(1e-6),
    bias_initializer='zeros',
    name='Intermediate_ReadoutDense_Scenes_{}'.format(layer_num))(x)
    x = tf.keras.layers.Activation('softmax', name='Softmax_intermediate_{}'.format(layer_num))(x)
    return x



def b_d_net_readouts(input_tensor, classes_scenes):
    '''Defines a B-D model with readouts after every pool layer
    '''

    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    readout_0 = intermediate_readout(x,classes_scenes,0)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    readout_1 = intermediate_readout(x,classes_scenes,1)
    x = b_layer(x, 192, 3, 5, pooling=False)
    x = b_layer(x, 256, 3, 6)
    readout_2 = intermediate_readout(x,classes_scenes,2)
    x = b_layer(x, 256, 3, 7, pooling=False)
    x = b_layer(x, 512, 3, 8)
    readout_3 = intermediate_readout(x,classes_scenes,3)
    x = b_layer(x, 512, 3, 9, pooling=False)
    x = b_layer(x, 1024, 3, 10)
    readout_4 = intermediate_readout(x,classes_scenes,4)
    x = b_layer(x, 1024, 3, 11, pooling=False)
    x = b_layer(x, 2048, 1, 12) 
    readout_5 = intermediate_readout(x,classes_scenes,5)
    x = b_layer(x, 2048, 1, 13, pooling=False)
    output_tensor = readout(x, classes_scenes)

    return tf.keras.Model(inputs=input_tensor, outputs=(readout_0,readout_1,readout_2,readout_3,readout_4,readout_5,output_tensor)) #not sure

def b_d_net_readout_0(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 0
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    readout_0 = intermediate_readout(x,classes_scenes,0)    
    x = b_layer(x, 128, 5, 2)
    # readout_0 = intermediate_readout(x,classes_scenes,0)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_0)

def b_d_net_readout_1(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 1
    '''
   
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    readout_1 = intermediate_readout(x,classes_scenes,1)
    # x = b_layer(x, 192, 3, 4)
    # readout_1 = intermediate_readout(x,classes_scenes,1)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_1)

def b_d_net_readout_2(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 2
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    x = b_layer(x, 192, 3, 5, pooling=False)
    readout_2 = intermediate_readout(x,classes_scenes,2)
    # x = b_layer(x, 256, 3, 6)
    # readout_2 = intermediate_readout(x,classes_scenes,2)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_2)

def b_d_net_readout_3(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 3
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    x = b_layer(x, 192, 3, 5, pooling=False)
    x = b_layer(x, 256, 3, 6)
    x = b_layer(x, 256, 3, 7, pooling=False)
    readout_3 = intermediate_readout(x,classes_scenes,3)
    # x = b_layer(x, 512, 3, 8)
    # readout_3 = intermediate_readout(x,classes_scenes,3)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_3)

def b_d_net_readout_4(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 4
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    x = b_layer(x, 192, 3, 5, pooling=False)
    x = b_layer(x, 256, 3, 6)
    x = b_layer(x, 256, 3, 7, pooling=False)
    x = b_layer(x, 512, 3, 8)
    x = b_layer(x, 512, 3, 9, pooling=False)
    readout_4 = intermediate_readout(x,classes_scenes,4)    
    # x = b_layer(x, 1024, 3, 10)
    # readout_4 = intermediate_readout(x,classes_scenes,4)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_4)

def b_d_net_readout_5(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after pool layer 5
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    x = b_layer(x, 192, 3, 5, pooling=False)
    x = b_layer(x, 256, 3, 6)
    x = b_layer(x, 256, 3, 7, pooling=False)
    x = b_layer(x, 512, 3, 8)
    x = b_layer(x, 512, 3, 9, pooling=False)
    x = b_layer(x, 1024, 3, 10)
    x = b_layer(x, 1024, 3, 11, pooling=False)
    readout_5 = intermediate_readout(x,classes_scenes,5)
    # x = b_layer(x, 2048, 1, 12) 
    # readout_5 = intermediate_readout(x,classes_scenes,5)

    return tf.keras.Model(inputs=input_tensor, outputs=readout_5)

def b_d_net_readout_final(input_tensor, classes_scenes):
    '''Defines a B-D model with readout after last layer
    '''
    
    x = b_layer(input_tensor, 96, 7, 0, pooling=False)
    x = b_layer(x, 96, 7, 1, pooling=False)
    x = b_layer(x, 128, 5, 2)
    x = b_layer(x, 128, 5, 3, pooling=False)
    x = b_layer(x, 192, 3, 4)
    x = b_layer(x, 192, 3, 5, pooling=False)
    x = b_layer(x, 256, 3, 6)
    x = b_layer(x, 256, 3, 7, pooling=False)
    x = b_layer(x, 512, 3, 8)
    x = b_layer(x, 512, 3, 9, pooling=False)
    x = b_layer(x, 1024, 3, 10)
    x = b_layer(x, 1024, 3, 11, pooling=False)
    x = b_layer(x, 2048, 1, 12) 
    x = b_layer(x, 2048, 1, 13, pooling=False)
    output_tensor = readout(x, classes_scenes)
    
    return tf.keras.Model(inputs=input_tensor, outputs=output_tensor)