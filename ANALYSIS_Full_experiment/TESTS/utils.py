import tensorflow as tf

def non_trainable(model: tf.keras.Model):
   for layer in model.layers:
            if 'Readout' in layer.name:
                layer._trainable = True
            else:
                layer._trainable = False
           
        
def assertion_block(model: tf.keras.Model):
    for layer in model.layers:  
            if 'Readout' in layer.name:
                assert layer.trainable == True
            else:
                assert layer.trainable == False           
def train_readout(x_train, y_train, x_val, y_val, batch_size, model):
    base_learning_rate = 0.0001
    assertion_block(model)
    model.compile(optimizer=tf.keras.optimizers.Adam(lr=base_learning_rate),
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(),
                  metrics=['accuracy'])
    print(model.summary())
    
    history = model.fit(
        x_train,
        y_train,
        shuffle=True,
        batch_size=batch_size,
        epochs=1, 
        validation_data=(x_val,y_val)
    )        
    
    #return history