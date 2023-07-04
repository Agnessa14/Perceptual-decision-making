import tensorflow as tf

def non_trainable(model: tf.keras.Model):
   for layer in model.layers:
            if 'Readout' in layer.name:
                layer.trainable = True
            else:
                layer.trainable = False
            return layer    