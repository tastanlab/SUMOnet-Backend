from django.http import JsonResponse
from sumonetML.model.architecture import SUMOnet
from sumonetML.utils.encodings import Encoding
from sumonetML.utils.data_pipe import Data
import numpy as np
def hello(request):
    
    
    data = Data()
    X_train, y_train, X_test, y_test = data.load_sumonet_experiment_data()
    
    y_train = np.asarray(y_train)
    y_train = (y_train[:,None] == np.arange(2)).astype(int)
    SUMOnet3_model = SUMOnet()
    SUMOnet3_model.load_weights()
    
    encoder = Encoding(encoderType='blosum62') ## Firstly we need to encode our test data
    X_test_encoded = encoder.encode_data(X_test)
    y_preds = SUMOnet3_model.predict(X_test_encoded)
    
    return JsonResponse({'message': y_preds.tolist()})

