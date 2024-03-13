import pandas as pd
from sumonet.model.architecture import SUMOnet
from sumonet.utils.encodings import Encoding
from sumonet.utils.data_pipe import Data


#* The create_dataframe function creates a pandas DataFrame from the protein_id, protein_seq, k_position, predicted_probs, and predicted_labels.
def create_dataframe(protein_id, protein_seq, k_position, predicted_probs, predicted_labels):

    data_dict = {'protein_id':protein_id,
                 'protein_seq':protein_seq,
                 'lysine_position':k_position,
                 'nonsumoylation_class_probs':predicted_probs[:,0],
                 'sumoylation_class_probs':predicted_probs[:,1],
                 'predicted_labels':predicted_labels}

    return pd.DataFrame(data_dict)

#* The prediction_outputs function takes the protein_id, protein_seq, k_position, and predicted_probs as input and returns the DataFrame created by the create_dataframe function.
def prediction_outputs(protein_id, protein_seq, k_position, predicted_probs):

    
    predicted_labels = predicted_probs.argmax(-1)
    return create_dataframe(protein_id, protein_seq, k_position, predicted_probs, predicted_labels)

#* The load_models function loads the SUMOnet model.
def load_models():

    my_model = SUMOnet()
    my_model.load_weights()
    return my_model

#* The make_prediction function takes the protein_ids, protein_seqs, and k_positions as input and returns the DataFrame created by the prediction_outputs function.
def make_prediction(protein_ids, protein_seqs, k_positions):
    encoder = Encoding()
    
    X_train = encoder.encode_data(protein_seqs)
    
    
    
    my_model = load_models()
   
    predicted_probs = my_model.predict(X_train)
    df = prediction_outputs(protein_ids, protein_seqs, k_positions, predicted_probs)
    return df
    

#* The seqIOParser function takes the record as input and returns the protein_ids, protein_seqs, and k_positions.
def seqIOParser(record):

    data_processes = Data()
    protein_ids, protein_seqs, k_positions = [], [], []
    
     
    mers, k_position = data_processes.find_mers_with_K(str(record.seq))
        
    protein_seqs += mers
    k_positions += k_position
    protein_ids += [record.id] * len(mers)
        
    return protein_ids, protein_seqs, k_positions 




