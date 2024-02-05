from django.http import JsonResponse
from sumonetML.model.architecture import SUMOnet
from sumonetML.utils.encodings import Encoding
from sumonetML.utils.data_pipe import Data
from rest_framework.decorators import api_view
from sumonet.serializers import UniprotSerializer, ProteinSequenceSerializer
from rest_framework import status
from rest_framework.response import Response
import re
import json
import pandas as pd

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
    


def validateProteinSequence(seq, alphabet='protein'):
    
    alphabets = {'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}


    if alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False


@api_view(['POST'])
def uniprotPrediction(request):
    
    serializer = UniprotSerializer(data=request.data)
    
    if serializer.is_valid():
        data_processes = Data()
        uniprot_id = serializer.data['uniprot_id']            
        lysine_position = serializer.data['lysine_position']
        protein_seq = data_processes.retrive_protein_sequence_with_uniprotid(uniprot_id)
            
        if protein_seq == None:
            return Response({'error': 'No sequence found with this UniprotID.'}, status=status.HTTP_400_BAD_REQUEST)
        
        if uniprot_id == '' or uniprot_id == None:
            return Response({'error': 'UniprotID must be entered.'}, status=status.HTTP_400_BAD_REQUEST)
                
        if lysine_position:
            lysine_position = int(lysine_position)
            protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq,uniprot_id,lysine_position)
        else:
            protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq,uniprot_id)

        df = make_prediction(protein_ids, protein_seqs, k_positions)
            
        df = df.head()
            
        return Response(df.to_dict(), status=status.HTTP_200_OK)
            
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        
        
@api_view(['POST'])
def proteinSequence(request):
    serializer = ProteinSequenceSerializer(data=request.data)
    
    if serializer.is_valid():
        data_processes = Data()
        protein_seq = serializer.data['protein_seq']
        
        if protein_seq == '' or protein_seq == None or protein_seq == []:
            return Response({'error': 'Protein sequence must be entered.'}, status=status.HTTP_400_BAD_REQUEST)
        
        jsonList = []
        uniprotId_pattern = re.compile(r">sp\|([A-Z0-9]+)\|")
        
        for i in range(len(protein_seq)):
            
            if validateProteinSequence(protein_seq[i]) == False:
                return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
            
            protein_ids, protein_seqs, k_positions = data_processes.protein_sequence_input(protein_seq[i].split())
        
            df = make_prediction(protein_ids, protein_seqs, k_positions)
        
            match = uniprotId_pattern.search(protein_seq[i])

            # Extract the UniProt ID from the match
            if match:
                uniprot_id = match.group(1)
        
            result_list = [
                {
                    "protein_id": uniprot_id,
                    "peptide_seq": protein_seq,
                    "lysine_position": lysine_position,
                    "nonsumoylation_class_probs": nonsumoylation_class_probs,
                    "sumoylation_class_probs": sumoylation_class_probs,
                    "predicted_labels": predicted_labels
                }
                for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels'])
            ]
            
            jsonList.append(result_list)

        # Convert the list to JSON without escaping
        #json_result = json.dumps(result_list, indent=2)
        
        return Response(jsonList, content_type='application/json', status=status.HTTP_200_OK)
    
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
 
        