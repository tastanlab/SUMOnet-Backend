from django.http import JsonResponse, StreamingHttpResponse
from sumonetML.model.architecture import SUMOnet
from sumonetML.utils.encodings import Encoding
from sumonetML.utils.data_pipe import Data
from rest_framework.decorators import api_view, parser_classes
from sumonet.serializers import UniprotSerializer, ProteinSequenceSerializer
from rest_framework import status
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser
import re
import json
from io import StringIO
import pandas as pd
from Bio import SeqIO

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
    
    #alphabets = {'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}

                
    if seq[0:3] == ">sp": #alphabets[alphabet].search(seq) is not None:
         return True
    else:
         return False

def parse_sequence(seq):
    return ">" + seq.strip()  # Strip any leading/trailing whitespace

def parse_fasta_to_protein_seq(fasta):
    
    # Remove all '\n' and '\r' characters from the fasta content
    fasta_content_cleaned = fasta.replace('\n', '').replace('\r', '')
 
  
    # Regular expression pattern
    SV_pattern = r"(SV=\d+)" #birden fazla space olma durumu handle edilebilir.

    # Find all occurrences of the pattern
    for match in re.finditer(SV_pattern, fasta_content_cleaned):
        match_str = match.group(0)
        # Check if the match is followed by a space
        if match.end() < len(fasta_content_cleaned) and fasta_content_cleaned[match.end()] != ' ':
            # If not, add a space after the match
            fasta_content_cleaned = fasta_content_cleaned[:match.end()] + ' ' + fasta_content_cleaned[match.end():]
        #print("Already spaced:", fasta_content_cleaned)

    # Split the content based on ">"
    sequences = fasta_content_cleaned.split(">")[1:]
    
    # Store each sequence in a list
    parsed_sequences = []
    
    # Parse each sequence and store in the list
    for seq in sequences:
        parsed_sequences.append(parse_sequence(seq))
    
    # Return the result as a dictionary
    #print(parsed_sequences)
    return {"protein_seq": parsed_sequences}

def parse_protein_sequence(protein_sequence):
    
    #print(json_object)
    
    # Remove all '\n' and '\r' characters from the fasta content
    protein_sequence_cleaned = protein_sequence.replace('\n', '').replace('\r', '')
 
  
    # Regular expression pattern
    SV_pattern = r"(SV=\d+)" #birden fazla space olma durumu handle edilebilir.

    # Find the first occurrence of the pattern
    match = re.search(SV_pattern, protein_sequence_cleaned)

    if match:
        # Check if the match is followed by a space
        if match.end() < len(protein_sequence_cleaned) and protein_sequence_cleaned[match.end()] != ' ':
            # If not, add a space after the match
            modified_string = protein_sequence_cleaned[:match.end()] + ' ' + protein_sequence_cleaned[match.end():]
        else:
            modified_string = protein_sequence_cleaned
    else:
        modified_string = protein_sequence_cleaned

    return modified_string

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
            
        result = [
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
            
        return Response(result, status=status.HTTP_200_OK)
            
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        

@api_view(['POST'])
def proteinSequence(request):
    #print(request.data)
    protein_seq = request.data.get('protein_seq', '')
    #print(protein_seq)
    result = parse_fasta_to_protein_seq(protein_seq) # fonksiyonun adı değismeli
    print(result)

    serializer = ProteinSequenceSerializer(data=result)
    
    if serializer.is_valid():
        protein_sequences = serializer.validated_data['protein_seq']
        

        if not protein_sequences:
            return Response({'error': 'No protein sequences provided.'}, status=status.HTTP_400_BAD_REQUEST)

        if len(protein_sequences) > 5: # burayı update etmek gerek
            return Response({'error': 'You can enter up to five proteins. For entries with more than five proteins, use the FASTA file option.'}, status=status.HTTP_400_BAD_REQUEST)
        
        data_processes = Data()
        jsonList = []
        uniprotId_pattern = re.compile(r">sp\|([A-Z0-9]+)\|")
        
        for protein_seq in protein_sequences:
            if not validateProteinSequence(protein_seq, alphabet='protein'):
                return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
            
            protein_sequence = parse_protein_sequence(protein_seq)
            protein_ids, protein_seqs, k_positions = data_processes.protein_sequence_input(protein_sequence.split())
            df = make_prediction(protein_ids, protein_seqs, k_positions)
        
            match = uniprotId_pattern.search(protein_sequence)

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

        return Response(jsonList, content_type='application/json', status=status.HTTP_200_OK)
    
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
 

@api_view(['POST'])
@parser_classes([MultiPartParser])
def fastaFile(request):
    file_obj = request.FILES.get('file')
    if file_obj:
        file_name, file_extension = file_obj.name.split('.')

        if str(file_extension.lower()) != 'fasta':
            return Response({'error': 'Invalid file extension. Only fasta files are allowed.'}, status=status.HTTP_400_BAD_REQUEST)

        print("here")
        records = file_obj.read().decode('utf-8')

        # You can use a custom filename, or keep the original filename
        custom_filename = f'{file_name}.fasta'
        content_disposition = f'attachment; filename={custom_filename}'

        response = Response({'fasta': records}, status=status.HTTP_200_OK)
        response['Content-Disposition'] = content_disposition
        # Parse fasta content to protein sequences
        print("RECORDS: ", records)
        result = parse_fasta_to_protein_seq(records)
        # print("RECORDS: ", records)
        print("RESULT", result)
        serializer = ProteinSequenceSerializer(data=result)  # Verinin parse edilmiş kısmı direkt sokuluyor.

        if serializer.is_valid():
            data_processes = Data()
            protein_seq = serializer.data['protein_seq']
            protein_seq_len = len(protein_seq)
            # print(protein_seq)

            if protein_seq == '' or protein_seq == None or protein_seq == []:
                return Response({'error': 'Protein sequence must be entered.'}, status=status.HTTP_400_BAD_REQUEST)

            jsonList = []
            uniprotId_pattern = re.compile(r">sp\|([A-Z0-9]+)\|")

            for i in range(len(protein_seq)):
                # ! Validate protein sequence
                if validateProteinSequence(protein_seq[i], alphabet='protein') == False:  # it simply checks the >sp part.
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
                    for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs,
                    predicted_labels in
                    zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'],
                        df['sumoylation_class_probs'], df['predicted_labels'])
                ]

                jsonList.append(result_list)

            return Response(jsonList, content_type='application/json', status=status.HTTP_200_OK)

        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
    
    else:
        return Response({"message": "An error occured."}, status=status.HTTP_400_BAD_REQUEST)

 

def seqIOParser(record):
    data_processes = Data()
    protein_ids, protein_seqs, k_positions = [], [], []
    
     
    mers, k_position = data_processes.find_mers_with_K(str(record.seq))
        
    protein_seqs += mers
    k_positions += k_position
    protein_ids += [record.id] * len(mers)
        
    return protein_ids, protein_seqs, k_positions    
        
        
        
@api_view(['POST'])
def denemeProteinSeq(request):
    protein_seq = request.data.get('protein_seq', '')

    
    fasta_file = StringIO(protein_seq)

    data_processes = Data()
    jsonList = []
    idS, seqS, positionS = [], [], []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_ids, protein_seqs, k_positions = seqIOParser(record)
        idS.extend(protein_ids)
        seqS.extend(protein_seqs)
        positionS.extend(k_positions)
        
          
        
    print(idS)
    df = make_prediction(idS, seqS, positionS)  
    result_list = [
            {
                "protein_id": protein_id,
                "peptide_seq": protein_seq,
                "lysine_position": lysine_position,
                "nonsumoylation_class_probs": nonsumoylation_class_probs,
                "sumoylation_class_probs": sumoylation_class_probs,
                "predicted_labels": predicted_labels
            }
            for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels'])
        ]
    return Response({"data": result_list}, status=status.HTTP_200_OK)


@api_view(['POST'])
@parser_classes([MultiPartParser])

def denemeFastaFile(request):
    file_obj = request.FILES.get('file')
    if file_obj:
        file_name, file_extension = file_obj.name.split('.')

        if str(file_extension.lower()) != 'fasta':
            return Response({'error': 'Invalid file extension. Only fasta files are allowed.'}, status=status.HTTP_400_BAD_REQUEST)

        records = file_obj.read().decode('utf-8')
        fasta_file = StringIO(records)
        data_processes = Data()
        jsonList = []
        idS, seqS, positionS = [], [], []
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_ids, protein_seqs, k_positions = seqIOParser(record)
            idS.extend(protein_ids)
            seqS.extend(protein_seqs)
            positionS.extend(k_positions)
        
        df = make_prediction(idS, seqS, positionS)  
        result_list = [
                {
                    "protein_id": protein_id,
                    "peptide_seq": protein_seq,
                    "lysine_position": lysine_position,
                    "nonsumoylation_class_probs": nonsumoylation_class_probs,
                    "sumoylation_class_probs": sumoylation_class_probs,
                    "predicted_labels": predicted_labels
                }
                for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels'])
            ]
        return Response({"data": result_list}, status=status.HTTP_200_OK)
    else:
        return Response({"message": "An error occured."}, status=status.HTTP_400_BAD_REQUEST)