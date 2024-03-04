from sumonet.utils.data_pipe import Data
from rest_framework.decorators import api_view, parser_classes
from sumonetWeb.serializers import UniprotSerializer
from rest_framework import status
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser
from io import StringIO
from Bio import SeqIO
from .helpers import make_prediction, seqIOParser

@api_view(['POST'])
def uniprotPrediction(request):
    
    serializer = UniprotSerializer(data=request.data)
    
    if serializer.is_valid():
        data_processes = Data()
        uniprot_id = serializer.data['uniprot_id']  

        try:
            lysine_position = serializer.data['lysine_position']

            lysine_position = int(lysine_position) # This line will raise a ValueError if the conversion fails

        except KeyError:
            lysine_position = None

        except ValueError: #++
            # Return a 400 Bad Request response if lysine_position is not an integer
            return Response({'error': 'Lysine position must be an integer.'}, status=status.HTTP_400_BAD_REQUEST)
        
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

        try: #++
            df = make_prediction(protein_ids, protein_seqs, k_positions)
        except TypeError:
            return Response({'error': 'Invalid Lysine Position'}, status=status.HTTP_400_BAD_REQUEST)

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
            
        #return Response(result, status=status.HTTP_200_OK)
        #return Response({"data": result, "length": len(result)}, status=status.HTTP_200_OK)
        return Response({"data": result}, status=status.HTTP_200_OK) #++
            
    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
           
        
        
        
@api_view(['POST'])
def proteinSequence(request):
    protein_seq = request.data.get('protein_seq', '')


    if protein_seq == '' or protein_seq == None:
        return Response({'error': 'Protein sequence must be entered.'}, status=status.HTTP_400_BAD_REQUEST)

    fasta_file = StringIO(protein_seq)
    idS, seqS, positionS = [], [], []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_ids, protein_seqs, k_positions = seqIOParser(record)
        idS.extend(protein_ids)
        seqS.extend(protein_seqs)
        positionS.extend(k_positions)
        
        
    try:
        df = make_prediction(idS, seqS, positionS) 

    except ValueError:
        return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
     
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
    return Response({"data": result_list, "length": len(result_list)}, status=status.HTTP_200_OK)


@api_view(['POST'])
@parser_classes([MultiPartParser])
def fastaFile(request):
    file_obj = request.FILES.get('file')
    if file_obj:
        file_name, file_extension = file_obj.name.split('.')

        if str(file_extension.lower()) != 'fasta':
            return Response({'error': 'Invalid file extension. Only fasta files are allowed.'}, status=status.HTTP_400_BAD_REQUEST)

        records = file_obj.read().decode('utf-8')

        if records == '' or records == None:
            return Response({'error': 'This file does not contain protein sequence(s).'}, status=status.HTTP_400_BAD_REQUEST)
        


        fasta_file = StringIO(records)
        idS, seqS, positionS = [], [], []
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_ids, protein_seqs, k_positions = seqIOParser(record)
            idS.extend(protein_ids)
            seqS.extend(protein_seqs)
            positionS.extend(k_positions)
        
        try:

            df = make_prediction(idS, seqS, positionS)  

        except ValueError:
            return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
        
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
        return Response({"data": result_list, "length": len(result_list)}, status=status.HTTP_200_OK)
    else:
        return Response({"error": "An error occured."}, status=status.HTTP_400_BAD_REQUEST)