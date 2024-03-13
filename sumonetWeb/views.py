
from sumonetML.sumonet.utils.data_pipe import Data
from rest_framework.decorators import api_view, parser_classes
from rest_framework import status
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser
from io import StringIO
from Bio import SeqIO
from .helpers import make_prediction, seqIOParser


#! @desc: This function takes a uniprot id and lysine position as input 
#! and returns the sumoylation prediction for the lysine position in the protein sequence.
#! However, if user does not give the lysine position, it will return the sumoylation 
#! prediction for all lysine positions in the protein sequence.

#! route: POST /uniprot-prediction/

#! @access: Public

@api_view(['POST'])
def uniprotPrediction(request):
   
    uniprot_id = request.data.get('uniprot_id', '')
    lysine_position = request.data.get('lysine_position', '')
    data_processes = Data()
    
    if uniprot_id == '' or uniprot_id == None:
        return Response({'error': 'Uniprot ID must be entered.'}, status=status.HTTP_400_BAD_REQUEST)
    
    else:
        protein_seq = data_processes.retrive_protein_sequence_with_uniprotid(uniprot_id)

        if protein_seq == '' or protein_seq == None:
            return Response({'error': 'No protein sequence found with this uniprot id.'}, status=status.HTTP_400_BAD_REQUEST)

        

        if lysine_position != '' and lysine_position != None:
            try:
                lysine_position = int(lysine_position)
            except KeyError:
                return Response({'error': 'Invalid Lysine Position.'}, status=status.HTTP_400_BAD_REQUEST)
            except ValueError:
                return Response({'error': 'Lysine position must be an integer.'}, status=status.HTTP_400_BAD_REQUEST)
            except IndexError:
                return Response({'error': 'Invalid Lysine Position.'}, status=status.HTTP_400_BAD_REQUEST)
            
            
        try:
            if lysine_position != '' and lysine_position != None:
                protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq, uniprot_id, lysine_position)
                
            else:
                protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq, uniprot_id)

            df = make_prediction(protein_ids, protein_seqs, k_positions)

        except IndexError:
            return Response({'error': 'Invalid Lysine Position. Index out of range'}, status=status.HTTP_400_BAD_REQUEST)
        #except ValueError:
            #return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
        except TypeError:
            return Response({'error': 'No data found with this lysine position and uniprot id.'}, status=status.HTTP_400_BAD_REQUEST)

        
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
        return Response({"data": result, "length": len(result)}, status=status.HTTP_200_OK) #++
    
    return Response({'error': 'An error occured.'}, status=status.HTTP_400_BAD_REQUEST)
        
        

#! @desc: This function takes a protein sequence as input and returns the 
#! sumoylation prediction for all lysine positions in the protein sequence.
#! You must be careful with the input format. The input must be in fasta format. 
#! The protein sequence must be in the second line of the fasta file.

#! route: POST /protein-sequence-prediction/

#! @access: Public

@api_view(['POST'])
def proteinSequence(request):
    protein_seq = request.data.get('protein_seq', '')


    if protein_seq == '' or protein_seq == None:
        return Response({'error': 'Protein sequence must be entered.'}, status=status.HTTP_400_BAD_REQUEST)

    fasta_file = StringIO(protein_seq)

    idS, seqS, positionS = [], [], []
    for record in SeqIO.parse(fasta_file, "fasta"):

        if (len(record.seq) < 15):
            return Response({'error': 'Protein sequence must be at least 15 amino acids long.'}, status=status.HTTP_400_BAD_REQUEST)
     
        
        protein_ids, protein_seqs, k_positions = seqIOParser(record)
        idS.extend(protein_ids)
        seqS.extend(protein_seqs)
        positionS.extend(k_positions)
        
    
    try:
        df = make_prediction(idS, seqS, positionS) 

    except ValueError:
        return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
    except KeyError:
         return Response({'error': 'Invalid protein sequence.'}, status=status.HTTP_400_BAD_REQUEST)
     
    result_list = [
            {
                "protein_id": protein_id.split('|')[1] if '|' in protein_id else protein_id,
                "peptide_seq": protein_seq,
                "lysine_position": lysine_position,
                "nonsumoylation_class_probs": nonsumoylation_class_probs,
                "sumoylation_class_probs": sumoylation_class_probs,
                "predicted_labels": predicted_labels
            }
            for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels'])
        ]
    return Response({"data": result_list, "length": len(result_list)}, status=status.HTTP_200_OK)



#! @desc: This function takes a fasta file as input and returns the 
#! sumoylation prediction for all lysine positions in the protein sequence.
#! You must be careful with the input format. The protein sequence must be in the second line of the fasta file.
#! Please leave at least one line of space between the protein sequence and the next protein sequence(s).

#! route: POST /fasta-file-prediction/

#! @access: Public

@api_view(['POST'])
@parser_classes([MultiPartParser])
def fastaFile(request):
    file_obj = request.FILES.get('file')
    if file_obj:
        file_name, file_extension = file_obj.name.split('.')

       
        if str(file_extension.lower()) != "fasta":
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
        return Response({"error": "Please attach a file."}, status=status.HTTP_400_BAD_REQUEST)
