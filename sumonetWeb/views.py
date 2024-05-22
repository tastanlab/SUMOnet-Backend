
from sumonetML.sumonet.utils.data_pipe import Data
from rest_framework.decorators import api_view, parser_classes
from rest_framework import status
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser
from io import StringIO
from Bio import SeqIO
from .helpers import make_prediction, seqIOParser
import re


#! @desc: This function takes a uniprot id and lysine position as input 
#! and returns the sumoylation prediction for the lysine position in the protein sequence.
#! However, if user does not give the lysine position, it will return the sumoylation 
#! prediction for all lysine positions in the protein sequence.

#! route: POST /uniprot-prediction/

#! @access: Public

def is_valid_uniprot_id(uniprot_id):
    return re.fullmatch(r'[a-zA-Z0-9 ,]*', uniprot_id) is not None

def An_Uniprot_Id_Predictor(lysine_position, protein_seq, uniprot_id, threshold):
    data_processes = Data()

    try:
        if lysine_position != '' and lysine_position is not None:
            protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq, uniprot_id, lysine_position)
        else:
            protein_ids, protein_seqs, k_positions = data_processes.uniprot_id_input(protein_seq, uniprot_id)

        df = make_prediction(protein_ids, protein_seqs, k_positions)

    except IndexError:
        return {'error': 'Invalid Lysine Position. Index out of range!'}
    except ValueError:
        return {'error': 'Invalid Lysine Position. Lysine Position must be positive!'}
    except TypeError:
        return {'error': 'No data found with this lysine position and Uniprot ID!'}

    result = []
    for protein_id, protein_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels']):
        if sumoylation_class_probs >= threshold:
            result.append({
                "protein_id": uniprot_id,
                "peptide_seq": protein_seq,
                "lysine_position": lysine_position,
                "nonsumoylation_class_probs": nonsumoylation_class_probs,
                "sumoylation_class_probs": sumoylation_class_probs,
                "predicted_labels": predicted_labels
            })

    return result



@api_view(['POST'])
def uniprotPrediction(request):
   
    at_least_one_valid_uniprot_id = False
    uniprot_id = request.data.get('uniprot_id', '')
    lysine_position = request.data.get('lysine_position', '')
    threshold = request.data.get('threshold', 0.5)

    try:
        threshold = float(threshold)
    except ValueError:
        return Response({'error': 'Threshold must be a number.'}, status=status.HTTP_400_BAD_REQUEST)
    
    if not (0 <= threshold <= 1):
        return Response({'error': 'Threshold must be between 0 and 1.'}, status=status.HTTP_400_BAD_REQUEST)

    data_processes = Data()
    
    if uniprot_id == '' or uniprot_id == None:
        return Response({'error': 'Uniprot ID must be entered!'}, status=status.HTTP_400_BAD_REQUEST)
    
    elif not is_valid_uniprot_id(uniprot_id):
        return Response({'error': 'You can only use comma for entering Uniprot ID list!'}, status=status.HTTP_400_BAD_REQUEST)

    elif "," in uniprot_id:
        if is_valid_uniprot_id(uniprot_id): # Logic is working in the previous elif line. However, this line was put as a precaution. It is for character checker not uniprot id.
            if lysine_position == "" or lysine_position == None :
                result_list = []
                invalid_uniprot_id_list = []
                uniquness_list = []
                uniprot_id_list = uniprot_id.replace(' ', '').split(",")
                print(uniprot_id_list)
                for uid in uniprot_id_list:
                    if uid not in uniquness_list:
                        uniquness_list.append(uid)
                        protein_seq = data_processes.retrive_protein_sequence_with_uniprotid(uid)
                        if (protein_seq == '' or protein_seq == None):
                            invalid_uniprot_id_list.append(uid)
                            pass
                        else:
                            at_least_one_valid_uniprot_id = True
                            result = An_Uniprot_Id_Predictor(lysine_position, protein_seq, uid, threshold)
                            result_list.append(result)
                flattened_list = [item for sublist in result_list for item in sublist]
                if at_least_one_valid_uniprot_id:
                    return Response({"data": flattened_list, "length": len(flattened_list), "invalid_idS": invalid_uniprot_id_list}, status=status.HTTP_200_OK)
                elif not at_least_one_valid_uniprot_id:
                    return Response({'error': 'No protein sequence found with this Uniprot ID!'}, status=status.HTTP_400_BAD_REQUEST)

            elif not lysine_position == "":
                return Response({'error': 'You cannot enter lysine position for multiple Uniprot ID!'}, status=status.HTTP_400_BAD_REQUEST)
        else:
            return Response({'error': 'You can only use comma for entering Uniprot ID list!'}, status=status.HTTP_400_BAD_REQUEST)
            
    else:
        protein_seq = data_processes.retrive_protein_sequence_with_uniprotid(uniprot_id)

        if protein_seq == '' or protein_seq == None:
            return Response({'error': 'No protein sequence found with this Uniprot ID!'}, status=status.HTTP_400_BAD_REQUEST)

        

        if lysine_position != '' and lysine_position != None:
            try:
                lysine_position = int(lysine_position)
            except KeyError:
                return Response({'error': 'Invalid Lysine Position!'}, status=status.HTTP_400_BAD_REQUEST)
            except ValueError:
                return Response({'error': 'Lysine position must be an integer!'}, status=status.HTTP_400_BAD_REQUEST)
            except IndexError:
                return Response({'error': 'Invalid Lysine Position!'}, status=status.HTTP_400_BAD_REQUEST)

            if lysine_position < 0: # This line checks the negativity of the integer other exceptions hold for only precaution. 
                return Response({'error': 'Invalid Lysine Position. Lysine Position must be positive!'}, status=status.HTTP_400_BAD_REQUEST)
            
            if lysine_position > len(protein_seq):
                return Response({'error': 'Invalid Lysine Position. Index out of range!'}, status=status.HTTP_400_BAD_REQUEST)

        result = An_Uniprot_Id_Predictor(lysine_position, protein_seq, uniprot_id, threshold)
        if 'error' in result:
            return Response(result, status=status.HTTP_400_BAD_REQUEST)

        return Response({"data": result, "length": len(result)}, status=status.HTTP_200_OK)
    
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
    invalid_responses = []  # Store invalid responses

    if not protein_seq:
        invalid_responses.append({'error': 'Protein sequence must be entered!'})
        return Response({"data": [], "errors": invalid_responses}, status=status.HTTP_400_BAD_REQUEST)

    fasta_file = StringIO(protein_seq)
    idS, seqS, positionS = [], [], []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) < 15:
            invalid_responses.append({'error': 'Protein sequence must be at least 15 amino acids long!', 'ids': [record.id]})
            continue

        protein_ids, protein_seqs, k_positions = seqIOParser(record)
        idS.extend(protein_ids)
        seqS.extend(protein_seqs)
        positionS.extend(k_positions)


        try:
            df = make_prediction(idS, seqS, positionS) 
            
        except (ValueError, KeyError) :
            invalid_responses.append({'error': "Invalid Protein Sequence.", 'ids': record.id})

    result_list = []
    if idS:  # Proceed only if there are valid sequences to predict
        
        result_list = [
            {
                "protein_id": protein_id.split('|')[1] if '|' in protein_id else protein_id,
                "peptide_seq": peptide_seq,
                "lysine_position": lysine_position,
                "nonsumoylation_class_probs": nonsumoylation_class_probs,
                "sumoylation_class_probs": sumoylation_class_probs,
                "predicted_labels": predicted_labels
            }
            for protein_id, peptide_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(
                df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels']
            )
        ]
       

    if  len(result_list) == 0 and  len(invalid_responses) == 0:  # If no valid data and no prior errors, assume all input was invalid
        return Response({"error": "No valid protein sequences provided."}, status=status.HTTP_400_BAD_REQUEST)

    elif len(result_list) == 0 and len(invalid_responses) > 0:
        return Response({"error": "No valid protein sequences provided."}, status=status.HTTP_400_BAD_REQUEST)

    # Always return a response including both results and any errors
    # response_status = status.HTTP_200_OK if result_list else status.HTTP_400_BAD_REQUEST
    return Response({"data": result_list, "length": len(result_list), "errors": invalid_responses}, status=status.HTTP_200_OK)



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
    invalid_responses = []  # Store invalid responses
    if file_obj:
        file_name, file_extension = file_obj.name.split('.')

        if file_extension.lower() != "fasta" and file_extension.lower() != "txt":
            return Response({'error': 'File must be in .fasta or .txt format!'}, status=status.HTTP_400_BAD_REQUEST)

       
        if str(file_extension.lower()) == "fasta" or str(file_extension.lower()) == "txt" :

            records = file_obj.read().decode('utf-8')

            if records == '' or records == None:
                return Response({'error': 'This file does not contain protein sequence(s)!'}, status=status.HTTP_400_BAD_REQUEST)
            #print(records)
            


            fasta_file = StringIO(records)
            idS, seqS, positionS = [], [], []
            for record in SeqIO.parse(fasta_file, "fasta"):

                #print("record", record.seq)

                if (len(record.seq) < 15):
                    invalid_responses.append({'error': 'Protein sequence must be at least 15 amino acids long!', 'ids': [record.id]})
                    continue

                

                protein_ids, protein_seqs, k_positions = seqIOParser(record)
                idS.extend(protein_ids)
                seqS.extend(protein_seqs)
                positionS.extend(k_positions)
            
                try:

                    df = make_prediction(idS, seqS, positionS)  

                except (ValueError, KeyError) :
                    invalid_responses.append({'error': "Invalid Protein Sequence.", 'ids': record.id})

           
            
            result_list = []
            if idS:  # Proceed only if there are valid sequences to predict
                    
                result_list = [
                    {
                        "protein_id": protein_id.split('|')[1] if '|' in protein_id else protein_id,
                        "peptide_seq": peptide_seq,
                        "lysine_position": lysine_position,
                        "nonsumoylation_class_probs": nonsumoylation_class_probs,
                        "sumoylation_class_probs": sumoylation_class_probs,
                        "predicted_labels": predicted_labels
                    }
                    for protein_id, peptide_seq, lysine_position, nonsumoylation_class_probs, sumoylation_class_probs, predicted_labels in zip(
                        df['protein_id'], df['protein_seq'], df['lysine_position'], df['nonsumoylation_class_probs'], df['sumoylation_class_probs'], df['predicted_labels']
                    )
                ]
        if  len(result_list) == 0 and  len(invalid_responses) == 0:  # If no valid data and no prior errors, assume all input was invalid
            return Response({"error": "No valid protein sequences provided."}, status=status.HTTP_400_BAD_REQUEST)

        elif len(result_list) == 0 and len(invalid_responses) > 0:
            return Response({"error": "No valid protein sequences provided."}, status=status.HTTP_400_BAD_REQUEST)
        
        return Response({"data": result_list, "length": len(result_list), "errors": invalid_responses}, status=status.HTTP_200_OK)
    else:
        return Response({'error': 'File must be in .fasta or .txt format!'}, status=status.HTTP_400_BAD_REQUEST)
