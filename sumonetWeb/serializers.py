# Probably this classes will not be needed since we have nothing to serialize
from rest_framework import serializers

class UniprotSerializer(serializers.Serializer):
    uniprot_id = serializers.CharField(max_length=100)
    lysine_position = serializers.CharField(required=False ,max_length=100)

    default_error_messages = {
        'invalid': 'Invalid Lysine Position',
    }
    
    
class ProteinSequenceSerializer(serializers.Serializer):
    protein_seq = serializers.ListField(child=serializers.CharField(max_length=1000))
    