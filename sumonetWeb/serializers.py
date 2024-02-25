from rest_framework import serializers

class UniprotSerializer(serializers.Serializer):
    uniprot_id = serializers.CharField(max_length=100)
    lysine_position = serializers.CharField(required=False ,max_length=100)
    
    
class ProteinSequenceSerializer(serializers.Serializer):
    protein_seq = serializers.ListField(child=serializers.CharField(max_length=1000))
    