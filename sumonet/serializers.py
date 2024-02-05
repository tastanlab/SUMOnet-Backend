from rest_framework import serializers
from .models import SUMO

class SUMOSerializer(serializers.ModelSerializer):
    class Meta:
        model = SUMO
        fields = ['id', 'protein_id', 'peptide_sequence', 'lysine_position', 'nonsumoylation_class_probs', 'sumoylation_class_probs', 'predicted_labels']
        