from django.db import models


class SUMO(models.Model):
    protein_id = models.CharField(max_length=100)
    peptide_sequence = models.CharField(max_length=100)
    lysine_position = models.IntegerField()
    nonsumoylation_class_probs = models.DecimalField(max_digits=10, decimal_places=9)
    sumoylation_class_probs = models.DecimalField(max_digits=10, decimal_places=9)
    predicted_labels = models.IntegerField()
    