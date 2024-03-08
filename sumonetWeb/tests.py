from django.conf import settings
from django.test import TestCase
from rest_framework.test import APIClient
from rest_framework import status
from django.urls import reverse

class TestSumonet(TestCase):

    def setUp(self):
        self.client = APIClient()

    def test_EmptyUniprotID(self):
        url = reverse('uniprot_prediction')  # Assuming you have a name for your URL pattern
        data = {'uniprot_id': ''}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})


    def test_NoUniprotID(self):
        url = reverse('uniprot_prediction')  # Assuming you have a name for your URL pattern
        data = {'lysine_position': '20'}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})


    def test_NoLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': 'O00566'}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertIsInstance(response.data, dict)


    def test_InvalidStringLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': 'O00566', 'lysine_position': 'aaaa'}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Lysine position must be an integer.'})

    def test_InvalidIntegersLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': 'O00566', 'lysine_position': 99998444}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Invalid Lysine Position. Index out of range'})