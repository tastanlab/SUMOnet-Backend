from django.conf import settings
from django.test import TestCase
from rest_framework.test import APIClient
from rest_framework import status
from django.urls import reverse

class TestSumonet(TestCase):

    def setUp(self):
        self.client = APIClient()
        
    #! Uniprot Prediction Test Cases

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

    def test_EmptyLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': 'O00566', 'lysine_position': ''}

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
        
        
    def test_EmptyUniprotIdAndLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})
        
        
    def test_NoUniprotIdAndLysinePosition(self):
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': '', 'lysine_position': ''}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})
        
        
    #! Uniprot Prediction Test Cases End 
    
    
    #! Protein Prediction Test Cases
    
    def test_EmptyProteinSequence(self):
        url = reverse('protein_sequence_prediction')
        data = {'protein_id': ''}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Protein sequence must be entered.'})
        
        
    def test_NoProteinSequence(self):
        url = reverse('protein_sequence_prediction')
        data = {}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Protein sequence must be entered.'})
        

    def test_InvalidProteinSequence(self):
        url = reverse('protein_sequence_prediction')
        data = {'protein_seq': 'AAAAA'}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Invalid protein sequence.'})
        
        
    def test_ValidProteinSequence(self):
        url = reverse('protein_sequence_prediction')
        data = {'protein_seq': '>sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHG\r\nSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGS\r\nEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKL\r\nEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDI\r\nDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEE\r\nAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSF\r\nEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEE\r\nTTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKL\r\nNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAP\r\nVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLE\r\nKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQIND\r\nAKKTEKKKKKRQDISVHKLKL\r\n\r\n>sp|Q9UER7|DAXX_HUMAN Death domain-associated protein 6 OS=Homo sapiens OX=9606 GN=DAXX PE=1 SV=2\r\nMATANSIIVLDDDDEDEAAAQPGPSHPLPNAASPGAEAPSSSEPHGARGSSSSGGKKCYK\r\nLENEKLFEEFLELCKMQTADHPEVVPFLYNRQQRAHSLFLASAEFCNILSRVLSRARSRP\r\nAKLYVYINELCTVLKAHSAKKKLNLAPAATTSNEPSGNNPPTHLSLDPTNAENTASQSPR\r\nTRGSRRQIQRLEQLLALYVAEIRRLQEKELDLSELDDPDSAYLQEARLKRKLIRLFGRLC\r\nELKDCSSLTGRVIEQRIPYRGTRYPEVNRRIERLINKPGPDTFPDYGDVLRAVEKAAARH\r\nSLGLPRQQLQLMAQDAFRDVGIRLQERRHLDLIYNFGCHLTDDYRPGVDPALSDPVLARR\r\nLRENRSLAMSRLDEVISKYAMLQDKSEEGERKKRRARLQGTSSHSADTPEASLDSGEGPS\r\nGMASQGCPSASRAETDDEDDEESDEEEEEEEEEEEEEATDSEEEEDLEQMQEGQEDDEEE\r\nDEEEEAAAGKDGDKSPMSSLQISNEKNLEPGKQISRSSGEQQNKGRIVSPSLLSEEPLAP\r\nSSIDAESNGEQPEELTLEEESPVSQLFELEIEALPLDTPSSVETDISSSRKQSEEPFTTV\r\nLENGAGMVSSTSFNGGVSPHNWGDSGPPCKKSRKEKKQTGSGPLGNSYVERQRSVHEKNG\r\nKKICTLPSPPSPLASLAPVADSSTRVDSPSHGLVTSSLCIPSPARLSQTPHSQPPRPGTC\r\nKTSVATQCDPEEIIVLSDSD\r\n\r\n\r\n'}
        
        response = self.client.post(url, data, format="json")
        
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertIsInstance(response.data, dict)
        

    #! Protein Prediction Test Cases End
    