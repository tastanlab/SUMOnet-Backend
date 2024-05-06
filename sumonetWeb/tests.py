from django.conf import settings
from django.test import TestCase
from rest_framework.test import APIClient
from rest_framework import status
from django.urls import reverse
from django.core.files.uploadedfile import SimpleUploadedFile
# pip install ddt

class TestSumonet(TestCase):

    def setUp(self):
        self.client = APIClient()
        
    #! Uniprot Prediction Test Cases

    def test_EmptyUniprotID(self): # Test for empty uniprot ID returns 400
        url = reverse('uniprot_prediction') 
        data = {'uniprot_id': ''}
        
        response = self.client.post(url, data, format='json')
        
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})

    def test_EmptyUniprotIDwithLysinePosition(self): # test for empty uniprot ID with Lysine Postition returns 400
        url = reverse('uniprot_prediction')  
        data = [{'lysine_position': '20'}, {'lysine_position': '25'}, 
        {'lysine_position': '30'}, {'lysine_position': '35'}, {'lysine_position': '40'}]
        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})
    
    def test_EmptyUniprotIdAndLysinePosition(self): # test for empty JSON returns 400
        url = reverse('uniprot_prediction')
        data = {}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})
    
    def test_NoUniprotIdAndLysinePosition(self): # Test for both empty returns 400
        url = reverse('uniprot_prediction')
        data = {'uniprot_id': '', 'lysine_position': ''}

        response = self.client.post(url, data, format='json')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Uniprot ID must be entered.'})
  
    def test_NoLysinePosition(self): # test for valid uniprot ID's without lysine position
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566'}, {'uniprot_id': 'A0JNW5'}, {'uniprot_id': 'E0CX11'},
        {'uniprot_id': 'A0A1B0GTW7'}, {'uniprot_id': 'O94762'} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)
    
    def test_EmptyLysinePosition(self): # test for valid uniprot ID and empty lysine position
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566', 'lysine_position': ''}, {'uniprot_id': 'A0JNW5' , 'lysine_position': ''}, {'uniprot_id': 'E0CX11' , 'lysine_position': ''},
        {'uniprot_id': 'A0A1B0GTW7', 'lysine_position': ''}, {'uniprot_id': 'O94762' , 'lysine_position': ''} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)

    def test_ValidUniprodIdLysinePositionPairs(self): # test for valid uniprot ID and empty lysine position
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566', 'lysine_position': '20'}, {'uniprot_id': 'A0JNW5' , 'lysine_position': '24'}, {'uniprot_id': 'E0CX11' , 'lysine_position': '36'},
        {'uniprot_id': 'A0A1B0GTW7', 'lysine_position': '197'}, {'uniprot_id': 'O94762' , 'lysine_position': '110'} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)
    

    def test_ValidUniprotList(self):
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': "O00566, A0JP26"}, {'uniprot_id': 'O00566, A0JP26 ,  O14763'}, {'uniprot_id': 'O00566 , O14763' },
        {'uniprot_id': 'O00566 , O14763, Q6GZV6'}, {'uniprot_id': 'O94762'} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)
  
    def test_ValidUniprotList(self):
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': "O00566, O00566"}, {'uniprot_id': '000566, O00566'}, {'uniprot_id': 'O00566, 000566' },
        {'uniprot_id': 'O00566, AoJP26, O14763'}, {'uniprot_id': 'O00566, AoJP26,  014763'}, {'uniprot_id': 'AoJP26, O00566, 014763'}, {'uniprot_id': 'AoJP26, 014763 ,O00566'} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)
                
                response_json = response.json()
                
                self.assertIn('invalid_idS', response_json)
                self.assertIsInstance(response_json['invalid_idS'], list)
    
    def test_ValidUniprotListWithLysinePosition(self): # negative lysine position test
        url = reverse('uniprot_prediction')
        data = [
        { "uniprot_id": "O00566, A0JP26", "lysine_position": 20 },
        { "uniprot_id": "O00566, A0JP26", "lysine_position": "1,2,3" },
        { "uniprot_id": "O00566, A0JP26", "lysine_position": "20,17" },
        { "uniprot_id": "O00566, A0JP26", "lysine_position": "aaaaa,bbbb" },
        { "uniprot_id": "O00566, A0JP26", "lysine_position": "aaa;bbb" },
        ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': "You cannot enter lysine position for multiple Uniprot ID."}) 

 
    def test_InvalidStringLysinePosition(self): # Invalid lysine position for string
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566', 'lysine_position': 'aaaa'}, {'uniprot_id': 'O00566', 'lysine_position': 'abcdef'}, {'uniprot_id': 'O00566', 'lysine_position': 'ghty'} ]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Lysine position must be an integer.'})

       
    def test_InvalidIntegersLysinePosition(self): # out of index test for lysine position
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566', 'lysine_position': 1000000000}, {'uniprot_id': 'O00566', 'lysine_position': 9999999}]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Invalid Lysine Position. Index out of range'})

    def test_NegativeIntegersLysinePosition(self): # negative lysine position test
        url = reverse('uniprot_prediction')
        data = [{'uniprot_id': 'O00566', 'lysine_position': -1} , {'uniprot_id': 'O00566', 'lysine_position': -3}, {'uniprot_id': 'O00566', 'lysine_position': -2}, {'uniprot_id': 'O00566', 'lysine_position': -4}, {'uniprot_id': 'O00566', 'lysine_position': -100000000}]

        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('uniprot_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Invalid Lysine Position. Lysine Position must be positive.'})  
         

      
      
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
        data = [{'protein_seq': 'AAAAA'}, {'protein_seq': 'BBBBB'}, {'protein_seq': '|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPP10 PE=1 SV=2\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHGSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGSEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKLEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDIDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEEAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSFEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEETTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKLNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAPVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLEKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQINDAKKTEKKKKKRQDISVHKLKL'}]
        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('protein_sequence_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Invalid protein sequence.'}) 
         
    def test_ValidProteinSequence(self):
        url = reverse('protein_sequence_prediction')
        data = [{'protein_seq': '>sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHG\r\nSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGS\r\nEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKL\r\nEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDI\r\nDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEE\r\nAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSF\r\nEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEE\r\nTTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKL\r\nNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAP\r\nVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLE\r\nKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQIND\r\nAKKTEKKKKKRQDISVHKLKL\r\n\r\n>sp|Q9UER7|DAXX_HUMAN Death domain-associated protein 6 OS=Homo sapiens OX=9606 GN=DAXX PE=1 SV=2\r\nMATANSIIVLDDDDEDEAAAQPGPSHPLPNAASPGAEAPSSSEPHGARGSSSSGGKKCYK\r\nLENEKLFEEFLELCKMQTADHPEVVPFLYNRQQRAHSLFLASAEFCNILSRVLSRARSRP\r\nAKLYVYINELCTVLKAHSAKKKLNLAPAATTSNEPSGNNPPTHLSLDPTNAENTASQSPR\r\nTRGSRRQIQRLEQLLALYVAEIRRLQEKELDLSELDDPDSAYLQEARLKRKLIRLFGRLC\r\nELKDCSSLTGRVIEQRIPYRGTRYPEVNRRIERLINKPGPDTFPDYGDVLRAVEKAAARH\r\nSLGLPRQQLQLMAQDAFRDVGIRLQERRHLDLIYNFGCHLTDDYRPGVDPALSDPVLARR\r\nLRENRSLAMSRLDEVISKYAMLQDKSEEGERKKRRARLQGTSSHSADTPEASLDSGEGPS\r\nGMASQGCPSASRAETDDEDDEESDEEEEEEEEEEEEEATDSEEEEDLEQMQEGQEDDEEE\r\nDEEEEAAAGKDGDKSPMSSLQISNEKNLEPGKQISRSSGEQQNKGRIVSPSLLSEEPLAP\r\nSSIDAESNGEQPEELTLEEESPVSQLFELEIEALPLDTPSSVETDISSSRKQSEEPFTTV\r\nLENGAGMVSSTSFNGGVSPHNWGDSGPPCKKSRKEKKQTGSGPLGNSYVERQRSVHEKNG\r\nKKICTLPSPPSPLASLAPVADSSTRVDSPSHGLVTSSLCIPSPARLSQTPHSQPPRPGTC\r\nKTSVATQCDPEEIIVLSDSD\r\n\r\n\r\n'}, {'protein_seq': ">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPP10 PE=1 SV=2\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHGSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGSEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKLEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDIDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEEAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSFEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEETTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKLNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAPVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLEKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQINDAKKTEKKKKKRQDISVHKLKL"}]
        
        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('protein_sequence_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_200_OK)
                self.assertIsInstance(response.data, dict)
    def test_ProteinSequenceSmallerThan15(self):
        url = reverse('protein_sequence_prediction')
        data = [{'protein_seq': '>sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nNKSIVDDKFEEE'}, {"protein_seq" : ">sp|Q9UER7|DAXX_HUMAN Death domain-associated protein 6 OS=Homo sapiens OX=9606 GN=DAXX PE=1 SV=2\r\nMATANSIIVL"}]
        
        for test_case in data:
            with self.subTest(test_case=test_case):
                url = reverse('protein_sequence_prediction')
                response = self.client.post(url, test_case, format='json')

                self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
                self.assertEqual(response.data, {'error': 'Protein sequence must be at least 15 amino acids long.'}) 

    #! Protein Prediction Test Cases End
    

    #! File Upload Test Cases
    def test_FileUploadWithProteinSequenceSmaller_15(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPP10 PE=1 SV=2\nMNAQPWWRRTT"
        file = SimpleUploadedFile('test.txt', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Protein sequence must be at least 15 amino acids long.'}) 

    def test_FileUploadForTxt(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHG\r\nSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGS\r\nEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKL\r\nEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDI\r\nDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEE\r\nAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSF\r\nEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEE\r\nTTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKL\r\nNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAP\r\nVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLE\r\nKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQIND\r\nAKKTEKKKKKRQDISVHKLKL"
        file = SimpleUploadedFile('test.txt', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertIsInstance(response.data, dict)

    def test_FileUploadForFasta(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHG\r\nSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGS\r\nEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKL\r\nEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDI\r\nDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEE\r\nAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSF\r\nEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEE\r\nTTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKL\r\nNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAP\r\nVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLE\r\nKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQIND\r\nAKKTEKKKKKRQDISVHKLKL"
        file = SimpleUploadedFile('test.fasta', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertIsInstance(response.data, dict)
    
    def test_FileUploadForPdf(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2\r\nMAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKILENGRIHG\r\nSPLQKLVIENFDDEQIWQQLELQNEPILQYFQNAVSETINDEDISLLPESEEQEREEDGS\r\nEIEADDKEDLEDLEEEEVSDMGNDDPEMGERAENSSKSDLRKSPVFSDEDSDLDFDISKL\r\nEQQSKVQNKGQGKPREKSIVDDKFFKLSEMEAYLENIEKEEERKDDNDEEEEDIDFFEDI\r\nDSDEDEGGLFGSKKLKSGKSSRNLKYKDFFDPVESDEDITNVHDDELDSNKEDDEIAEEE\r\nAEELSISETDEDDDLQENEDNKQHKESLKRVTFALPDDAETEDTGVLNVKKNSDEVKSSF\r\nEKRQEKMNEKIASLEKELLEKKPWQLQGEVTAQKRPENSLLEETLHFDHAVRMAPVITEE\r\nTTLQLEDIIKQRIRDQAWDDVVRKEKPKEDAYEYKKRLTLDHEKSKLSLAEIYEQEYIKL\r\nNQQKTAEEENPEHVEIQKMMDSLFLKLDALSNFHFIPKPPVPEIKVVSNLPAITMEEVAP\r\nVSVSDAALLAPEEIKEKNKAGDIKTAAEKTATDKKRERRKKKYQKRMKIKEKEKRRKLLE\r\nKSSVDQAGKYSKTVASEKLKQLTKTGKASFIKDEGKDKALKSSQAFFSKLQDQVKMQIND\r\nAKKTEKKKKKRQDISVHKLKL"
        file = SimpleUploadedFile('test.pdf', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Invalid file extension. Only fasta and txt files are allowed.'}) 

    def test_FileUploadForTxtInvalid(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b"xxx"
        file = SimpleUploadedFile('test.txt', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Invalid protein sequence.'}) 

    def test_FileUploadForFastaInvalid(self):
        url = reverse('fasta_file_prediction')  

        file_includes = b"xxx"
        file = SimpleUploadedFile('test.fasta', file_includes, content_type='text/plain')

        response = self.client.post(url, {'file': file}, format='multipart')

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(response.data, {'error': 'Invalid protein sequence.'}) 

#! File Upload Test Cases End