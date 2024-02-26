# Tutorial

## Setting up the environment

To use the necessary dependencies, we need to set up the virtual environment.

### For Windows users

Write the following line to the terminal: `python -m venv env`

### For Unix/macOS users

Write the following line to the terminal: `python3 -m venv env`

## Activating the virtual environment

To install our dependencies, we need to activate the virtual environment. (You must be in the same directory with the env folder. You can check that up from terminal.)

### For Windows users

Write the following line to the terminal: `env\Scripts\activate`

### For Unix/macOS users

Write the following line to the terminal: `source env/bin/activate`

## Installing the dependencies

After you set up the virtual environment and activate it, you have to install the necessary dependencies below:

- django
- djangorestframework
- numpy
- pandas
- scikit-learn
- joblib
- tensorflow
- keras
- requests
- biopython
- loguru
- sumonet

To install all of the dependencies, write the following line to the terminal: `pip install django djangorestframework numpy pandas scikit-learn joblib tensorflow keras requests biopython loguru sumonet`

## Routes

We prepared different routes for different kind of inputs.

1. UniprotID Prediction

   - For this functionality, use http://127.0.0.1:8000/uniprot-prediction/. This route needs Request Body which consists of uniprot_id(required, string) and lysine_position(required, integer).

   - Example:

   {
   "uniprot_id": "O00566",
   "lysine_position": 20
   }

2. Protein Sequence Prediction

   - For this functionality, use http://127.0.0.1:8000/protein-sequence-prediction/. This route needs Request Body which consists of only protein_seq(required, string).

   - Example:

   {
   "protein_seq": ">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2 MAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKIL..."
   }

3. FASTA File Prediction

   - For this functionality, use http://127.0.0.1:8000/fasta-file-prediction/. This route needs a file which is a fasta file.
