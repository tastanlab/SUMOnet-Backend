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

## Starting the server

To start the server, make sure that your are in the same directory with `manage.py` file. After that, write the following line to the terminal: `python manage.py runserver`.

> [!important]
> After starting the web server, there can be warnings about unapplied migrations. They are not important, you do not need to be distracted by that error.

## Database

This project does not contain any database connection.

## Routes

We prepared different routes for different kind of inputs.

1. UniprotID Prediction

   - For this functionality, use http://127.0.0.1:8000/uniprot-prediction/. This route needs Request Body which consists of uniprot_id(required, string) and lysine_position(not required, integer).

   - If you do not write lysine position, it will give all the possible lysine positions.

   - Example:

   {
   "uniprot_id": "O00566",
   "lysine_position": 20
   }

2. Protein Sequence Prediction

   - For this functionality, use http://127.0.0.1:8000/protein-sequence-prediction/. This route needs Request Body which consists of only protein_seq(required, string).

   - The format is highly important in this route. Please write the requences as it is written in this example: https://rest.uniprot.org/uniprotkb/O00566.fasta. If you are going to add more protein sequence, please leave at least one line space.

   - Example:

   {
   "protein_seq": ">sp|O00566|MPP10_HUMAN U3 small nucleolar ribonucleoprotein protein MPP10 OS=Homo sapiens OX=9606 GN=MPHOSPH10 PE=1 SV=2 MAPQVWRRRTLERCLTEVGKATGRPECFLTIQEGLASKFTSLTKVLYDFNKIL..."
   }

3. FASTA File Prediction

   - For this functionality, use http://127.0.0.1:8000/fasta-file-prediction/. This route needs a file which is a fasta file.

   - The format is highly important in this route. Please write the requences as it is written in this example: https://rest.uniprot.org/uniprotkb/O00566.fasta. If you are going to add more protein sequence, please leave at least one line space.

For more information about the routes, please take a look at sumonetWeb/views.py file. This file contains routes for functionalities.

## Contact

If you have any questions or problem, please reach eustun@sabanciuniv.edu(mailto:eustun@sabanciuniv.edu) or aturel@sabanciuniv.edu(mailto:aturel@sabanciuniv.edu)
