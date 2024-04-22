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

After you set up the virtual environment and activate it, you have to install the necessary dependencies by following the steps:

To install all of the dependencies, write the following line to the terminal: `pip install -r requirements.txt`

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

## API Documentation

Detailed explanation for request type and response type can be seen from the following link: https://docs.google.com/document/d/1vUGhEV5oa-8WZYv_gfbh6oWbaFuhxcrQgFWU-8qHF_0/edit?usp=sharing

## Tests

### Unit Tests

This Python script contains a series of test cases for a Django REST API, focusing on validation of Uniprot ID handling. It uses Django's testing framework along with the Django REST Framework's APIClient for sending requests and evaluating responses. The tests ensure proper response codes and error messages are returned under various conditions, such as empty or incorrect Uniprot IDs.

Key Features:

Test cases and expected outputs can be accessed from the folllowing link: https://docs.google.com/document/d/1zeTQoopuXNkh80EMZUu0ro7CFHZPM1Z6nzKmzTIqxp8/edit?usp=sharing

Requirements:

Django
Django REST Framework
(Optional) Additional dependencies may be required, such as ddt for data-driven testing.
Usage:

The script is designed to be run as part of Django's test suite via the "python manage.py test sumonetWeb" command.
This script is essential for developers working with this specific REST API, ensuring robustness and reliability through comprehensive testing.

### Postman Tests

Postive and negative path tests are automated with Postman. Tests can be accessible from https://www.postman.com/sumonet/workspace/my-workspace/collection/34442047-794ed123-d395-4e6a-b7cc-b41694e49e0a?action=share&source=copy-link&creator=34442047

## Contact

If you have any questions or problem, please reach [eustun@sabanciuniv.edu](mailto:eustun@sabanciuniv.edu) or [aturel@sabanciuniv.edu](mailto:aturel@sabanciuniv.edu)
