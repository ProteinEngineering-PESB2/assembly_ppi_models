import pandas as pd
import os
from scipy.fft import fft
import numpy as np
import argparse

#Apply Fast Fourier transform
def apply_FFT(row, number_sample):
    T = 1.0 / float(number_sample)
    x = np.linspace(0.0, number_sample * T, number_sample)
    yf = fft(row)
    xf = np.linspace(0.0, 1.0 / (2.0 * T), number_sample // 2)
    yf = np.abs(yf[0:number_sample // 2])
    return xf, yf

#encoding dataset using FFT
def encoding_using_fft(dataset):

    matrix_response = []
    for i in range(len(dataset)):
        row = [dataset[key][i] for key in dataset.columns]
        domain, encoding = apply_FFT(row, len(row))
        matrix_response.append(encoding)

    header = ["p_{}".format(i) for i in range(len(matrix_response[0]))]
    df_export = pd.DataFrame(matrix_response, columns=header)
    return df_export

#encoding sequences using property sequences
def encoding_sequence(sequence, dict_encoder, value_encoding):

    array_response = []
    sequence = sequence.upper()
    for residue in sequence:
        try:
            value_encode = value_encoding[dict_encoder[residue]]
            array_response.append(value_encode)
        except:
            pass
    return array_response

#encoding dataset using protein sequences
def encoding_dataset(value_encoding, dataset_to_encode, dict_encoder):

    matrix_encoder = []
    length_values = []
    response_data = []
    for i in range(len(dataset_to_encode)):

        #name columns define in dataset
        sequence = dataset_to_encode['seq'][i]
        response = dataset_to_encode['response'][i]
        vector_encoder =  encoding_sequence(sequence, dict_encoder, value_encoding)
        matrix_encoder.append(vector_encoder)
        length_values.append(len(vector_encoder))
        response_data.append(response)

    #apply zero-padding
    max_lenght = max(length_values)
    for i in range(len(matrix_encoder)):
        length_data = matrix_encoder[i]
        for j in range(len(length_data), max_lenght):
            matrix_encoder[i].append(0)

    header = ["p_{}".format(i) for i in range(0, max_lenght)]
    df_data = pd.DataFrame(matrix_encoder, columns=header)
    df_data['response'] = response_data
    return df_data

#console params sections, use argparse to process input console line
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dataset", help="input dataset in csv format to encoding sequences", required=True)
parser.add_argument("-e", "--encoders", help="input dataset with AAIndex encoders", required=True)
parser.add_argument("-p", "--path", help="Path export to save encoding results", required=True)
parser.add_argument("-t", "--task", help="Name objective task", required=True)
parser.add_argument("-c", "--combine", help="Name of combination, example: properties_fft, properties", required=True)
parser.add_argument("-f", "--fft", help="Use FFT 1, 0 not use FFT", choices=[0, 1], default=0, type=int)
parser.add_argument("-i", "--property_id", help="AAIndex property ID", default='ANDN920101')

#compile
args = parser.parse_args()

#get values in variables
dataset = pd.read_csv(args.dataset)
dataset_encoders = pd.read_csv(args.encoders)
path_export = args.path
task_to_solve = args.task
combination_to_use = args.combine
use_fft = args.fft
property_id = args.property_id

#support data
dict_encoder = {}
for i in range(len(dataset_encoders['residue'])):
    dict_encoder.update({dataset_encoders['residue'][i]:i})

print("Creating directory")
command = "mkdir {}{}_{}".format(path_export, task_to_solve, combination_to_use)
print(command)
os.system(command)

print("Start encoding with property: ", property_id)

response_encoding = encoding_dataset(dataset_encoders[property_id], dataset, dict_encoder)

if use_fft == 1:
    data_to_fft = response_encoding.drop(columns=['response'])
    response_fft = encoding_using_fft(data_to_fft)
    response_fft['response'] = response_encoding['response']
    response_encoding = response_fft

command = "mkdir {}{}_{}\\{}".format(path_export, task_to_solve, combination_to_use, property_id)
os.system(command)
name_export = "{}{}_{}\\{}\\encoding_data.csv".format(path_export, task_to_solve, combination_to_use, property_id)
response_encoding.to_csv(name_export, index=False)
