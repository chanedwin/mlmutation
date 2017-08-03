from keras.models import *

from ANNgeneratematrixes import *

import os
from h5py import *

def main(paths):
    input, model, reference, output, threshold = load_references(paths)
    generate_matrixes(input, model, reference, output, threshold)


def predone_train_neural_net(length_of_caller_outputs, my_x_dataset, model):
    X_1, X_2, X_3, X_4, X_5 = prep_input_samples(length_of_caller_outputs, my_x_dataset)
    final_predictions = model.predict([X_1, X_2, X_3, X_4, X_5])
    return final_predictions

# INPUT MUST BE A DIRECTORY!

def generate_matrixes(input, model, reference, outputpath, threshold):
    length_of_caller_outputs, list_of_called_samples, vcf_list = generate_input(input, reference)
    my_x_dataset = check_predicted_with_truth(list_of_called_samples)
    calculated_prediction_actual = predone_train_neural_net(length_of_caller_outputs, my_x_dataset, model)
    list_of_records = create_list_of_records(calculated_prediction_actual, vcf_list, threshold)
    vcf_reader = vcf.Reader(filename=original_vcf_reader)
    vcf_writer = vcf.Writer(open(outputpath + "/truevcf.vcf", 'w'), vcf_reader)
    for record in list_of_records:
        vcf_writer.write_record(record)


def create_list_of_records(calculated_prediction_actual, vcf_dictionary, threshold):
    list_of_records = []
    if len(calculated_prediction_actual) != len(vcf_dictionary) :
        raise Exception("vcf list should be same length as calculated predictions")
    for i in range(len(calculated_prediction_actual)):
        if calculated_prediction_actual[i] >= threshold:
            vcf_dictionary[i].INFO['NN_prediction'] = calculated_prediction_actual[i][0]
            list_of_records.append(vcf_dictionary[i])
    return list_of_records


def load_references(paths):
    input = paths['input']
    os.chdir(paths['model'])
    model = load_model(paths['name'])
    os.chdir(paths['input'])
    reference = paths['reference']
    output = paths['output']
    threshold = float(paths['threshold'])
    return input, model, reference, output, threshold


def prep_input_samples(array_sizes, x_training_data):
    count = 0
    X_fb = np.array(map(lambda x: x[count:array_sizes[0]], x_training_data))
    count += array_sizes[0]
    X_hc = np.array(map(lambda x: x[count:count + array_sizes[1]], x_training_data))
    count += array_sizes[1]
    X_ug = np.array(map(lambda x: x[count:count + array_sizes[2]], x_training_data))
    count += array_sizes[2]
    X_pindel = np.array(map(lambda x: x[count:count + array_sizes[3]], x_training_data))
    count += array_sizes[3]
    X_st = np.array(map(lambda x: x[count:count + array_sizes[4]], x_training_data))
    count += array_sizes[4]
    print "x data set is", x_training_data[:100]
    return X_fb, X_hc, X_ug, X_pindel, X_st


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-r', '--reference', help="give directories with files")
    parser.add_argument('-m', '--model', help="give directories with files")
    parser.add_argument('-i', '--input', help="give directories with files")
    parser.add_argument('-o', '--output', help="give directories with files")
    parser.add_argument('-n', '--name', help="give directories with files")
    parser.add_argument('-t', '--threshold', help="give directories with files")
    input_path = parser.parse_args()
    input_path = vars(input_path)
    main(input_path)
