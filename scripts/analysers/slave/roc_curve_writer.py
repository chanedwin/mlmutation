import argparse
import cPickle as pickle
import sys
import csv
import numpy as np
import vcf
from imblearn.over_sampling import SMOTE
from keras.layers import Dense, Dropout, Activation, Merge
from keras.layers.advanced_activations import LeakyReLU
from keras.layers.normalization import BatchNormalization
from keras.models import Sequential
from keras.callbacks import *
from keras.optimizers import RMSprop
from keras.models import load_model

from sklearn.metrics import *
from sklearn.model_selection import train_test_split

STEP_INCREMENT = 10

RECURSION_LIMIT = 0.0002

VERBOSE = 1

vcf_file_name = "/ANN/truevcf.vcf"

roc_file_name = "/roc.csv"

keras_model_name = "/ANN/model"

model_truth_name = "/ANN/modeltruths.txt"

model_predictions_name = "/ANN/modelpredictions.txt"

seed = 1337
np.random.seed(seed)  # for reproducibility

original_vcf_reader = "/data/backup/metacaller/stage/data/version6.3a/hc.vcf.normalisedtrain.vcf"


###6.9.16 First attempt at building ANN

### OUTPUT GOAL : a dictionary of lists, where the dictionary keys are mutations, and the list contains a matrix containing information of all six callers

### To create the keys for the dictionary, all VCF files are first unioned to find a set of unique dictionary keys

### Secondly, a check matrix is created such that it is known which mutations are called in which caller, and which are not (basic 1,0 matrix)

### Finally, for each record in each file, a checkscript first checks if that mutation is inside the caller, if it is, append data, if not, append all 0s

# this method takes in a path and returns training matrixes for the ANN
# The path should contain n caller vcf files and 1 truth file
# vcf files should be labelled with vcf and truth file should be labelled with truth
# no other file should be present in the folder
def main_gather_input_execute_prep_output(array_sizes, dict_of_truth_input, fullmatrix_sample, fullmatrix_truth,
                                          list_of_samples_input, save_location, vcf_dictionary):
    calculated_prediction_actual, calculated_truth_actual = train_neural_net(20, 20, fullmatrix_sample,
                                                                             fullmatrix_truth,
                                                                             save_location, array_sizes)
    get_all_relevant_scores(calculated_prediction_actual, calculated_truth_actual, dict_of_truth_input,
                            list_of_samples_input, vcf_dictionary, save_location)


def count_false_negative(calculated_prediction_actual, calculated_truth_actual):
    count_false_negative = 0
    for i in range(len(calculated_prediction_actual)):
        if calculated_prediction_actual[i] == 0 and calculated_truth_actual[i] == 1:
            count_false_negative += 1
    return count_false_negative


def get_all_relevant_scores(calculated_prediction_actual, calculated_truth_actual, dict_of_truth_input,
                            list_of_samples_input, vcf_list, outputpath):
    list_of_x_variables = []
    list_of_precision_scores = []
    list_of_recall_scores = []
    for i in np.linspace(0,1,101):
        precision_score,recall_score,f1_score = get_scores(calculated_prediction_actual, calculated_truth_actual, i, list_of_samples_input,
                   dict_of_truth_input, VERBOSE)
        list_of_x_variables.append(i)
        list_of_precision_scores.append(precision_score)
        list_of_recall_scores.append(recall_score)
    promise_vcf_file(list_of_x_variables, list_of_precision_scores, list_of_recall_scores, outputpath)


def promise_vcf_file(list_of_x_variables, list_of_precision_scores, list_of_recall_scores, outputpath):
    with open(outputpath + roc_file_name, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(list_of_x_variables)):
            spamwriter.writerow([list_of_x_variables[i]]+[list_of_precision_scores[i]]+[list_of_recall_scores[i]])


def load_references(input_paths):
    input_paths = vars(input_paths)
    fullmatrix_sample = np.load(input_paths['input'][0])
    fullmatrix_truth = np.load(input_paths['input'][1])
    with open(input_paths['input'][3], 'rb') as fp1:
        list_of_samples_input = pickle.load(fp1)
    with open(input_paths['input'][4], 'rb') as fp2:
        dict_of_truth_input = pickle.load(fp2)
    array_sizes = np.load(input_paths['input'][5])
    with open(input_paths['input'][6], 'rb') as fp3:
        vcf_dictionary = pickle.load(fp3)
    orig_stdout = sys.stdout
    f = file(str(input_paths['input'][3]) + '.txt', 'w')
    sys.stdout = f
    return array_sizes, dict_of_truth_input, fullmatrix_sample, fullmatrix_truth, list_of_samples_input, input_paths, vcf_dictionary


def remove_duplicated_false_negative(prediction_list, truth_list, false_negatives):
    count = 0
    removal_list = []
    for i in range(len(prediction_list)-1,-1,-1):
        if count == false_negatives :
            break
        if prediction_list[i] == 0 and truth_list[i] == 1 :
            removal_list.insert(0,i)
            count += 1
    for index in removal_list :
        prediction_list.pop(index)
        truth_list.pop(index)
    return prediction_list,truth_list


def get_scores(actual_predictions, actual_truth, value, sample_list, truth_dictionary, verbose=0):
    temp_actual_truth = list(actual_truth)
    prediction = []
    for item in actual_predictions:
        if item > value:
            prediction.append(1)
        else:
            prediction.append(0)
    false_negatives = count_false_negative(actual_predictions,actual_truth)
    finalpredictionnumbers, finaltruthnumbers = add_negative_data(sample_list, truth_dictionary,
                                                                  prediction, temp_actual_truth)
    finalpredictionnumbers, finaltruthnumbers = remove_duplicated_false_negative(finalpredictionnumbers,
                                                                                 finaltruthnumbers, false_negatives)
    final_f1_score = f1_score(finaltruthnumbers, finalpredictionnumbers)
    if verbose:
        print_scores(actual_truth, final_f1_score, finalpredictionnumbers, finaltruthnumbers, prediction, value)
    return precision_score(finaltruthnumbers, finalpredictionnumbers), recall_score(finaltruthnumbers, finalpredictionnumbers), final_f1_score


def print_scores(actual_truth, final_f1_score, finalpredictionnumbers, finaltruthnumbers, prediction, value):
    false_positive_before_adjust, true_negative_before_adjust = perf_measure(actual_truth, prediction)
    print "false positive is :", false_positive_before_adjust
    print "true negative is :", true_negative_before_adjust
    print "precision score is :", precision_score(actual_truth, prediction)
    print "recall score is :", recall_score(actual_truth, prediction)
    print "F1 score is : ", f1_score(actual_truth, prediction)
    final_false_positive, final_true_negative = perf_measure(finaltruthnumbers, finalpredictionnumbers)
    print "final false positive is :", final_false_positive
    print "final true negative is :", final_true_negative
    print "final precision score is :", precision_score(finaltruthnumbers, finalpredictionnumbers)
    print "final recall score is :", recall_score(finaltruthnumbers, finalpredictionnumbers)
    print "threshold is", value
    print "final F1 score is : ", final_f1_score


def print_details_of_score(actual_truth, finalpredictionnumbers, finaltruthnumbers, prediction):
    false_positive_before_adjust, true_negative_before_adjust = perf_measure(actual_truth, prediction)
    print "false positive is :", false_positive_before_adjust
    print "true negative is :", true_negative_before_adjust
    print "precision score is :", precision_score(actual_truth, prediction)
    print "recall score is :", recall_score(actual_truth, prediction)
    print "F1 score is : ", f1_score(actual_truth, prediction)
    final_false_positive, final_true_negative = perf_measure(finaltruthnumbers, finalpredictionnumbers)
    print "final false positive is :", final_false_positive
    print "final true negative is :", final_true_negative
    print "final precision score is :", precision_score(finaltruthnumbers, finalpredictionnumbers)
    print "final recall score is :", recall_score(finaltruthnumbers, finalpredictionnumbers)


def add_negative_data(list_of_samples, dict_of_truth, array_of_predicted, array_of_truth):
    dict_of_samples = generate_sample_dictionary(array_of_predicted, list_of_samples)
    list_of_truth = generate_list_of_truth(dict_of_truth)
    new_array_of_predicted = list(array_of_predicted)
    new_array_of_truth = list(array_of_truth)
    original_length = len(new_array_of_predicted)
    for item in list_of_truth:
        fillnegative(item, dict_of_samples, new_array_of_predicted, new_array_of_truth)
    print "number of false data samples are", (len(new_array_of_predicted) - original_length)
    return new_array_of_predicted, new_array_of_truth


def generate_list_of_truth(dict_of_truth):
    list_of_truth = []
    for key in dict_of_truth:
        mytuple = dict_of_truth[key]
        temptuple =[]
        for item in mytuple:
            temptuple.append(item)
        list_of_truth.append([key[0], key[1], key[2], temptuple])
    return list_of_truth

def generate_sample_dictionary(array_of_predicted, list_of_samples):
    dict_of_samples = {}
    for i in range(len(list_of_samples)):
        item = list_of_samples[i]
        if array_of_predicted[i] == 0:
            continue
        new_key = (item[0][0], item[0][1], item[0][2])
        new_value = item[0][3]
        if new_key not in dict_of_samples:
            dict_of_samples[new_key] = new_value
        else:
            dict_of_samples[new_key] = list(dict_of_samples[new_key])
            dict_of_samples[new_key].extend(new_value)
            dict_of_samples[new_key] = tuple(dict_of_samples[new_key])
            # print dict_of_samples[new_key]
    return dict_of_samples


def perf_measure(y_actual, y_hat):
    true_positive = 0
    false_positive = 0
    false_negative = 0
    true_negative = 0

    for i in range(len(y_hat)):
        if y_actual[i] == 1 and y_hat[i] == 1:
            true_positive += 1
    for i in range(len(y_hat)):
        if y_hat[i] == 1 and y_actual[i] == 0:
            false_positive += 1
    for i in range(len(y_hat)):
        if y_actual[i] == 1 and y_hat[i] == 0:
            false_negative += 1
    for i in range(len(y_hat)):
        if y_hat[i] == 0 and y_actual[i] == 0:
            true_negative += 1

    print "true positives :", true_positive
    print "false positives :", false_positive
    print "false negatives :", false_negative
    print "true negatives :", true_negative

    true_positive = float(true_positive)
    false_positive = float(false_positive)
    false_negative = float(false_negative)
    if false_positive == 0 and true_positive == 0:
        false_positive_rate = 0
    else:
        false_positive_rate = false_positive / (false_positive + true_positive)
    if false_negative == 0 and true_positive == 0:
        true_negative_rate = 0
    else:
        true_negative_rate = false_negative / (false_negative + true_positive)

    return false_positive_rate, true_negative_rate


def fillnegative(tuple1, sampledict, arrayofsamples, arrayoftruths):
    tuple2 = (tuple1[0], tuple1[1], tuple1[2])
    if tuple2 in sampledict:
        for ALT in tuple1[3]:
            if ALT in sampledict[tuple2]:
                return
    arrayofsamples.append(0)
    arrayoftruths.append(1)



def train_neural_net(mybatch_size, mynb_epoch, myX_train, myy_train, location, array_sizes):
    fb_size, hc_size, ug_size, pindel_size, st_size = get_sizes(array_sizes)
    X_resampled, y_resampled = do_smote_resampling(myX_train, myy_train)
    X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled,
                                                        test_size=0.33, random_state=seed)
    X_fb, X_hc, X_ug, X_pindel, X_st = prep_input_samples(array_sizes, X_train)
    X_fb_test, X_hc_test, X_ug_test, X_pindel_test, X_st_test = prep_input_samples(array_sizes, X_test)
    batch_size = mybatch_size
    nb_epoch = mynb_epoch

    fb_branch = Sequential()
    develop_first_layer_matrixes(fb_branch, fb_size)

    hc_branch = Sequential()
    develop_first_layer_matrixes(hc_branch, hc_size)

    ug_branch = Sequential()
    develop_first_layer_matrixes(ug_branch, ug_size)

    pindel_branch = Sequential()
    develop_first_layer_matrixes(pindel_branch, pindel_size)

    st_branch = Sequential()
    develop_first_layer_matrixes(st_branch, st_size)

    final_model = Sequential()
    final_model.add(Merge([fb_branch, hc_branch, ug_branch, pindel_branch, st_branch], mode='concat', concat_axis=1))
    final_model.add(LeakyReLU(alpha=0.05))
    final_model.add(Dense(24, activation='linear'))
    final_model.add(LeakyReLU(alpha=0.05))
    final_model.add(Dense(1, activation='linear'))
    final_model.add(Activation('sigmoid'))
    print (final_model.summary())
    rmsprop = RMSprop(lr=0.000001, rho=0.9, epsilon=1e-08, decay=0.0)
    final_model.compile(loss='binary_crossentropy',
                        optimizer=rmsprop,
                        metrics=['accuracy'])
    filepath = location + "/best_weights.hdf5"
    checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')
    callbacks_list = [checkpoint]
    final_model.fit([X_fb, X_hc, X_ug, X_pindel, X_st], y_train, batch_size=batch_size, nb_epoch=nb_epoch,
                    validation_split=0.2, verbose=2, callbacks= callbacks_list)
    final_model = load_model(location + "/best_weights.hdf5")
    scores = final_model.evaluate([X_fb_test, X_hc_test, X_ug_test, X_pindel_test, X_st_test], y_test)
    print scores
    X_fb_pred, X_hc_pred, X_ug_pred, X_pindel_pred, X_st_pred = prep_input_samples(array_sizes, myX_train)
    final_prediction_array_probabilities = final_model.predict([X_fb_pred, X_hc_pred, X_ug_pred, X_pindel_pred,
                                                                X_st_pred])
    final_prediction_array_probabilities = np.squeeze(final_prediction_array_probabilities)

    save_model_details(final_model, final_prediction_array_probabilities, myy_train, location)

    return final_prediction_array_probabilities, myy_train


def do_smote_resampling(myX_train, myy_train):
    sm = SMOTE(kind='regular')
    where_are_NaNs = np.isnan(myX_train)
    myX_train[where_are_NaNs] = 0
    X_resampled, y_resampled = sm.fit_sample(myX_train, myy_train)
    return X_resampled, y_resampled


def save_model_details(final_model, save_model_probabilities, trutharray, location):
    name1 = location + model_predictions_name
    name2 = location + model_truth_name
    name3 = location + keras_model_name
    np.save(name1, save_model_probabilities)
    np.save(name2, trutharray)
    final_model.save(name3)


def develop_first_layer_matrixes(neural_net_branch, branch_size):
    neural_net_branch.add(BatchNormalization(input_shape=(branch_size,), axis=1))
    neural_net_branch.add(Dense(24, activation='linear'))
    neural_net_branch.add(LeakyReLU(alpha=0.05))
    neural_net_branch.add(Dense(24, activation='linear'))
    neural_net_branch.add(LeakyReLU(alpha=0.05))
    neural_net_branch.add(Dropout(0.2))
    neural_net_branch.add(Dense(24, activation='linear'))
    neural_net_branch.add(LeakyReLU(alpha=0.05))
    neural_net_branch.add(Dropout(0.2))
    neural_net_branch.add(Dense(6, activation='linear'))
    neural_net_branch.add(LeakyReLU(alpha=0.05))


def get_sizes(array_sizes):
    fb_size = array_sizes[0]
    hc_size = array_sizes[1]
    ug_size = array_sizes[2]
    pindel_size = array_sizes[3]
    st_size = array_sizes[4]
    return fb_size, hc_size, ug_size, pindel_size, st_size


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
    return X_fb, X_hc, X_ug, X_pindel, X_st


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="train neural net")
    parser.add_argument('-i', '--input', help="give directories with files", nargs='+')
    input_path = parser.parse_args()
    array_sizes, dict_of_truth_input, fullmatrix_sample, fullmatrix_truth, \
    list_of_samples_input, paths, vcf_dictionary = load_references(input_path)
    main_gather_input_execute_prep_output(array_sizes, dict_of_truth_input, fullmatrix_sample, fullmatrix_truth,
                                          list_of_samples_input, paths, vcf_dictionary)
