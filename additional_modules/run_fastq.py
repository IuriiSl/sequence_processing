from typing import Union
q_score = {33: 0, 34: 1, 35: 2, 36: 3, 37: 4, 38: 5, 39: 6, 40: 7, 41: 8, 42: 9, 43: 10,
           44: 11, 45: 12, 46: 13, 47: 14, 48: 15, 49: 16, 50: 17, 51: 18, 52: 19, 53: 20,
           54: 21, 55: 22, 56: 23, 57: 24, 58: 25, 59: 26, 60: 27, 61: 28, 62: 29, 63: 30,
           64: 31, 65: 32, 66: 33, 67: 34, 68: 35, 69: 36, 70: 37, 71: 38, 72: 39, 73: 40}
fastq_alphabet = {'A', 'G', 'C', 'T', 'U', 'a', 'g', 'c', 't', 'u'}


def bounds_convert(**kwargs) -> dict:
    """
This function converts bounds in several keyword arguments of different type to dict with lower and upper bounds
    :param kwargs: gc_bounds:int, float or tuple; length_bounds: int or tuple, quality_threshold: int or float
    :return: dictionary containing quality threshold and lower and upper bounds of gc-content and length
    """
    dict_of_function = kwargs
    dict_of_bound = {'quality_threshold': dict_of_function['quality_threshold']}
    if type(dict_of_function['gc_bounds']) is int or type(dict_of_function['gc_bounds']) is float:
        dict_of_bound['lower_gc_bound'] = 0
        dict_of_bound['upper_gc_bound'] = dict_of_function['gc_bounds']
    else:
        dict_of_bound['lower_gc_bound'], dict_of_bound['upper_gc_bound'] = dict_of_function['gc_bounds']
    if type(dict_of_function['length_bounds']) is int or type(dict_of_function['length_bounds']) is float:
        dict_of_bound['lower_length_bound'] = 0
        dict_of_bound['upper_length_bound'] = dict_of_function['length_bounds']
    else:
        dict_of_bound['lower_length_bound'], dict_of_bound['upper_length_bound'] = dict_of_function['length_bounds']
    return dict_of_bound


def gc_filter(seqs: dict, lower_gc_bound: Union[int, float] = 0, upper_gc_bound: Union[int, float] = 80) -> dict:
    """
This function filters reads based on the GC-content
    :param seqs: dictionary consisting of fastq sequences
    :param lower_gc_bound: int or float
    :param upper_gc_bound: int or float
    :return: dictionary where sequences matching filter replaced by T
    """
    seqs_gc_filtered = {}
    for keys, values in seqs.items():
        gc_sum = 0
        sequence, quality = values
        for nucl in sequence:
            if nucl == 'G' or nucl == 'C':
                gc_sum += 1
        gc_percent = gc_sum / len(sequence) * 100
        if lower_gc_bound <= gc_percent <= upper_gc_bound:
            seqs_gc_filtered[keys] = 'T'
        else:
            seqs_gc_filtered[keys] = values
    return seqs_gc_filtered


def length_filter(seqs: dict, lower_length_bound: int = 0, upper_length_bound: int = 2**32) -> dict:
    """
This function filters reads based on the length of their nucleotide sequence
    :param seqs: dictionary consisting of fastq sequences
    :param lower_length_bound: lower bound of interval for filtering.
    :param upper_length_bound: upper bound of interval for filtering.
    :return: dictionary where sequences matching filter replaced by T
    """
    seqs_len_filtered = {}
    for keys, values in seqs.items():
        sequence, quality = values
        if lower_length_bound <= len(sequence) <= upper_length_bound:
            seqs_len_filtered[keys] = 'T'
        else:
            seqs_len_filtered[keys] = values
    return seqs_len_filtered


def quality_filter(seqs: dict, quality_threshold: Union[int, float] = 0) -> dict:
    """
This function filters reads based on phred-score
    :param seqs: dictionary consisting of fastq sequences
    :param quality_threshold: threshold value of average read quality for filtering. Default = 0
    :return: dictionary where sequences matching filter replaced by T
    """
    seqs_qual_filtered = {}
    for keys, values in seqs.items():
        q_score_sum = 0
        sequence, quality = values
        for symbol in quality:
            q_score_sum += q_score[ord(symbol)]
        mean_quality = q_score_sum / len(quality)
        if mean_quality > quality_threshold:
            seqs_qual_filtered[keys] = 'T'
        else:
            seqs_qual_filtered[keys] = values
    return seqs_qual_filtered
