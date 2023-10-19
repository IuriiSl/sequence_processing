from typing import Union
dna_rna_complement_dict = {'A': 'U', 'G': 'C', 'C': 'G', 'T': 'A',
                           'a': 'u', 'g': 'c', 'c': 'g', 't': 'a'}
dna_rna_transcribe = {'A': 'A', 'G': 'G', 'C': 'C', 'T': 'U',
                                 'a': 'a', 'g': 'g', 'c': 'c', 't': 'u'}
DNA_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
RNA_COMPLEMENT = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C',
                  'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}


def reverse(sequence: Union[str, list]) -> Union[str, list]:
    """
Returns the reversed sequence
    :param sequence: str or list of sequence of nucleotides
    :return: str or list of reversed sequence
    """
    reverse_sequence = []
    for nucleotides in sequence:
        sequence_str = ''.join(nucleotides)
        reverse_str = sequence_str[::-1]
        reverse_sequence.append(reverse_str)
    if len(sequence) > 1:
        return reverse_sequence
    reverse_str = ''.join(reverse_sequence)
    return reverse_str


def transcribe(sequence: Union[str, list]) -> Union[str, list]:
    """
Returns the transcribed sequence
    :param sequence: str or list of sequence of nucleotides
    :return: str or list of transcribed sequence
    """
    transcribed_seq = []
    for nucleotides in sequence:
        for nucleotide in nucleotides:
            transcribed_seq.append(dna_rna_transcribe[nucleotide])
    if len(sequence) > 1:
        transcribed_seq = [transcribed_seq[i:i + 3] for i in range(0, len(transcribed_seq), 3)]
        transcribed_seq = [''.join(sub_list) for sub_list in transcribed_seq]
        return transcribed_seq
    transcribed_str = ''.join(transcribed_seq)
    return transcribed_str


def complement(sequence: Union[str, list]) -> Union[str, list]:
    """
Returns the complementary sequence
    :param sequence: str or list of sequence of nucleotides
    :return: str or list of complementary sequence
    """
    if ('U' or 'u') in sequence:
        dict_complement = RNA_COMPLEMENT
    else:
        dict_complement = DNA_COMPLEMENT
    complement_seq = []
    for nucleotides in sequence:
        for nucleotide in nucleotides:
            complement_seq.append(dict_complement[nucleotide])
    if len(sequence) > 1:
        complement_seq = [complement_seq[i:i + 3] for i in range(0, len(complement_seq), 3)]
        complement_seq = [''.join(sub_list) for sub_list in complement_seq]
        return complement_seq
    complement_str = ''.join(complement_seq)
    return complement_str


def reverse_complement(sequence: Union[str, list]) -> Union[str, list]:
    """
Returns the reverse complementary sequence
    :param sequence: str or list of sequence of nucleotides
    :return: str or list of reverse complementary sequence
    """
    if ('U' or 'u') in sequence:
        dict_complement = RNA_COMPLEMENT
    else:
        dict_complement = DNA_COMPLEMENT
    reverse_complement_seq = []
    for nucleotides in sequence:
        for nucleotide in nucleotides[::-1]:
            reverse_complement_seq.append(dict_complement[nucleotide])
    if len(sequence) > 1:
        reverse_complement_seq = [reverse_complement_seq[i:i + 3] for i in range(0, len(reverse_complement_seq), 3)]
        reverse_complement_seq = [''.join(sub_list) for sub_list in reverse_complement_seq]
        return reverse_complement_seq
    reverse_complement_str = ''.join(reverse_complement_seq)
    return reverse_complement_str


