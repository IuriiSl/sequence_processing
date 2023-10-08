from typing import Union
def run_aminoacid_seq(sequence: str, function: str = 'summary', record_type: int = 1, percent: bool = False):
    """
    Performs the following list of operations:
    - count - counts number of amino acid
    - translate - converts one record type to another
    - determine_charge - counts number or percent of amino acid with different charges
    - determine_polarity - counts number or percent of amino acid with different polarity
    - convert_amino_acid_seq_to_dna - takes an amino acid sequence as input
        and returns the optimal DNA sequence for E.coli
    - count_possible_number_of_disulfide_bonds - counting the number of possible combinations of two different cysteines
        to form a disulfide bond
    - count_molecular_weight - takes an amino acid sequence as input
        and returns the molecular weight of the protein
    - summary - returns results of all functions (default)

    Arguments:
    - sequence:str - sequence for function
    - function:str - name of the function you need to perform. You can use: 'count', 'translate', 'summary'(default)
    - record_type:int - record type of your sequence. You can use 1(default) for single letter type
        or 3 for three letter type
    - percent: bool - for determine_charge and determine_polarity shows result in percent (True) or in number (False).
        Default - False

    Return:
    - count - int
    - translate - str
    - determine_charge - dict
    - determine_polarity - dict
    - convert_amino_acid_seq_to_dna - str
    - count_possible_number_of_disulfide_bonds - int
    - count_molecular_weight - int
    - summary - dict

    """
    import additional_modules.run_aminoacid_seq
    if record_type == 1:
        for amino_acid in sequence:
            if amino_acid not in additional_modules.run_aminoacid_seq.aminoacid_dict.values():
                raise ValueError(amino_acid + ' in your sequence is not amino acid')
    elif record_type == 3:
        tuple_sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        for amino_acid in tuple_sequence:
            if amino_acid not in additional_modules.run_aminoacid_seq.aminoacid_dict.keys():
                raise ValueError(amino_acid + ' is not amino acid')
        sequence = ''.join(tuple_sequence)
    if function == 'count':
        return additional_modules.run_aminoacid_seq.count(sequence, record_type)
    elif function == 'translate':
        return translate(sequence, record_type)
    elif function == 'summary':
        return additional_modules.run_aminoacid_seq.summary(sequence, record_type, percent)
    elif function == 'determine_polarity':
        if record_type == 3:
            return additional_modules.run_aminoacid_seq.determine_polarity(translate(sequence, record_type), percent)
        return additional_modules.run_aminoacid_seq.determine_polarity(sequence, percent)
    elif function == 'determine_charge':
        if record_type == 3:
            return additional_modules.run_aminoacid_seq.determine_charge(translate(sequence, record_type), percent)
        return additional_modules.run_aminoacid_seq.determine_charge(sequence, percent)
    elif function == 'count_possible_number_of_disulfide_bonds':
        if record_type == 3:
            return additional_modules.run_aminoacid_seq.count_possible_number_of_disulfide_bonds(translate(sequence, record_type))
        return additional_modules.run_aminoacid_seq.count_possible_number_of_disulfide_bonds(sequence)
    elif function == 'count_molecular_weight':
        if record_type == 3:
            return additional_modules.run_aminoacid_seq.count_molecular_weight(translate(sequence, record_type))
        return additional_modules.run_aminoacid_seq.count_molecular_weight(sequence)
    elif function == 'convert_amino_acid_seq_to_dna':
        if record_type == 3:
            return additional_modules.run_aminoacid_seq.convert_amino_acid_seq_to_dna(translate(sequence, record_type))
        return additional_modules.run_aminoacid_seq.convert_amino_acid_seq_to_dna(sequence)


def run_dna_rna_tools(*args) -> Union[str, list]:
    """
Performs the following procedures with nucleotide sequences:
    transcribe — returns the transcribed sequence
    reverse — returns the reversed sequence
    complement — returns the complementary sequence
    reverse_complement — returns the reverse complementary sequence
    :param args: str or list of sequence of nucleotides and name of function
    :return:str or list of processed sequence
    """
    import additional_modules.run_dna_rna_tools
    *sequence, function = args
    if type(sequence) == list:
        for nucleotides in sequence:
            for nucleotide in nucleotides:
                if nucleotide in additional_modules.run_dna_rna_tools.dna_rna_complement_dict.keys() and nucleotides in additional_modules.run_dna_rna_tools.dna_rna_complement_dict.values():
                    raise ValueError('use correct sequence')
                if ('U' or 'u') in nucleotides and (function == 'transcribe' or function == 'complement_transcribe'):
                    raise ValueError('use correct sequence')
    else:
        for nucleotides in sequence:
            if nucleotides in additional_modules.run_dna_rna_tools.dna_rna_complement_dict.keys() and nucleotides in additional_modules.run_dna_rna_tools.dna_rna_complement_dict.values():
                raise ValueError('use correct sequence')
            if ('U' or 'u') in nucleotides and (function == 'transcribe' or function == 'complement_transcribe'):
                raise ValueError('use correct sequence')
    if function == 'reverse':
        return additional_modules.run_dna_rna_tools.reverse(sequence)
    elif function == 'transcribe':
        return additional_modules.run_dna_rna_tools.transcribe(sequence)
    elif function == 'complement':
        return additional_modules.run_dna_rna_tools.complement(sequence)
    elif function == 'reverse_complement':
        return additional_modules.run_dna_rna_tools.reverse_complement(sequence)


def run_fastq(seqs: dict, gc_bounds: Union[int, float, tuple] = (0, 100),
              length_bounds: Union[int, float, tuple] = (0, 2**32), quality_threshold: int = 0) -> dict:
    """
This function filters reads based on the length of their nucleotide sequence, GC-content and phred33-score
    :param seqs: dictionary consisting of fastq sequences
    :param gc_bounds: GC-content interval (percentage) for filtering. Default = 0, 100
    :param length_bounds: length interval for filtering. Default = 0, 2**32
    :param quality_threshold: threshold value of average read quality for filtering. Default = 0
    :return: dictionary consisting of filtered fastq sequences matching all filters
    """
    import additional_modules.run_fastq
    for keys, values in seqs.items():
        sequence, quality = values
        for nucleotide in sequence:
            if nucleotide not in additional_modules.run_fastq.fastq_alphabet:
                raise ValueError(f"{nucleotide} in {sequence} is not nucleotide")
    dict_bounds = additional_modules.run_fastq.bounds_convert(gc_bounds=gc_bounds, length_bounds=length_bounds,
                                 quality_threshold=quality_threshold)
    seqs_gc_filtered = additional_modules.run_fastq.gc_filter(seqs, lower_gc_bound=dict_bounds['lower_gc_bound'],
                                 upper_gc_bound=dict_bounds['upper_gc_bound'])
    seqs_len_filtered = additional_modules.run_fastq.length_filter(seqs, lower_length_bound=dict_bounds['lower_length_bound'],
                                      upper_length_bound=dict_bounds['upper_length_bound'])
    seqs_qual_filtered = additional_modules.run_fastq.quality_filter(seqs, quality_threshold=dict_bounds['quality_threshold'])
    seqs_filtered = {}
    for keys, values in seqs.items():
        if seqs_gc_filtered[keys] == 'T' and seqs_len_filtered[keys] == 'T' and seqs_qual_filtered[keys] == 'T':
            seqs_filtered[keys] = values
    return seqs_filtered
