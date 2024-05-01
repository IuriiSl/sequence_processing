import os
import pytest
import numpy as np

from bioinformatics_toolkit import DNASequence, AminoAcidSequence, filter_fastq
from custom_random_forest import RandomForestClassifierCustom
from bio_files_processor import convert_multiline_fasta_to_oneline, OpenFasta, FastaRecord


def test_dna_len():
    inp = 'AUGC'
    target = 4
    result = DNASequence(inp).__len__()
    assert target == result


def test_dna_transcribe():
    inp = 'ATGCG'
    target = 'UACGC'
    result = str(DNASequence(inp).transcribe())
    assert target == result


def test_molecular_weight():
    inp = 'ANCQE'
    target = 635
    result = AminoAcidSequence(inp).count_molecular_weight()
    assert target == result


def test_random_forest_fit():
    inp_1 = np.array([0.09673332, -0.38946408,  0.25724783,  0.02541841,  0.25841846,
                      0.44867822, -0.42542166, -1.28446403,  0.33217226, -0.47122772,
                      0.18490433,  2.02586933,  0.64312296, -1.6976305,  0.76982948,
                      -0.0491514,  1.45639268,  0.39745317, -0.91571756,  1.2780103])
    inp_2 = np.array([1, 0, 1])
    with pytest.raises(IndexError):
        RandomForestClassifierCustom().fit(inp_1, inp_2, n_jobs=2)


@pytest.fixture
def tmp_file():
    file_path = 'tmp.fasta'
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)


def test_convert_multiline_fasta(tmp_file):
    inp = './data/example_multiline_fasta.txt'
    convert_multiline_fasta_to_oneline(inp, tmp_file)
    assert os.path.exists(tmp_file)


def test_write_file_filter_fastq(tmp_file):
    inp = './data/example_fastq.fastq'
    filter_fastq(inp, tmp_file)
    assert os.path.exists(tmp_file)


def test_filter_fastq(tmp_file):
    inp = './data/example_fastq.fastq'
    target = 4
    filter_fastq(inp, tmp_file, upper_length_bound=2)
    with open(tmp_file, 'r') as file:
        line_count = 0
        for _ in file.readlines():
            line_count += 1
    assert target == line_count


def test_open_fasta_output_type():
    inp = './data/example_fasta.fasta'
    target = FastaRecord
    with OpenFasta(inp) as fasta_file:
        result = type(fasta_file.read_record())
    assert target == result
