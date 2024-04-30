import sys
import requests
import os
import io
from typing import Union
from Bio import SeqIO, SeqUtils
from abc import ABC, abstractmethod
from dotenv import load_dotenv
from datetime import datetime, timedelta
from bs4 import BeautifulSoup
from dataclasses import dataclass
from typing import List


def length_filter(record, lower_length_bound: int, upper_length_bound: int) -> bool:
    """
    Applies a length filter to a SeqRecord object.
    :param record: The SeqRecord object to be filtered.
    :param lower_length_bound: The lower bound for the length of the sequence, inclusive.
    :param upper_length_bound: The upper bound for the length of the sequence, inclusive.
    :return: True or False: If the length of the sequence in the SeqRecord falls within the specified
        bounds, the True is returned. Otherwise, False is returned.
    """
    return lower_length_bound <= len(record.seq) <= upper_length_bound


def gc_filter(record, lower_gc_bound: float, upper_gc_bound: float) -> bool:
    """
    Applies a GC content filter to a SeqRecord object.
    :param record: The SeqRecord object to be filtered.
    :param lower_gc_bound: The lower bound for the GC content of the sequence, inclusive.
    :param upper_gc_bound: The upper bound for the GC content of the sequence, inclusive.
    :return: True or False: If the GC content of the sequence in the SeqRecord falls within the specified
        bounds, the True is returned. Otherwise, False is returned.
    """
    return lower_gc_bound <= SeqUtils.GC(record.seq) <= upper_gc_bound


def quality_filter(record, quality_threshold: Union[int, float] = 0) -> bool:
    """
    Applies a quality score filter to a SeqRecord object.
    :param record: The SeqRecord object to be filtered.
    :param quality_threshold:
    :return: True or False: True if the SeqRecord object passes the quality filter, False otherwise.
    """
    return (quality_threshold < sum(record.letter_annotations["phred_quality"]) /
            len(record.letter_annotations["phred_quality"]))


def filter_fastq(file_path: str, output_file: str = None,
                 lower_gc_bound: int = 0, upper_gc_bound: int = 30,
                 lower_length_bound: int = 0, upper_length_bound: int = 2**32,
                 quality_threshold: Union[int, float] = 0) -> None:
    """
    Filters a FASTQ file based on GC content, sequence length, and quality score.
    :param file_path: The path to the input FASTQ file.
    :param output_file: The path to the output FASTQ file. If not provided, it will be named
            '{input_file_name}_filter.fastq'.
    :param lower_gc_bound: The lower bound for the GC content of the sequences, inclusive.
            Defaults to 0.
    :param upper_gc_bound: The upper bound for the GC content of the sequences, inclusive.
            Defaults to 30.
    :param lower_length_bound: The lower bound for the length of the sequences, inclusive.
            Defaults to 0.
    :param upper_length_bound: The upper bound for the length of the sequences, inclusive.
            Defaults to 2^32 (approximately 4.3 billion).
    :param quality_threshold: The minimum quality score threshold for the sequences.
            Defaults to 0.
    :return: None: The function writes the filtered FASTQ records to the output file.
    """
    records = SeqIO.parse(file_path, "fastq")
    if not output_file:
        output_filename = os.path.basename(file_path)
        output_file = output_filename.replace('.fastq', '_filter.fastq')
    with open(output_file, "w") as output:
        for record in records:
            if length_filter(record, lower_length_bound, upper_length_bound) and \
                    gc_filter(record, lower_gc_bound, upper_gc_bound) and \
                    quality_filter(record, quality_threshold):
                SeqIO.write(record, output, "fastq")


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, item: int):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class SequenceFunction(BiologicalSequence):
    def __init__(self, seq: str):
        self.seq = seq
        self.length = len(self.seq)

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, item: int) -> str:
        if 0 <= item < len(self.seq):
            return self.seq[item]
        else:
            raise IndexError("Your index is incorrect")

    def __str__(self) -> str:
        return str(self.seq)

    def check_alphabet(self) -> bool:
        return set(self.seq) <= set(self.alphabet)


class NucleicAcidSequence(SequenceFunction):
    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    def complement(self):
        complement_seq = self.seq.translate(str.maketrans(self.complement_dict))
        return type(self)(complement_seq)

    def gc_content(self) -> float:
        gc_count = self.seq.count('C') + self.seq.count('G')
        return gc_count/len(self.seq)*100


class DNASequence(NucleicAcidSequence):
    alphabet = ('A', 'T', 'G', 'C')
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    def transcribe(self):
        transcribed_seq = self.seq.translate(str.maketrans("ATGC", "UACG"))
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence):
    alphabet = ('A', 'U', 'G', 'C')
    complement_dict = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, seq: str) -> None:
        super().__init__(seq)


class AminoAcidSequence(SequenceFunction):
    alphabet = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    amino_acid_weights = {
        'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121,
        'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131,
        'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115,
        'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117
    }

    def __init__(self, seq) -> None:
        super().__init__(seq)

    def count_molecular_weight(self) -> int:
        molecular_weight = sum(self.amino_acid_weights.get(aa, 0) for aa in self.seq)
        return molecular_weight


load_dotenv('Bot_token.env')
TG_API_TOKEN = os.getenv('TG_API_TOKEN')
LINK = 'https://api.telegram.org/bot' + TG_API_TOKEN


def telegram_logger(chat_id: int):
    def decorator(func):
        def inner_function():
            message_link = LINK + '/sendMessage'
            stdout_backup = sys.stdout
            stderr_backup = sys.stderr
            sys.stdout = io.StringIO()
            sys.stderr = io.StringIO()
            start = datetime.now()
            try:
                func()
            except Exception as e:
                stdout_contents = sys.stdout.getvalue()
                stderr_contents = sys.stderr.getvalue()
                file_contents = stdout_contents + stderr_contents
                file_name = func.__name__ + '.log'
                error_type = type(e)
                error_type = str(error_type).split("'")[1]
                error_text = str(e)
                if file_contents:
                    files = {'document': (file_name, io.BytesIO(file_contents.encode()))}
                    requests.post(LINK + '/sendDocument',
                                  params={'chat_id': chat_id, 'caption': f'Function `{func.__name__}` failed with an exception `{error_type}: {error_text}`',
                                          'parse_mode': 'Markdown'},
                                  files=files)
                else:
                    requests.post(message_link,
                                  params={'chat_id': chat_id,
                                          'text': f'Function `{func.__name__}` failed with an exception `{error_type}: {error_text}`',
                                          'parse_mode': 'Markdown'})
                raise e

            else:
                stdout_contents = sys.stdout.getvalue()
                stderr_contents = sys.stderr.getvalue()
                file_contents = stdout_contents + stderr_contents
                file_name = func.__name__ + '.log'
                end = datetime.now()
                execution_time = end - start
                if execution_time < timedelta(days=1):
                    if file_contents:
                        files = {'document': (file_name, io.BytesIO(file_contents.encode()))}
                        requests.post(LINK + '/sendDocument',
                                      params={'chat_id': chat_id, 'caption': f'Function `{func.__name__}` successfully finished in `{end - start}`',
                                              'parse_mode': 'Markdown'},
                                      files=files)
                    else:
                        requests.post(message_link,
                                      params={'chat_id': chat_id,
                                              'text': f'Function `{func.__name__}` successfully finished in `{end - start}`',
                                              'parse_mode': 'Markdown'})
                else:
                    days = execution_time.days
                    hours, remainder = divmod(execution_time.seconds, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    if file_contents:
                        files = {'document': (file_name, io.BytesIO(file_contents.encode()))}
                        requests.post(LINK + '/sendDocument',
                                      params={'chat_id': chat_id,
                                              'caption': f'Function `{func.__name__}` successfully finished in `{days} days, {hours:02}:{minutes:02}:{seconds:02}`',
                                              'parse_mode': 'Markdown'},
                                      files=files)
                    else:
                        requests.post(message_link,
                                      params={'chat_id': chat_id,
                                              'text': f'Function `{func.__name__}` successfully finished in `{days} days, {hours:02}:{minutes:02}:{seconds:02}`',
                                              'parse_mode': 'Markdown'})
            finally:
                sys.stdout = stdout_backup
                sys.stderr = stderr_backup
        return inner_function
    return decorator


@dataclass
class GenscanOutput:
    """
    Represents the output of the run_genscan function.
    Attributes:
        status (int): The status of the run_genscan output.
        cds_list (List[str]): A list of protein sequences.
        intron_list (List[Dict]): A list of dictionaries representing intron information.
        exon_list (List[Dict]): A list of dictionaries representing exon information.
    """
    status: int
    cds_list: List[str]
    intron_list: List[dict]
    exon_list: List[dict]


def run_genscan(sequence: str = None, sequence_file: str = None, organism: str = "Vertebrate",
                exon_cutoff: float = 1.00, sequence_name: str = ""):
    """
    API for the service http://hollywood.mit.edu/GENSCAN.html.
    :param sequence: The DNA sequence to be analyzed.
    :param sequence_file: The path to the sequence file to be analyzed.
    :param organism: The organism type for Genscan analysis. "Vertebrate", "Arabidopsis" and "Maize" are available. Defaults to "Vertebrate".
    :param exon_cutoff: The exon cutoff value. Defaults to 1.00.
    :param sequence_name: The name of the sequence (optional)
    :return: Dataclass object containing the analysis results including status, coding sequences, intron information, and exon information.
    """
    url = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi"
    if sequence_file:
        file = open(sequence_file, 'rb')
        form_data = {
            "-o": organism,
            "-e": exon_cutoff,
            "-n": sequence_name,
            "-u": file.read(),
            "-p": "Predicted peptides only"
        }
        file.close()
    else:
        form_data = {
            "-o": organism,
            "-e": exon_cutoff,
            "-n": sequence_name,
            "-s": sequence,
            "-p": "Predicted peptides only"
        }
    response = requests.post(url, data=form_data)
    soup = BeautifulSoup(response.content, 'html.parser')
    status = response.status_code

    element = soup.find_all('pre')[0]
    element_list = element.text.split('\n')
    filtered_list = list(filter(lambda x: x.strip(), element_list))
    exon_list = []
    current_sequence = []
    cds_list = []
    for row in filtered_list:
        if row.startswith(' '):
            row_list = row.split()
            exon_number = int(row_list[0].split('.')[1])
            exon_type = row_list[1]
            exon_start = int(row_list[3])
            exon_end = int(row_list[4])
            exon_dict = {exon_number: {'Type': exon_type,
                                       'Begin': exon_start,
                                       'End': exon_end}}
            exon_list.append(exon_dict)
        if row.startswith('>'):
            if current_sequence:
                cds_list.append(''.join(current_sequence))
                current_sequence = []
        else:
            if row.isupper():
                current_sequence.append(row)
    cds_list.append(''.join(current_sequence))

    excluded_list = ['Prom', 'Term', 'PlyA']
    filtered_exon_list = [item for item in exon_list if list(item.values())[0]['Type'] not in excluded_list]
    intron_list = []
    for i in range(1, len(filtered_exon_list)):
        previous_key = list(filtered_exon_list[i - 1].keys())[0]
        next_key = list(filtered_exon_list[i].keys())[0]
        intron_dict = {
            i: {'Begin': filtered_exon_list[i - 1][previous_key]['End'] + 1, 'End': filtered_exon_list[i][next_key]['Begin'] - 1}}
        intron_list.append(intron_dict)

    return GenscanOutput(status, cds_list, intron_list, exon_list)
