import os
from dataclasses import dataclass


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts multiline fasta to oneline
    :param input_fasta: path to your file. If you use Windows, add r before path
    :param output_fasta: name of your file or path to it. Defaults to None.
    :return: None
    """
    if not output_fasta:
        output_filename, ext = os.path.splitext(os.path.basename(input_fasta))
        if not ext or ext != '.fasta':
            ext = '.fasta'
        output_fasta = f'{output_filename}_parsed{ext}'
    line_count = 0
    line_number = 0
    with open(input_fasta) as input_file:
        for _ in input_file.readlines():
            line_count += 1
    with open(input_fasta) as input_file:
        seq = []
        seq_read = False
        with open(output_fasta, mode='w') as output_file:
            for line in input_file.readlines():
                line_number += 1
                if line_number == line_count:
                    if not seq_read:
                        output_file.write(line)
                    else:
                        seq.append(line.strip('\n'))
                        output_file.write(''.join(seq) + '\n')
                elif line.startswith('>'):
                    if not seq_read:
                        output_file.write(line)
                    else:
                        output_file.write(''.join(seq) + '\n' + line)
                        seq = []
                        seq_read = False
                else:
                    seq.append(line.strip('\n'))
                    seq_read = True


@dataclass
class FastaRecord:
    id_: str
    description: str
    seq: str

    def __repr__(self):
        return f'id = {self.id_}\ndescription = {self.description}\nsequence = {self.seq}'


class OpenFasta:
    def __init__(self, file_path: str, mode: str = 'r'):
        self.file_path = file_path
        self.mode = mode
        self.seqs = []
        self.generator = self.generate_records()

    def __enter__(self):
        self.fasta_file = open(self.file_path, mode=self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fasta_file.close()

    def generate_records(self):
        record_id = None
        record_description = ''
        record_seq = ''

        for line in self.fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if record_id:
                    yield FastaRecord(record_id, record_description, record_seq)
                split_line = line[1:].split(maxsplit=1)
                if len(split_line) == 2:
                    record_id, record_description = split_line
                else:
                    record_id = ''.join(split_line)
                    record_description = None
                record_seq = ''
            else:
                record_seq += line
        if record_id:
            yield FastaRecord(record_id, record_description, record_seq)

    def __iter__(self):
        return self.generator

    def read_record(self):
        return next(self.generator)

    def read_records(self):
        return [record for record in self.generator]
