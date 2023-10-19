def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = 'input_fasta_name_parsed.fasta') -> None:
    """
    Converts multiline fasta to oneline
    :param input_fasta: str - path to your file. If you use Windows, add r before path
    :param output_fasta: str - name of your file or path to it. Default 'input_fasta_name_parsed.fasta'
    :return: None
    """
    if output_fasta == 'input_fasta_name_parsed.fasta':
        output_filename = input_fasta.split('\\')[-1]
        output_fasta = output_filename.replace('.fasta', '_parsed.fasta')
    line_count = 0
    line_number = 0
    with open(input_fasta) as input_file:
        for line in input_file.readlines():
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
