# Bioinformatics Box

This module contains the learning results of the **Python** course at the **Bioinformatics Institute**.

## Overview:

### `bioinformatics_toolkit.py`

* Classes for working with biological sequences, such as DNA, RNA, and amino acid sequences
* FastQ filter using Biopython
* Python API for the service [GENSCAN](http://hollywood.mit.edu/GENSCAN.html)
* `Telegram_logger` decorator, which logs launches of decorated functions and sends messages to telegram

### `bio_files_processor.py`

* Function for converting multiline fasta into oneline
* Class `FastaRecord` for convenient display of fasta sequences
* Сontext manager for working with fasta files, providing the ability to get each sequence in `FastaRecord` format

### `custom_random_forest.py`

* Custom `RandomForestClassifier` with the ability to determine the number of processes running in parallel

### `test_bioinformatics_box.py`

* Tests that check the module's work

Also, module contains `Showcases.ipynb` in which you can look at examples of module's work.

## Author: I.K. Slepov
