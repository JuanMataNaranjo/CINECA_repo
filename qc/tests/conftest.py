from pytest import fixture
from ..fastqc import Fastqc
from ..qc.bwa import Bwa


@fixture
def read_bwa():
    return Bwa(path='data/bwa/', sample='SRR062634', table_path='data/fastq.csv')


@fixture
def read_fastqc():
    return Fastqc('data/fastqc/', sample='SRR062634')
