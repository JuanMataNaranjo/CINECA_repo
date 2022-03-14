
def test_parameters(read_fastqc):
    """ basic test to make sure the parameters are read as expected """
    # Make sure these parameters are not False
    assert read_fastqc.log_file_1 != False
    assert read_fastqc.log_file_2 != False

    # Check that the file is long a specific length
    assert len(read_fastqc.log_file_1) == 21
    assert len(read_fastqc.log_file_2) == 21



