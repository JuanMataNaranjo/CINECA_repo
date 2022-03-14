def test_parameters(read_bwa):
    """ basic test to make sure the parameters are read as expected """

    # The sample should be paired
    assert read_bwa.paired == True

    # We should have been able to read the log file
    assert read_bwa.log_file != None

    # The intermediate parameters need to read something
    assert len(read_bwa.mem_pestat) > 0
    assert len(read_bwa.mem_process_seqs) > 0
    assert len(read_bwa.process) > 0

    # Make sure the log files has been read correctly
    assert read_bwa.log_file[0][:-1] == '[M::bwa_idx_load_from_disk] read 3171 ALT contigs'
