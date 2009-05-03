#!/usr/bin/env python
"""\
 TessSweet: Runs the entire PyQuante Test Suite. Has the
 results from runs on 2002-09-08.

 The user has the option of setting do_xx flags, which also run the
 test jobs of type xx, or the do_only_xx jobs, which *only* run
 those jobs of type xx. These flags are passed in as keyword arguments. 

 The test suite is called TessSweet because my wife, Tess, suggested
 it. Woe be unto the man who denies his pregnant wife.
"""

import logging
import h2,he,he_dft,h2o,h2o_mindo,oh_mindo,h2o_dft,ne,no_uhf,h2_cis,\
       h2_mp2,h2_dft,lih_dft,h_dft,li_dft,no_dft,\
       h2_ft_dft,li_ft_dft,be_oep
import time
import datetime

def main(**opts):
    echo = opts.get('echo',False)
    logging.basicConfig(#filename='testsuite.log',
                        level=logging.INFO,
                        format="%(message)s",
                        filemode='w')
    logging.info("Starting Python Test Suite")
    logging.info(time.asctime())

    t1 = time.time()
    nfailed = 0

    do_hf = opts.get('do_hf',True)
    do_only_hf = opts.get('do_only_hf',False)
    do_dft = opts.get('do_dft',True)
    do_only_dft = opts.get('do_only_dft',False)
    do_mindo = opts.get('do_mindo',True)
    do_only_mindo = opts.get('do_only_mindo',False)
    do_other = opts.get('do_other',True)
    do_only_other = opts.get('do_only_other',False)
    do_uhf = opts.get('do_uhf',False) # F b/c it's slow
    do_only_uhf = opts.get('do_only_uhf',False)

    hf_jobs = [h2,he,h2o,ne]
    dft_jobs = [h_dft,h2_dft,h2_ft_dft,he_dft,li_dft,li_ft_dft,lih_dft,
                h2o_dft,no_dft,be_oep]
    mindo_jobs = [h2o_mindo,oh_mindo]
    other_jobs = [h2_cis,h2_mp2]
    uhf_jobs = [no_uhf]

    joblist = []
    if do_only_hf:
        joblist = hf_jobs
    elif do_only_dft:
        joblist = dft_jobs
    elif do_only_mindo:
        joblist = mindo_jobs
    elif do_only_other:
        joblist = other_jobs
    elif do_only_uhf:
        joblist = uhf_jobs
    else:
        if do_hf: joblist.extend(hf_jobs)
        if do_dft: joblist.extend(dft_jobs)
        if do_mindo: joblist.extend(mindo_jobs)
        if do_other: joblist.extend(other_jobs)
        if do_uhf: joblist.extend(uhf_jobs)

    for job in joblist:
        logging.info("Running %-15s:" % job.name)
        en = job.main()
        error = abs(en-job.energy)
        if error < 1.e-4:
            logging.info("--- E=%12.6f Worked ---" % en)
        else:
            logging.info("*** Warning: E=%12.6f should be %12.6f ***" %
                      (en,job.energy))
            nfailed += 1
    t2 = time.time()
    logging.info("Total time for test suite = %f" % (t2-t1))
    logging.info("%d out of %d test cases failed" % (nfailed,len(joblist)))
    return

def profmain():
    import cProfile,pstats
    cProfile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__': main()
    
