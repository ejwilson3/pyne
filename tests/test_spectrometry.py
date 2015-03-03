"""Spectrometry tests """
import nose 
from nose.tools import assert_equal, assert_true

import numpy as np

import warnings


from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne import gammaspec

gspec1 = gammaspec.read_spe_file('test.spe')
gspec2 = gammaspec.read_dollar_spe_file("gv_format_spect.spe")
eff_coeff = [-2.818615042612040000, -0.727352820018942000, -0.039579888648190400,
             -0.059230525466409600, 0.023772637347443000, 0.032530647507267100]


def test_read_dollar_spe():
    assert_equal(gspec2.real_time, 209)
    assert_equal(gspec2.live_time, 199)
    assert_equal(gspec2.start_time, "11:43:41  ")
    assert_equal(gspec2.start_date, "08/01/2014")
    assert_equal(gspec2.dead_time, 9.300003)
    assert_equal(gspec2.det_id, "2")
    assert_equal(gspec2.det_descp, "DSPEC1")
    assert_equal(gspec2.start_chan_num, 0)
    assert_equal(gspec2.num_channels, 1024)
    
def test_read_spe():
    assert_equal(gspec1.real_time, 209.100006)
    assert_equal(gspec1.live_time, 199.800003)
    assert_equal(gspec1.start_time, "11:43:41  ")
    assert_equal(gspec1.start_date, "01-Aug-2014")
    assert_equal(gspec1.dead_time, 9.300003)
    assert_equal(gspec1.det_id, "2")
    assert_equal(gspec1.det_descp, "DSPEC1")
    assert_equal(gspec1.start_chan_num, 0)
    assert_equal(gspec1.num_channels, 1024)
    
def test_calib():
    assert_equal(gammaspec.calc_e_eff(1, eff_coeff, 1), 0.059688551591347033)

def test_str():
    s = str(gspec1)
    assert_true(len(s) > 0)

if __name__ == "__main__":
    nose.runmodule()
