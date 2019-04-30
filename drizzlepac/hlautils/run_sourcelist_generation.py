""" script that directly calls sourcelist_generation.create_sourcelists() using an obs_info_dict hardcoded to
acs_10265_01 values. This will speed development of sourcelist_generation.py becasue it allows the user to just run it
without having to run runhlaprocessing.py first.

"""
import sourcelist_generation
import collections


collections.OrderedDict()
obs_info_dict = collections.OrderedDict()

obs_info_dict = \
    [('filter product 00',
      {'info': '10265 01S ACS WFC F606W',
       'files': ['j92c01b4q_flc.fits', 'j92c01b5q_flc.fits', 'j92c01b7q_flc.fits', 'j92c01b9q_flc.fits'],
       'product filenames': {'image': 'hst_10265_01S_ACS_WFC_F606W.fits',
                             'source catalog': 'hst_10265_01S_ACS_WFC_F606W.cat'},
       'subproduct #0 filenames': {'image': 'hst_10265_01S_ACS_WFC_F606W_j92c01b4q_00.fits',
                                   'source catalog': 'hst_10265_01S_ACS_WFC_F606W_j92c01b4q_00.cat'},
       'subproduct #1 filenames': {'image': 'hst_10265_01S_ACS_WFC_F606W_j92c01b5q_01.fits',
                                   'source catalog': 'hst_10265_01S_ACS_WFC_F606W_j92c01b5q_01.cat'},
       'subproduct #2 filenames': {'image': 'hst_10265_01S_ACS_WFC_F606W_j92c01b7q_02.fits',
                                   'source catalog': 'hst_10265_01S_ACS_WFC_F606W_j92c01b7q_02.cat'},
       'subproduct #3 filenames': {'image': 'hst_10265_01S_ACS_WFC_F606W_j92c01b9q_03.fits',
                                   'source catalog': 'hst_10265_01S_ACS_WFC_F606W_j92c01b9q_03.cat'}}),
     ('total detection product 00',
      {'info': '10265 01S ACS WFC',
       'files': ['j92c01b4q_flc.fits', 'j92c01b5q_flc.fits', 'j92c01b7q_flc.fits', 'j92c01b9q_flc.fits'],
       'product filenames': {'image': 'hst_10265_01S_ACS_WFC.fits',
                             'source catalog': 'hst_10265_01S_ACS_WFC.cat'}})]
sourcelist_generation.create_sourcelists(obs_info_dict)
