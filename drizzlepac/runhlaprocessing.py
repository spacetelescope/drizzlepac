#!/usr/bin/env python

"""This script duplicates the functionality of the HLA pipeline.

"""
from drizzlepac import generate_final_product_filenames
from drizzlepac import runastrodriz

def run_processing():
    """
    Main calling subroutine.

    Parameters
    ----------

    Returns
    --------
    """

    # 1: Interpret input csv file as an astropy table with defined column names (HLA-211)

    # 2: Apply rules to determine what exposures need to be combined into separate products (HLA-211 or a new ticket if necessary)
    # PLACEHOLDER UNTIL STEP 2 IS COMPLETED
    obs_info_dict = {}
    obs_info_dict["single exposure product 00"] = "50 A1S WFC3 IR F110W ia1s70jrq" #test proposal_id padding
    obs_info_dict["single exposure product 01"] = "11150 A1S WFC3 UVIS F110W ia1s70jtq"
    obs_info_dict["single exposure product 02"] = "11150 A1S WFC3 IR F110W ia1s70jvq"
    obs_info_dict["single exposure product 03"] = "11150 A1S WFC3 IR F110W ia1s70jwq"
    obs_info_dict["single exposure product 04"] = "11150 A1S WFC3 IR F160W ia1s70jkq"
    obs_info_dict["single exposure product 05"] = "11150 A1S WFC3 IR F160W ia1s70jmq"
    obs_info_dict["single exposure product 06"] = "11150 A1S WFC3 IR F160W ia1s70joq"
    obs_info_dict["single exposure product 07"] = "11150 A1S WFC3 IR F160W ia1s70jpq"
    obs_info_dict["single exposure product 08"] = "10182 A1S ACS HRC PR200LPOL120UV j90za1hyq" #determine maximum generated name length
    obs_info_dict["filter product 00"] = "11150 A1S WFC3 IR F110W"
    obs_info_dict["filter product 01"] = "11150 A1S WFC3 IR F160W"
    obs_info_dict["total detection product 00"] = "11150 A1S WFC3 IR"
    obs_info_dict['multivisit mosaic product 00'] = "1234567 ACS WFC F606W"

    # 3: For each defined product...
    ctr=1
    for obs_category in obs_info_dict.keys():
    #   3.1: generate an output name
        product_filename_dict = generate_final_product_filenames.run_generator(obs_category, obs_info_dict[obs_category])
        for key in product_filename_dict.keys():
            print("----> ", len(product_filename_dict[key]), key,
                  product_filename_dict[key])
        input()
        ctr += 1
    #   3.2: Run astrodrizzle on inputs which define the new product using parameters defined by HLA along with the newly defined output name

    #   3.3: Create source catalog from newly defined product (HLA-204)

    #   3.4: (OPTIONAL) Determine whether there are any problems with alignment or photometry of product

    # 4: (OPTIONAL/TBD) Create trailer file for new product to provide information on processing done to generate the new product.

    # 5: Return exit code for use by calling Condor/OWL workflow code: 0 (zero) for success, 1 for error condition

if __name__ == '__main__':
    run_processing()