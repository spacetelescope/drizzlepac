#!/usr/bin/env python

"""This script is development space for code that will be most likely abosorbed into a larger script later on.

"""
import pdb
def run_generator(product_category,obs_info):
    """
    This is the main calling subroutine. It decides which filename generation subroutine should be run based on the
    input product_category, and then passes the information stored in input obs_info to the subroutine so that the
    appropriate filenames can be generated.

    Parameters
    ----------
    product_category : string
        The type of final output product which filenames will be generated for
    :param obs_info: string
        A string containing space-separated items that will be used to generate the filenames.

    :return: dictionary
        A dictionary containing the generated filenames.
    """
    category_generator_mapping = {'single exposure product': single_exposure_product_filename_generator,
                                  'filter product': filter_product_filename_generator,
                                  'total detection product': total_detection_product_filename_generator,
                                  'multivisit mosaic product': multivisit_mosaic_product_filename_generator}

    for key in category_generator_mapping.keys():
        if product_category.startswith(key):
            generator_name = category_generator_mapping[key]
            category_num = product_category.replace(key+" ","")
            break

    obs_info = obs_info.split(" ")
    if key != "multivisit mosaic product":
        obs_info[0] = "{}{}".format("0"*(5-len(obs_info[0])),obs_info[0])

    product_filename_dict=generator_name(obs_info,category_num)

    return(product_filename_dict)
# ----------------------------------------------------------------------------------------------------------------------

def single_exposure_product_filename_generator(obs_info,nn):
    proposal_id = obs_info[0]
    visit_id = obs_info[1]
    instrument = obs_info[2]
    detector = obs_info[3]
    filter = obs_info[4]
    ipppssoot = obs_info[5]

    product_filename_dict = {}
    product_filename_dict["image"] = "hst_{}_{}_{}_{}_{}_{}_{}.fits".format(proposal_id,visit_id,instrument,detector,filter,ipppssoot,nn)
    product_filename_dict["source catalog"]= product_filename_dict["image"].replace(".fits",".cat")

    return(product_filename_dict)

# ----------------------------------------------------------------------------------------------------------------------

def filter_product_filename_generator(obs_info,nn):
    proposal_id = obs_info[0]
    visit_id = obs_info[1]
    instrument = obs_info[2]
    detector = obs_info[3]
    filter = obs_info[4]

    product_filename_dict = {}
    product_filename_dict["image"] = "hst_{}_{}_{}_{}_{}.fits".format(proposal_id,visit_id,instrument,detector,filter)
    product_filename_dict["source catalog"] = product_filename_dict["image"].replace(".fits",".cat")

    return(product_filename_dict)


# ----------------------------------------------------------------------------------------------------------------------

def total_detection_product_filename_generator(obs_info,nn):
    proposal_id = obs_info[0]
    visit_id = obs_info[1]
    instrument = obs_info[2]
    detector = obs_info[3]

    product_filename_dict = {}
    product_filename_dict["image"] = "hst_{}_{}_{}_{}.fits".format(proposal_id, visit_id, instrument, detector)
    product_filename_dict["source catalog"] = product_filename_dict["image"].replace(".fits",".cat")

    return (product_filename_dict)

# ----------------------------------------------------------------------------------------------------------------------

def multivisit_mosaic_product_filename_generator(obs_info,nn):
    group_num = obs_info[0]
    instrument = obs_info[1]
    detector = obs_info[2]
    filter = obs_info[3]

    product_filename_dict = {}
    product_filename_dict["image"] = "hst_mos_{}_{}_{}_{}.fits".format(group_num,instrument,detector,filter)
    product_filename_dict["source catalog"] = product_filename_dict["image"].replace(".fits",".cat")

    return (product_filename_dict)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    #FOR TESTING!
    obs_info_dict = {}
    obs_info_dict["single exposure product 00"] = "1150 A1S WFC3 IR F110W ia1s70jrq"
    obs_info_dict["single exposure product 01"] = "11150 A1S WFC3 IR F110W ia1s70jtq"
    obs_info_dict["single exposure product 02"] = "11150 A1S WFC3 IR F110W ia1s70jvq"
    obs_info_dict["single exposure product 03"] = "11150 A1S WFC3 IR F110W ia1s70jwq"
    obs_info_dict["single exposure product 04"] = "11150 A1S WFC3 IR F160W ia1s70jkq"
    obs_info_dict["single exposure product 05"] = "11150 A1S WFC3 IR F160W ia1s70jmq"
    obs_info_dict["single exposure product 06"] = "11150 A1S WFC3 IR F160W ia1s70joq"
    obs_info_dict["single exposure product 07"] = "11150 A1S WFC3 IR F160W ia1s70jpq"
    obs_info_dict["filter product 00"] = "11150 A1S WFC3 IR F110W"
    obs_info_dict["filter product 01"] = "11150 A1S WFC3 IR F160W"
    obs_info_dict["total detection product 00"] = "11150 A1S WFC3 IR"
    obs_info_dict['multivisit mosaic product 00'] = "1234567 ACS WFC F606W"

    ctr=1
    for item in obs_info_dict.keys():
        print("{}: {} {}".format(ctr,item,obs_info_dict[item]))
        product_filename_dict = run_generator(item,obs_info_dict[item])

        for key in product_filename_dict.keys():
            print("----> ", len(product_filename_dict[key]), key,
                  product_filename_dict[key])
        input()

        ctr+=1
