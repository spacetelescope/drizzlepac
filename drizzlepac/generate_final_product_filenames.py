#!/usr/bin/env python

"""This script is development space for code that will be most likely abosorbed into a larger script later on.

"""
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
            A stirng containing space-seperated items that will be used to generate the filenames.

        :return: dictionary
            A dictionary containing the generated filenames.
        """
        category_generator_mapping = {'single exposure product': single_exposure_product_filename_generator,
                                      'filter product': filter_product_filename_generator,
                                      'total detection product': total_detection_product_filename_generator,
                                      'multivisit mosaic product': multivisit_mosaic_product_filename_generator}
        for prod_cat in category_generator_mapping:
            if prod_cat.replace("")
            generator_name

# ----------------------------------------------------------------------------------------------------------------------

    def single_exposure_product_filename_generator(obs_info):
        print("single_exposure_product_filename_generator")

# ----------------------------------------------------------------------------------------------------------------------

    def filter_product_filename_generator(obs_info):
        print("filter_product_filename_generator")

# ----------------------------------------------------------------------------------------------------------------------

    def total_detection_product_filename_generator(obs_info):
        print("total_detection_product_filename_generator")

# ----------------------------------------------------------------------------------------------------------------------

    def multivisit_mosaic_product_filename_generator(obs_info):
        print("multivisit_mosaic_product_filename_generator")

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    obs_info_dict = {}
    obs_info_dict["single exposure product 00"]="11150 A1S WFC3 IR F110W ia1s70jrq"
    obs_info_dict["single exposure product 01"]="11150 A1S WFC3 IR F110W ia1s70jtq"
    obs_info_dict["single exposure product 02"]="11150 A1S WFC3 IR F110W ia1s70jvq"
    obs_info_dict["single exposure product 03"]="11150 A1S WFC3 IR F110W ia1s70jwq"
    obs_info_dict["single exposure product 04"]="11150 A1S WFC3 IR F160W ia1s70jkq"
    obs_info_dict["single exposure product 05"]="11150 A1S WFC3 IR F160W ia1s70jmq"
    obs_info_dict["single exposure product 06"]="11150 A1S WFC3 IR F160W ia1s70joq"
    obs_info_dict["single exposure product 07"]="11150 A1S WFC3 IR F160W ia1s70jpq"
    obs_info_dict["filter product 00"]="11150 A1S WFC3 IR F110W"
    obs_info_dict["filter product 01"]="11150 A1S WFC3 IR F160W"
    obs_info_dict["total detection product 00"]="11150 A1S WFC3 IR"
