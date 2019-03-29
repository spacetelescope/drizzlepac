#!/usr/bin/env python

from astropy.table import Table, Column

#Define information/formatted strings to be included in output dict
sep_str = 'single exposure product {:02d}'
fp_str = 'filter product {:02d}'
tdp_str = 'total detection product {:02d}'

"""This script is development space for code that will be most likely
    absorbed into a larger script later on.

"""
def interpret_input(results, exposure_limit=2.0):
    """
    Interpret the database query for a given visit to prepare the returned
    values for use in generating the names of all the expected output products.
    
    Input will have format of:
        ia1s70jtq_flt.fits,11150,A1S,70,70,149.232269,F110W,IR,/ifs/archive/ops/hst/public/ia1s/ia1s70jtq/ia1s70jtq_flt.fits
        ia1s70iwq_flt.fits,11150,A1S,70,70,0.91161000000000003,F160W,IR,/ifs/archive/ops/hst/public/ia1s/ia1s70iwq/ia1s70iwq_flt.fits
        
    Output dict will have format (as needed by further code for creating the
        product filenames) of:
    
        obs_info_dict["single exposure product 00"] = {'info': '11150 A1S WFC3 IR F110W IA1S70JTQ', 
                                                        'files':['ia1s70jtq_flt.fits']}
        obs_info_dict["single exposure product 01"] = {'info': '11150 A1S WFC3 IR F160W IA1S70JWQ', 
                                                        'files':['ia1s70jwq_flt.fits']}
                                                        
        obs_info_dict["filter product 00"] = {"info": '11150 A1S WFC3 IR F110W', 
                                              "files":['ia1s70jtq_flt.fits']}
        obs_info_dict["filter product 01"] = {"info": '11150 A1S WFC3 IR F160W', 
                                              "files":['ia1s70jwq_flt.fits']}
        obs_info_dict["total detection product 00"] = {'info': '11150 A1S WFC3 IR',
                                                       'files':['ia1s70jtq_flt.fits',
                                                                'ia1s70iwq_flt.fits']}
    
    """
    colnames = ['filename', 'proposal_id', 'program_id', 'obset_id', 'visit_id',
                'exptime', 'filters', 'detector', 'pathname']
    instrument_dict = {'i':'WFC3', 'j':'ACS','o':'STIS','u':'WFPC2','x':'FOC','w':'WFPC'}
    visit_table = Table.read(results, format='ascii.fast_no_header')
    # label each column with a descriptive name
    for col, cname in zip(visit_table.colnames, colnames):
        visit_table[col].name = cname
    # Add INSTRUMENT column
    instr= instrument_dict[visit_table['filename'][0][0]]
    # convert input to an Astropy Table for parsing
    visit_table.add_column(Column([instr]*len(visit_table)),name='instrument')
    # parse Table into a tree-like dict
    visit_tree = build_visit_tree(visit_table, exposure_limit=exposure_limit)
    # Now create final dict
    visit_dict = parse_visit_tree(visit_tree)
    
    return visit_dict

# Translate the database query on a visit into actionable lists of filenames
def build_visit_tree(visit_table, exposure_limit=2.0):
    """Convert visit table into a tree listing all products to be created."""

    # Each product will consist of the appropriate string as the key and
    # a dict of 'info' and 'files' information
    base_entry = {'info':"", 'files':[]}

    # Start interpreting the visit table
    visit_tree = {}
    for row in visit_table:
        # Check for multiple instruments
        instr = row['instrument']
        det = row['detector']
        filt = row['filters']
        expt = row['exptime']
        row_info, filename = create_row_info(row)
        if instr not in visit_tree:
            visit_tree[instr] = {}
            visit_tree[instr][det] = {}
            visit_tree[instr][det][filt] = {}
            visit_tree[instr][det][filt][expt] = [(row_info, filename)]
        else:
            instr_node = visit_tree[instr]
            if det not in instr_node:
                instr_node[det] = {}
                instr_node[det][filt] = {}
                instr_node[det][filt][expt] = [(row_info, filename)]

            else:
                det_node = instr_node[det]
                if filt not in det_node:
                    det_node[filt] = {}
                    det_node[filt][expt] = [(row_info, filename)]
                else:
                    filt_node = det_node[filt]
                    add_expt = False
                    for e in filt_node:
                        if 1./exposure_limit <= expt/e <= exposure_limit:
                            filt_node[e].append((row_info, filename))
                            add_expt = True
                            break
                    if not add_expt:
                        filt_node[expt] = [(row_info, filename)]
    return visit_tree

def create_row_info(row):
    """Build info string for a row from the visit table"""
    instr = row['instrument']
    propid = row['proposal_id']
    progid = row['program_id']
    det = row['detector']
    filt = row['filters']
    rootname = row['filename'][:row['filename'].find('_')]
    info_list = [propid, progid, instr, det, filt, rootname]
    row_info = ''
    for info in info_list:
        row_info += '{} '.format(info)
    return row_info.strip().upper(), row['filename']

def parse_visit_tree(tree):
    """Convert tree into products 
    
    Tree generated by `create_row_info()` will always have the following 
    levels:
      * Instrument
          * Detector 
              * filters
                  * exptime
                      * exposure
    Each exposure will have an entry dict with keys 'info' and 'filename'.
    
    Products created will be:
      * total detection product
      * filter products (one for each "exposure level")
      * single exposure product
    """
    # Initialize products dict
    visit_products = {}

    # for each level, define a product, starting with the instrument used...
    sep_indx = 0
    for inum, instr in enumerate(tree.keys()):
        # intialize now, populate after all inputs are known
        tdp = tdp_str.format(inum)
        visit_products[tdp] = {'info':"", 'files':[]}
        det_tree = tree[instr]
        # setup products for each detector used
        for detnum, det in enumerate(det_tree.keys()):
            filt_tree = det_tree[det]
            # find all filters user...
            for filtnum, filters in enumerate(filt_tree.keys()):
                # each filter may have been taken with multiple exptimes...
                exp_tree = filt_tree[filters]
                for expnum, exptime in enumerate(exp_tree.keys()):
                    #create consecutive filt indices that spans all exptimes 
                    filt_indx = filtnum + expnum 
                    # use this to create and populate filter products entry
                    fp = fp_str.format(filt_indx)
                    visit_products[fp] = {'info':"", 'files':[]}
                    # populate single exposure entry now as well
                    for filename in exp_tree[exptime]:
                        sep = sep_str.format(sep_indx) # keep 80 char wide code
                        visit_products[sep] = {'info':filename[0], 
                                              'files':filename[1]}
                                              
                        # Initialize `info` key for this filter product
                        if not visit_products[fp]['info']:
                            # use all but last entry in exposure level info
                            fp_info = " ".join(filename[0].split()[:-1])
                            visit_products[fp]['info'] = fp_info
                        # populate filter product with input filename
                        visit_products[fp]['files'].append(filename[1])
                        # Initialize `info` key for total detection product
                        if not visit_products[tdp]['info']:
                            tdp_info = " ".join(filename[0].split()[:-2])
                            visit_products[tdp]['info'] = tdp_info
                        # append exposure filename to input list for total detection product
                        visit_products[tdp]['files'].append(filename[1])
                        # Increment single exposure master index
                        sep_indx += 1
                        
    # Done... return dict
    return visit_products
                            
                                                
        


def run_generator(product_category,obs_info):
    """
    This is the main calling subroutine. It decides which filename generation subroutine should be run based on the
    input product_category, and then passes the information stored in input obs_info to the subroutine so that the
    appropriate filenames can be generated.

    Parameters
    ----------
    product_category : string
        The type of final output product which filenames will be generated for
    obs_info : string
        A string containing space-separated items that will be used to
        generate the filenames.

    Returns
    --------
    product_filename_dict : dictionary
        A dictionary containing the generated filenames.
    """
    category_generator_mapping = {'single exposure product': single_exposure_product_filename_generator,
                                  'filter product': filter_product_filename_generator,
                                  'total detection product': total_detection_product_filename_generator,
                                  'multivisit mosaic product': multivisit_mosaic_product_filename_generator}

    # Determine which name generator to use based on input product_category
    for key in category_generator_mapping.keys():
        if product_category.startswith(key):
            generator_name = category_generator_mapping[key]
            category_num = product_category.replace(key+" ","")
            break

    # parse out obs_info into a list
    obs_info = obs_info.split(" ")

    # pad 4-character proposal_id values with leading 0s so that proposal_id is
    # a 5-character sting.
    if key != "multivisit mosaic product": # pad
        obs_info[0] = "{}{}".format("0"*(5-len(obs_info[0])),obs_info[0])

    # generate and return filenames
    product_filename_dict=generator_name(obs_info,category_num)
    return(product_filename_dict)
# ----------------------------------------------------------------------------------------------------------------------

def single_exposure_product_filename_generator(obs_info,nn):
    """
    Generate image and sourcelist filenames for single-exposure products

    Parameters
    ----------
    obs_info : list
        list of items that will be used to generate the filenames: proposal_id,
        visit_id, instrument, detector, filter, and ipppssoot
    nn : string
        the single-exposure image number (NOTE: only used in
        single_exposure_product_filename_generator())

    Returns
    --------
    product_filename_dict : dictionary
        A dictionary containing the generated filenames.
    """
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
    """
    Generate image and sourcelist filenames for filter products

    Parameters
    ----------
    obs_info : list
        list of items that will be used to generate the filenames: proposal_id,
        visit_id, instrument, detector, and filter
    nn : string
        the single-exposure image number (NOTE: only used in
        single_exposure_product_filename_generator())

    Returns
    --------
    product_filename_dict : dictionary
        A dictionary containing the generated filenames.
    """
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
    """
    Generate image and sourcelist filenames for total detection products

    Parameters
    ----------
    obs_info : list
        list of items that will be used to generate the filenames: proposal_id,
        visit_id, instrument, and detector
    nn : string
        the single-exposure image number (NOTE: only used in
        single_exposure_product_filename_generator())

    Returns
    --------
    product_filename_dict : dictionary
        A dictionary containing the generated filenames.
    """
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
    """
    Generate image and sourcelist filenames for multi-visit mosaic products

    Parameters
    ----------
    obs_info : list
        list of items that will be used to generate the filenames: group_id,
        instrument, detector, and filter
    nn : string
        the single-exposure image number (NOTE: only used in
        single_exposure_product_filename_generator())

    Returns
    --------
    product_filename_dict : dictionary
        A dictionary containing the generated filenames.
    """
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

    ctr=1
    for item in obs_info_dict.keys():
        print("{}: {} {}".format(ctr,item,obs_info_dict[item]))
        product_filename_dict = run_generator(item,obs_info_dict[item])

        for key in product_filename_dict.keys():
            print("----> ", len(product_filename_dict[key]), key,
                  product_filename_dict[key])
        input()
        ctr+=1
