#from enum import unique
import urllib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy import coordinates as coord
from astropy import units as u
from astropy.table import Table, join, unique, Column
from astropy.io.ascii import convert_numpy
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import pandas as pd
from zero_point import zpt
import math

def get_nea():

    #Get exoplanet data from NASA exoplanet archive, specifiying different datasets for those with gaia ids and those without

    global NEA_ID
    NEA = NasaExoplanetArchive.query_criteria(table="pscomppars", select="*")
    NEA_ID = NEA[NEA['gaia_id']!='']
    NEA_NO_ID = NEA[NEA['gaia_id']=='']
    print(len(NEA_NO_ID), len(NEA_ID),len(NEA_NO_ID) + len(NEA_ID), len(NEA))
    return NEA, NEA_ID, NEA_NO_ID

def DR2_id_to_id(DR2_id):

    #Remove 'Gaia DR2' from the string of the Gaia identifier in the exoplanet archive, making it comparable to those in the Gaia archive
    id = str(DR2_id).replace("Gaia DR2 ", "")
    id = id.strip()
    return id

def col_transform(table, colname):

    #here set table to NEA_ID
    #Transform whole dr2 column to single string

    table[colname]= [DR2_id_to_id(table[colname][i]) for i in range(len(table))]
    print(len(table))
    return table

def select_metallicity(table, met):
    #Select subset of table from nasa exoplanet archive, depending on what type of metallicity is reported
    #This is important as Feh Must be set as Feh in BASTA inputs etc.

    if met == 'FeH':
        table = table[table['st_metratio']=='[Fe/H]']
    elif met == 'MH':
        table = table[table['st_metratio'] == '[M/H]']

    return table

def remove_weird_names(table):
    #THIS NEEDS  MORE WORK - MAKE SURE TYPE IS AN INT AND RAISE VALUE ERROR IF NOT

    # There remove names with gaia in them from the dataset - for now they are the first two transiting planets discovered by Gaia
    table = table[table['gaia_id'] != 'Gaia-1']
    table = table[table['gaia_id'] != 'Gaia-2']
    #if type(table['gaia_id'][i] !=)
    return table


def col_to_int64(table, col):

    #Change a column's datatype to int64
    #I don't know why I wrote this but I think it's becuase there is no in built astropy funct to do it
    #If/when change to pandas this is useless

    table[col]= [(table[col][i]).astype(np.int64) for i in range(len(table))]
    if [type((table[col][i])) == 'numpy.int64' for i in range(len(table))]:
        print('Column value type converted to int64')
    else:
        print('conversion failed')

    return table


def save_table(table):
    #save to NEA_ID - maybe need some more arguments here

    table.write('/Users/afw2/PycharmProjects/hoststarss/dr3/NEA_ID.csv', format = 'csv', overwrite=True)

def get_gaia(user, password, input_path, output_path):

    #Query the gaia archive with your credentials, crossmatch to exoplanet archive and get out and save a table with gaia and exoplanet values
    #Maybe save ADQL queries in txt files for clarity - this is really messy
    #Also maybe chache passowrd here...
    #This is the slowest function

    Gaia.login(user, password)
    gaia_table = Gaia.load_table('gaiaedr3.gaia_source')
    print(f"table = {gaia_table}")
    Gaia.delete_user_table("nea")
    Gaia.upload_table(upload_resource=input_path, table_name="nea", format ='csv')
    job2 = Gaia.launch_job_async(query='SELECT * FROM user_aweeks.nea AS nea, gaiadr3.gaia_source AS gaia WHERE nea.gaia_id = gaia.source_id AND gaia.parallax_over_error >5.0 AND gaia.ruwe < 1.4 ') 
    gaia_phs = job2.get_results()
    gaia_phs.pprint(max_lines=-1) 
    gaia_phs.write(output_path, format = 'csv', overwrite=True)

def join_2mass(user, password, gaia_path, joiner_path, output_path):

    #Join the table to the 2mass catalogue avaiable in gaia archive

    Gaia.login(user, password)
    gaia_table = Gaia.load_table('gaiaedr3.gaia_source')
    print(f"table = {gaia_table}")
    Gaia.delete_user_table("gaia_phs")
    Gaia.upload_table(upload_resource=gaia_path, table_name="gaia_phs", format ='csv')
    gaia_phs = Table.read(gaia_path, format = 'csv')

    tmass_job = Gaia.launch_job_async(query='SELECT gaia.gaia_id, gaia.ra_2, gaia.dec_2, gaia.parallax, ' + \
    'xmatch.original_ext_source_id, xmatch.clean_tmass_psc_xsc_oid, tmass.ra AS ra_tm, tmass.dec AS dec_tm, ' + \
    'err_maj, err_min, tmass.tmass_oid, tmass.designation, xmatch.angular_distance ' + \
    'FROM user_aweeks.nea_dr3 as gaia ' + \
    'JOIN gaiadr3.tmass_psc_xsc_best_neighbour as xmatch USING (source_id) ' + \
    'JOIN gaiadr3.tmass_psc_xsc_join as xjoin ' + \
	    'ON xmatch.original_ext_source_id = xjoin.original_psc_source_id ' + \
    'JOIN gaiadr1.tmass_original_valid AS tmass ' + \
	    'ON xjoin.original_psc_source_id = tmass.designation ' + \
    'ORDER BY gaia.dec ASC')


    nea_dr3_tm_joiner = tmass_job.get_results()
    nea_dr3_tm_joiner.pprint(max_lines=-1) 
    nea_dr3_tm_joiner.write(joiner_path, format = 'csv', overwrite=True)

    Gaia.delete_user_table("nea_dr3_tm_joiner")
    Gaia.upload_table(upload_resource=joiner_path, table_name="nea_dr3_tm_joiner", format ='csv')

    tmass_job2 = Gaia.launch_job_async(query='SELECT * FROM user_aweeks.nea_dr3_tm_joiner AS gntj, ' + \
    'gaiadr1.tmass_original_valid AS tm WHERE gntj.designation = tm.designation')

    nea_gaia_tm = tmass_job2.get_results()
    nea_gaia_tm.pprint(max_lines=-1)
    nea_gaia_tm = astropy.table.join(nea_gaia_tm, gaia_phs, join_type='inner', keys='gaia_id', table_names=['tm_join', 'phs_gaia'])
    nea_gaia_tm = astropy.table.unique(nea_gaia_tm, keys='pl_name')
    nea_gaia_tm.write(output_path, format = 'csv', overwrite=True)
    nea_gaia_tm = Table.read(output_path, format = 'csv')

    return nea_gaia_tm

def join_ap(user, password, input_path, output_path):

    #Join to the astrophysical parameters table from DR3 - essential to get gspspec teff, meh, logg etc, and elemental abundances

    Gaia.login(user, password)
    gaia_table = Gaia.load_table('gaiadr3.astrophysical_parameters')
    print(f"table = {gaia_table}")
    Gaia.delete_user_table("gaia_phs_tm")
    Gaia.upload_table(upload_resource=input_path, table_name="gaia_phs_tm", format ='csv')
    job3 = Gaia.launch_job_async(query='SELECT * FROM user_aweeks.gaia_phs_tm AS nea, gaiadr3.astrophysical_parameters AS gaia WHERE nea.gaia_id = gaia.source_id')
    gaia_phs = job3.get_results()
    gaia_phs.pprint(max_lines=-1)
    gaia_phs.write(output_path, format = 'csv', overwrite=True)

def join_ap_supp(user, password, input_path, output_path):
    #Join to the astrophysical parameters supplementary table from DR3 - essential to get gspspec teff, meh, logg etc, and elemental abundances

    Gaia.login(user, password)
    gaia_table = Gaia.load_table('gaiadr3.astrophysical_parameters')
    print(f"table = {gaia_table}")
    Gaia.delete_user_table("gaia_phs_tm")
    Gaia.upload_table(upload_resource=input_path, table_name="gaia_phs_tm_ap", format ='csv')
    job3 = Gaia.launch_job_async(query='SELECT * FROM user_aweeks.gaia_phs_tm AS nea, gaiadr3.astrophysical_parameters_supp AS gaia WHERE nea.gaia_id = gaia.source_id')
    gaia_phs = job3.get_results()
    gaia_phs.pprint(max_lines=-1)
    gaia_phs.write(output_path, format = 'csv', overwrite=True)

def sweetcat_spectroscopy(input_table):

    #crossmatch to sweetcat catalogue of semi-homogenous exoplaet hsot star parameters

    sweetCat_table_url = "https://sweetcat.iastro.pt/catalog/SWEETCAT_Dataframe.csv"
    converters={'gaia_dr2': [convert_numpy(np.int64)],'gaia_dr3': [convert_numpy(np.int64)] }
    sc = Table.read(sweetCat_table_url, encoding='UTF-8',format='csv', converters=converters)

    table_init = join(input_table, sc, keys_left='source_id', keys_right='gaia_dr2', table_names = [None, '_sc'])
    table = astropy.table.unique(table_init, keys = 'gaia_dr2')
    print(len(table_init), len(table))
    return table_init

def largest_error(lo,hi):

    #select the largest error of a tuple of and lower errors

    err = []
    for i in range(len(lo)):
        if lo[i] > hi[i]:
            error = lo[i]
        else:
            error = hi[i]
        err.append(error)
    #error_arr = np.array([error[i] for i in range(len(lo))])
    error_array = np.array(err)
    return error_array


def error_quadrature(lo, hi):

    #add hi and lo errors in quadrature - this is irrelevant
    err = []
    for i in range(len(lo)):
        error = np.sqrt((lo[i]**2 + hi[i]**2)/2)
    err.append(error)

    error_array = np.array(err)
    return error_array

def Teff_calibration(input, lo, hi, offset=False, systematic=False):

    #Applies global offset to Teff if needed
    #Adds systematic error for teff in quadrature to teh already claculated symmetrized average error

    teff_err_lo = input[lo]
    teff_err_hi = input[hi]
    #apply linear calibration offset
    if offset:
        input['Teff_calibrated'] = input[Teff_calibration()]

    #add and average upper and lower errors
    input['Teff_err_symm'] = np.sqrt((teff_err_lo**2 + teff_err_hi**2)/2)

    if systematic:
        #Add systematic uncertainty in quadrature
        input['Teff_err_final'] = np.sqrt((input['Teff_err_symm']**2 + symm**2)/2)

def calibrate_gspspec_MeH(input):

    #Applies calibration for gspspec Meh (coefficients from table 3, Recio-Blanco 2022)

    coeffs = [0.274, -0.1373, -0.0050, 0.0048]
    correction = np.poly1d((np.flip(coeffs)))(input['logg'])
    #input['MeH_calibrated'] = input['MeH'] + sum([p * input['logg']**(p) for p in coeffs])
    input['MeH_calibrated'] = input['MeH'] + correction

    print('Mean MeH shift: {}'.format(np.mean(correction)))

    return input

def basta_input_table(table,spec_source="gaia", no_tm_crossmatch=False):

    #needs work on adding asymmetric errors
    #could conditionally assign star type here
    #could conditionally assign science cases
    #2mass cuts here
    #photgmag cuts here
    #gspspec quality flag cuts here



    # Positions
    ra = table['ra_2_phs_gaia']
    dec = table['dec_2_phs_gaia']

    #Gaia Magnitude

    gaiamag = table['sy_gaiamag']
    gaiamag_err_lo = table['sy_gaiamagerr1']
    gaiamag_err_hi = table['sy_gaiamagerr2']

    starid = table['gaia_id']

    pllx_gaia = table['parallax_corrected']
    logg = table['st_logg']

    #if no corrected parallax, raise a reminder to correct for it

    pllx_err = table['parallax_error']

    if no_tm_crossmatch:

        # Set values for hjk magnitude from 2mass

        pllx_err = table['parallax_error']

        hm = table['sy_hmag']
        hm_err_hi = table['sy_hmagerr1']
        hm_err_lo = table['sy_hmagerr2']
        hm_err = largest_error(hm_err_hi, hm_err_lo)

        km = table['sy_kmag']
        km_err_hi = table['sy_kmagerr1']
        km_err_lo = table['sy_kmagerr2']
        km_err = largest_error(km_err_hi, km_err_lo)



        jm = table['sy_jmag']
        jm_err_hi = table['sy_jmagerr1']
        jm_err_lo = table['sy_jmagerr2']
        jm_err = largest_error(jm_err_hi, jm_err_lo)

    else:

        hm = table['h_m']
        hm_err = table['h_msigcom']

        km = table['ks_m']
        km_err = table['ks_msigcom']

        jm = table['j_m']
        jm_err = table['j_msigcom']

    #need to add systematic errors to these effective temperatures and feh in order to make them more realistic
    #do this for gaia case only
    #find in gaia literature

    if spec_source == "gaia_phot":

    # Include astrophysical values from gaia's photometry pipeline (general star classifier)

        teff_gaia = table['teff_gspphot']
        teff_gaia_hi = table['teff_gspphot_upper']
        teff_gaia_lo = table['teff_gspphot_lower']
        teff_gaia_err_hi = teff_gaia_hi - teff_gaia
        teff_gaia_err_lo = teff_gaia - teff_gaia_lo

        met_gaia = table['mh_gspspec']
        met_gaia_hi = table['mh_gspspec_upper']
        met_gaia_lo = table['mh_gspspec_lower']
        met_gaia_err_hi = met_gaia_hi - met_gaia
        met_gaia_err_lo = met_gaia - met_gaia_lo

        logg_gaia = table['logg_gspphot']
        logg_gaia_hi = table['logg_gspphot_upper']
        logg_gaia_lo = table['logg_gspphot_lower']
        logg_gaia_err_hi = logg_gaia_hi - logg_gaia
        logg_gaia_err_lo = logg_gaia - logg_gaia_lo

        basta_input_table = Table([starid, ra, dec, teff_gaia, largest_error(teff_gaia_err_hi, teff_gaia_err_lo),
                                   met_gaia, largest_error(met_gaia_err_hi, met_gaia_err_lo), logg,
                                   largest_error(logg_gaia_err_hi, logg_gaia_err_lo), pllx_gaia, pllx_err,
                                   jm, jm_err, km, km_err, hm, hm_err],


        # Make sure all st_met is [Feh/H], this is done in notebook but you should write a function for it!!!! (

                                  names=(
                                      '#starid', 'RA', 'DEC', 'Teff', 'Teff_err', 'MeH', 'MeH_err','logg', 'logg_err',
                                      'parallax',
                                      'parallax_err', 'Mj_2MASS', 'Mj_2MASS_err', 'Mk_2MASS', 'Mk_2MASS_err',
                                      'Mh_2MASS',
                                      'Mh_2MASS_err'), masked=True)

    if spec_source == "gaia_spec_msc":

        # Include spectroscopic values from gaia's spectroscopy pipeline - msc (multiple star classifier)

        met_gaia = table['mh_msc']
        met_gaia_hi = table['mh_msc_upper']
        met_gaia_lo = table['mh_msc_lower']
        met_gaia_err_hi = met_gaia_hi - met_gaia
        met_gaia_err_lo = met_gaia - met_gaia_lo

        teff_gaia = table['teff_msc1']
        teff_gaia_hi = table['teff_msc1_upper']
        teff_gaia_lo = table['teff_msc1_lower']
        teff_gaia_err_hi = teff_gaia_hi - teff_gaia
        teff_gaia_err_lo = teff_gaia - teff_gaia_lo

        logg_gaia = table['logg_msc1']
        logg_gaia_hi = table['logg_msc1_upper']
        logg_gaia_lo = table['logg_msc1_lower']
        logg_gaia_err_hi = logg_gaia_hi - logg_gaia
        logg_gaia_err_lo = logg_gaia - logg_gaia_lo



        basta_input_table = Table([starid, ra, dec, teff_gaia, largest_error(teff_gaia_err_hi, teff_gaia_err_lo),
                                   met_gaia, largest_error(met_gaia_err_hi, met_gaia_err_lo), logg,
                                    largest_error(logg_gaia_err_hi, logg_gaia_err_lo), pllx_gaia, pllx_err,
                                    jm, jm_err, km, km_err, hm, hm_err],

                                    names=(
                                        '#starid', 'RA', 'DEC', 'Teff', 'Teff_err', 'MeH', 'MeH_err', 'logg',
                                        'logg_err',
                                        'parallax',
                                        'parallax_err', 'Mj_2MASS', 'Mj_2MASS_err', 'Mk_2MASS', 'Mk_2MASS_err',
                                        'Mh_2MASS',
                                        'Mh_2MASS_err'), masked=True)

    if spec_source == "gaia_spec_gsp":

        # Include spectroscopic values from gaia's spectroscopy pipeline - gsp (general stellar parameteriser)

        met_gaia = table['mh_gspspec']
        met_gaia_hi = table['mh_gspspec_upper']
        met_gaia_lo = table['mh_gspspec_lower']
        met_gaia_err_hi = met_gaia_hi - met_gaia
        met_gaia_err_lo = met_gaia - met_gaia_lo

        teff_gaia = table['teff_gspspec']
        teff_gaia_hi = table['teff_gspspec_upper']
        teff_gaia_lo = table['teff_gspspec_lower']
        teff_gaia_err_hi = teff_gaia_hi - teff_gaia
        teff_gaia_err_lo = teff_gaia - teff_gaia_lo

        logg_gaia = table['logg_gspspec']
        logg_gaia_hi = table['logg_gspspec_upper']
        logg_gaia_lo = table['logg_gspspec_lower']
        logg_gaia_err_hi = logg_gaia_hi - logg_gaia
        logg_gaia_err_lo = logg_gaia - logg_gaia_lo



        basta_input_table = Table([starid, ra, dec, teff_gaia, largest_error(teff_gaia_err_hi, teff_gaia_err_lo),
                                   met_gaia, largest_error(met_gaia_err_hi, met_gaia_err_lo), logg,
                                    largest_error(logg_gaia_err_hi, logg_gaia_err_lo), pllx_gaia, pllx_err,
                                    jm, jm_err, km, km_err, hm, hm_err],

                                    names=(
                                        '#starid', 'RA', 'DEC', 'Teff', 'Teff_err', 'MeH', 'MeH_err', 'logg',
                                        'logg_err',
                                        'parallax',
                                        'parallax_err', 'Mj_2MASS', 'Mj_2MASS_err', 'Mk_2MASS', 'Mk_2MASS_err',
                                        'Mh_2MASS',
                                        'Mh_2MASS_err'), masked=True)

    elif spec_source == "sc":

        # Include spectroscopic values from sweetcat


        teff_sc = table['Teff']
        teff_sc_err = table['eTeff']
        feh_sc = table['[Fe/H]']
        feh_sc_err = table['e[Fe/H]']

        basta_input_table = Table([starid, ra, dec, teff_sc, teff_sc_err,
                                   feh_sc, feh_sc_err, pllx_gaia, pllx_err,
                                   jm, jm_err, km, km_err, hm, hm_err],

                                  names=(
                                  '#starid', 'RA', 'DEC', 'Teff', 'Teff_err', 'FeH', 'FeH_err',
                                  'parallax',
                                  'parallax_err', 'Mj_2MASS', 'Mj_2MASS_err', 'Mk_2MASS', 'Mk_2MASS_err', 'Mh_2MASS',
                                  'Mh_2MASS_err'), masked=True)

    elif spec_source == "nasa":

        # Include spectroscopic values from Nasa Exoplanet Archive

        feh = table['st_met']
        feh_err_hi = table['st_meterr1']
        feh_err_lo = table['st_meterr2']

        logg = table['st_logg']
        logg_err_hi = table['st_loggerr1']
        logg_err_lo = table['st_loggerr2']

        teff = table['st_teff']
        teff_err_hi = table['st_tefferr1']
        teff_err_lo = table['st_tefferr2']


        basta_input_table = Table([starid, ra, dec, teff, largest_error(teff_err_hi, teff_err_lo),
                                   feh, largest_error(feh_err_hi, feh_err_lo), pllx_gaia, pllx_err,
                                   jm, jm_err, km, km_err, hm, hm_err],

                                  names=(
                                  '#starid', 'RA', 'DEC','Teff', 'Teff_err', 'FeH', 'FeH_err',
                                  'parallax', 'parallax_err', 'Mj_2MASS', 'Mj_2MASS_err', 'Mk_2MASS', 'Mk_2MASS_err', 'Mh_2MASS',
                                  'Mh_2MASS_err'), masked=True)


    print(len(basta_input_table), 'Stars in BASTA table before cleaning')
    return basta_input_table

def remove_bad_stars(table,col_name, print=False):

    #For stars without a value (masked) in a columns, remove from sample

    missing_value_idx = table[col_name].mask.nonzero()[0]
    table.remove_rows(missing_value_idx)

    return table

def clean_table(table):

    #remove bad stars for all columns in table
    #make sure to apply this only once columns are selected, i.e. not if they have missing values for irrelevant columns

    cols = list(table.columns)
    del cols[0]
    for col in cols:
    #col = "'"+col+"'"
        clean_table = remove_bad_stars(table, col)
    print(len(clean_table))
    return clean_table

def remove_nans(table):

    #remove stars with any nans in columns

    has_nan = np.zeros(len(table,), dtype = bool)
    for col in table.itercols():
        if col.info.dtype.kind == 'f':
            has_nan |= np.isnan(col)
    no_nans_table = table[~has_nan]

    return no_nans_table

def parallax_correct(table):

    #Apply lindegren 2021 parallax offset code to correct DR3 parallaxes (around 0.17mas for most but polynomial fitted)


    phot_g_mean_mag = table['phot_g_mean_mag']
    nu_eff_used_in_astrometry = table['nu_eff_used_in_astrometry']
    pseudocolour = table['pseudocolour']
    ecl_lat = table['ecl_lat']
    astrometric_params_solved = table['astrometric_params_solved']

    zpt.load_tables()

    zpt_ = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)
    plx_w_offset = table['parallax_phs_gaia'] - zpt_

    table['parallax_corrected'] = plx_w_offset

    return table 

def remove_2mass_flags(table):

    #Remove stars without 2mass AAA flags in H,J,K photometry

    clean = table[table['ph_qual']=='AAA']
    print(len(table), 'stars tested', len(clean), 'stars with AAA photometry', len(table)-len(clean), 'stars discarded')
    return clean

def save_basta_table(basta_table, filename):
    basta_table.write('/Users/afw2/BASTA/Paper1/{}'.format(filename), overwrite = True,format = 'ascii')

#run basta here lol

def change_col_to_float(table, colname):

    #Change datatype of a column to float, inplace.
    newvar = Column(table[colname], dtype='float')
    return newvar

def mag_sanity(table):

    #Plot 2mass magnitudes in exo archive data and from gaia crossmatch to validate the crossmatch
  
    plt.scatter(table['j_m'], table['sy_jmag'])
    plt.scatter(table['ks_m'], table['sy_kmag'])
    plt.scatter(table['h_m'], table['sy_hmag'])
    plt.show()

def small_planet_sample(table, mass_threshold, radius_threshold):

    #Cut planet sample based on desired planet parameters in sample


    final_idx = (table['rv_flag']==1)*(table['tran_flag']==1)*(table['pl_bmassprov']=='Mass')*(table['pl_rade']<=4)
    intermediate = table[final_idx]

    print(len(intermediate), 'small planets with transit and rv observations found in sample')

    error_cut_idx  = (intermediate['pl_bmasseerr1']/intermediate['pl_bmasse'] <= mass_threshold) * (intermediate['pl_radeerr1']/intermediate['pl_rade'] <= radius_threshold)
    final = intermediate[error_cut_idx]

    final_unique = unique(final, keys='gaia_id')

    print(len(final), 'small planets with transit and rv observations found in sample with at least {} % precision in Mass and {} % in Radius of planet around {} stars'.format(mass_threshold*100, radius_threshold*100, len(final_unique)))

    return final_unique

def diagnostic_plots(input_table, result, path, filename):


    result1 = Table.read(path, format = 'ascii')
    print(len(input_table), len(result1))

    mass = Column(input_table['Mass'], dtype='float') 
    rad = Column(input_table['Rad'], dtype='float')
    pllx = Column(input_table['plx1'], dtype='float') 
    dist = 1000/pllx 
    logage = Column(input_table['logA'], dtype='float') 
    age = 10**logage / 1e6

    plt.rcParams["font.family"] = "serif" 

    fig, axs = plt.subplots(1,4, figsize=(15,2))

    axs[0].scatter(mass,result['massfin'], s=10, color = 'darkslateblue', alpha = 0.2)
    axs[0].set_xlabel(str(r'CKS Mass [$M_{\odot}$]'))
    axs[0].set_ylabel(str(r'BASTA Mass [$M_{\odot}$]'))  

    x = np.linspace(0,max(mass))
    y = x

    axs[0].plot(x,y, linestyle='dashed', color = 'lightgray')
    axs[0].set_aspect('equal', 'box')
    axs[0].margins(x=0, y=0)
    plt.tight_layout()

    axs[1].scatter(rad,result['radPhot'], s=10, color = 'darkslateblue', alpha = 0.2)
    axs[1].set_xlabel(str(r'CKS Radius [$R_{\odot}$]'))
    axs[1].set_ylabel(str(r'BASTA Radius [$R_{\odot}$]'))    

    x = np.linspace(0,max(rad))
    y = x

    axs[1].plot(x,y, linestyle='dashed', color = 'lightgray')
    axs[1].set_aspect('equal', 'box')
    axs[1].margins(x=0, y=0)
    plt.tight_layout()


    pllx = Column(input_table['plx1'], dtype='float') 
    dist = 1000/pllx 

    axs[2].scatter(dist,result['distance'], s=10, color = 'darkslateblue', alpha = 0.2)
    axs[2].set_xlabel(str(r'CKS distance [$pc$]'))
    axs[2].set_ylabel(str(r'BASTA distance [$pc$]'))    
    x = np.linspace(0,max(dist))
    y = x
    axs[2].plot(x,y, linestyle='dashed', color = 'lightgray')
    axs[2].set_aspect('equal', 'box')
    axs[2].margins(x=0, y=0)
    plt.tight_layout()

    logage = Column(input_table['logA'], dtype='float') 
    age = 10**logage / 1e6
    print(age)

    axs[3].scatter(age, result['age'], s=10, color = 'darkslateblue', alpha = 0.3)
    axs[3].set_xlabel(str(r'CKS age [$Myr$]'))
    axs[3].set_ylabel(str(r'Basta Age [$Myr$]'))    
    x = np.linspace(0,result['age'])
    y = x
    axs[3].plot(x,y, linestyle='dashed', color = 'lightgray')
    axs[3].set_aspect('equal', 'box')
    axs[3].margins(x=0, y=0)
    plt.tight_layout()
    plt.savefig('/Users/afw2/PycharmProjects/notebooks/reportfigs/cks_basta_15_02')

    plt.show()

    fig, axs = plt.subplots(1,4, figsize=(15,2))

    plt.rcParams["font.family"] = "serif" 
    axs[0].scatter(mass,mass/result['massfin'], s=10, color = 'darkslateblue', alpha = 0.2)
    axs[0].set_xlabel(str(r'CKS Mass [$M_{\odot}$]'))
    axs[0].set_ylabel(str(r'CKS/BASTA Mass[$M_{\odot}$]'))    
    x = np.linspace(0,max(mass))
    y = np.full(50,1)
    axs[0].plot(x,y, linestyle='dashed', color = 'lightgray')
    #axs[0].set_aspect('equal', 'box')
    axs[0].margins(x=0, y=0)
    plt.tight_layout()

    axs[1].scatter(rad,rad/result['radPhot'], s=10, color = 'darkslateblue', alpha = 0.2)
    axs[1].set_xlabel(str(r'CKS Radius [$R_{\odot}$]'))
    axs[1].set_ylabel(str(r'CKS/BASTA Radius [$R_{\odot}$]'))    
    x = np.linspace(0,max(rad))
    y = y = np.full(50,1)
    axs[1].plot(x,y, linestyle='dashed', color = 'lightgray')
    #axs[1].set_aspect('equal', 'box')
    axs[1].margins(x=0, y=0)
    plt.tight_layout()


    pllx = Column(input_table['plx1'], dtype='float') 
    dist = 1000/pllx 

    axs[2].scatter(dist,dist/result['distance'], s=10, color = 'darkslateblue', alpha = 0.05)
    axs[2].set_xlabel(str(r'CKS distance [$pc$]'))
    axs[2].set_ylabel(str(r'CKS/BASTA distance [$pc$]'))    
    x = np.linspace(0,max(dist))
    y = np.full(50,1)
    axs[2].plot(x,y, linestyle='dashed', color = 'lightgray')
    #axs[2].set_aspect('equal', 'box')
    axs[2].margins(x=0, y=0)
    plt.tight_layout()


    axs[3].scatter(age, age/result['age'], s=10, color = 'darkslateblue', alpha = 0.3)
    axs[3].set_xlabel(str(r'CKS age [$Myr$]'))
    axs[3].set_ylabel(str(r'KCS/BASTA Age [$Myr$]'))    
    x = np.linspace(0,max(age))
    y = np.full(50,1)
    axs[3].plot(x,y, linestyle='dashed', color = 'lightgray')
    #axs[3].set_aspect('equal', 'box')
    axs[3].margins(x=0, y=0)
    plt.tight_layout()
    plt.savefig('/Users/afw2/PycharmProjects/notebooks/reportfigs{}'.format(filename))

    plt.show()


def recalculate_planet_radius(subsample):
    
    #input sample should be table with basta and nasa values of stellar parameters
    # (i.e. results joined with phs input as in previous function)

    #jupiter and solar radii
    rj = 71492
    rs = 696340

    #normalising constant
    r_norm = rs/rj

    incl = subsample['pl_orbincl']
    ecc = subsample['pl_orbeccen']
    p = subsample['pl_orbper']

    rat = subsample['pl_ratror']
    r_1 = (rat*subsample['st_rad'])*r_norm
    r_basta = (rat*subsample['radPhot'])*r_norm

    return r_1, r_basta

def recalculate_planet_mass(subsample):
    m_sun = 1.98847e30
    m_jup = 1.898e27

    mp = subsample['pl_bmassj']
    incl = subsample['pl_orbincl']
    ecc = subsample['pl_orbeccen']
    p = subsample['pl_orbper']

    k = subsample['pl_rvamp']
    rat = subsample['pl_ratror']


    mp_mj = 4.191e-3 * subsample['st_mass'] ** (2 / 3) * p ** (1 / 3) * k * np.sqrt(1 - ecc ** 2) * 1 / (np.sin(incl))
    mp_mj_basta = 4.191e-3 * subsample['massfin'] ** (2 / 3) * p ** (1 / 3) * k * np.sqrt(1 - ecc ** 2) * 1 / (np.sin(incl))

    return mp_mj, mp_mj_basta
def same_gaia_sweetcat_sample(gaia_name, sc_name, save_name):
    path = '/Users/afw2/BASTA/Paper1/'

    gaia = Table.read(path+gaia_name, format = 'ascii')
    sc = Table.read(path+sc_name, format = 'ascii')

    sc_gaia_joined = astropy.table.join(sc, gaia, join_type='inner', keys='starid',
                                     table_names=['tm_join', 'phs_gaia'])
    sc_gaia_joined = astropy.table.unique(sc_gaia_joined, keys='starid')

    sc_gaia_joined.write('{}.ascii'.format(save_name), format = 'ascii', overwrite=True)

def plot_new_pparams(subsample, mp_mj, mp_mj_basta, r_basta, r_const, ecc, incl, r_1, p, k):
    fig, ax = plt.subplots(4, 2, figsize=(15, 25))

    ax[0, 0].scatter(subsample['pl_radj'], r_1, s=1.5, color='k')
    ax[0, 0].set_xlabel('NASA planet radius')
    ax[0, 0].set_ylabel('Recalculated planet radius')

    ax[1, 0].scatter(subsample['pl_radj'], subsample['pl_radj'] - r_1, s=1.5, color='k')
    ax[1, 0].set_xlabel('NASA planet radius')
    ax[1, 0].set_ylabel('Recalculated planet radius residual')
    x = np.linspace(0, max(subsample['pl_radj']))
    y = np.full(50, 0)
    ax[1, 0].plot(x, y)

    ax[0, 1].set_xlim(0.001, 0.2)
    ax[0, 1].set_ylim(0.001, 0.2)

    ax[0, 1].scatter(subsample['pl_bmassj'], mp_mj, s=1.5, color='k')
    ax[0, 1].set_xlabel('NASA planet mass')
    ax[0, 1].set_ylabel('Recalculated planet mass')

    ax[1, 1].scatter(subsample['pl_bmassj'], subsample['pl_bmassj'] - mp_mj, s=1.5, color='k')
    ax[1, 1].set_xlabel('NASA planet mass')
    ax[1, 1].set_ylabel('Recalculated planet mass residual')
    x = np.linspace(0, max(subsample['pl_bmassj']))
    y = np.full(50, 0)
    ax[1, 1].plot(x, y)

    ax[2, 1].set_xlim(0.001, 0.2)
    ax[2, 1].set_ylim(0.001, 0.2)

    ax[2, 1].scatter(subsample['pl_radj'], r_basta, s=1.5, color='k')
    ax[2, 1].set_xlabel('NASA planet radius')
    ax[2, 1].set_ylabel('Recalculated planet radius BASTA')

    ax[3, 1].scatter(subsample['pl_radj'], subsample['pl_radj'] - r_basta, s=1.5, color='k')
    ax[3, 1].set_xlabel('NASA planet radius')
    ax[3, 1].set_ylabel('Recalculated planet radius residual BASTA')
    x = np.linspace(0, max(subsample['pl_radj']))
    y = np.full(50, 0)
    ax[3, 1].plot(x, y)

    ax[2, 0].set_xlim(0.001, 0.2)
    ax[2, 0].set_ylim(0.001, 0.2)

    ax[2, 0].scatter(subsample['pl_bmassj'], mp_mj_basta, s=1.5, color='k')
    ax[2, 0].set_xlabel('NASA planet mass')
    ax[2, 0].set_ylabel('Recalculated planet mass BASTA')

    ax[3, 0].scatter(subsample['pl_bmassj'], subsample['pl_bmassj'] - mp_mj_basta, s=1.5, color='k')
    ax[3, 0].set_xlabel(r'NASA planet mass $R_e$')
    ax[3, 0].set_ylabel('Recalculated planet mass residual BASTA')
    x = np.linspace(0, max(subsample['pl_bmassj']))
    y = np.full(50, 0)
    ax[3, 0].plot(x, y)

    plt.show()

    fig, ax = plt.subplots(figsize=(15, 15))
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1)
    ax.scatter(subsample['pl_bmassj'], subsample['pl_radj'], s=20, c='k', label='NASA values')
    ax.set_xlabel(r'Planet mass [$M_J$]', fontsize=30)
    ax.set_ylabel(r'Planet Radius [$R_J$]', fontsize=30)
    # ax.set_yscale('log')
    ax.set_xscale('log')

    # ax.scatter(mp_mj, r_1, s=9, c='gray', label = 'Recalculated')
    # ax.set_yscale('log')
    # ax.set_xscale('log')

    ax.scatter(mp_mj_basta, r_basta, s=20, c='blue', label='Recalculated BASTA values')
    # ax.patch.set_facecolor('white')
    # ax.patch.set_alpha(1)
    # ax.set_yscale('log')
    ax.set_xscale('log')

    ax.legend(fontsize=30, loc='lower right')
    plt.show()

    plt.hist(subsample['pl_bmassj'] * r_const)
    plt.show()

