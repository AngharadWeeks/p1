#recalcualte planet masses and radii using new stellar parameters
from astropy.table import Table, join, unique
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from scipy.optimize import root
from astropy import units as u
import pandas as pd
import astropy.constants as c

#subsample = Table.read("/Users/afw2/PycharmProjects/hoststarss/DR3/sc_phs_72.csv", format = 'csv')


def new_r(df, radius_col_name):

    rj = c.R_jup
    rs = c.R_sun
    re = c.R_earth

#normalising constant
    r_norm = rs / re

    incl = df['pl_orbincl']
    ecc = df['pl_orbeccen']
    p = df['pl_orbper']

    rat = df['pl_ratror']
    r_p = df['pl_rade']
    r_p_err = df['pl_radeerr2']

    r_new = (rat* df[radius_col_name])*r_norm

    nun = []

    for i in range(len(r_p)):
        diff = abs(np.array(r_p)[i] - np.array(r_new)[i])
        if diff >= 0.01:
            nun.append([diff,i])
        else:
            pass
    #print((nun))
    #print(len(nun))

    idxs = []
    for an in nun:
        idxs.append(an[1])

    #print(idxs)

    #for i in range(len(idxs)):
    #    idx = idxs[i]
    #    print(idx)
    #    fig, ax = plt.subplots()
    #    ax.errorbar(np.array(r_p)[idx], np.array(r_new)[idx], yerr=np.array(r_p_err)[idx], c='r', marker='^', ms=4, capsize=3, alpha=0.4)
    #    ax.set_xlabel(r'NASA $R_p$')
    #    ax.set_ylabel(r'Recalculated $R_p$')
    #    plt.xlim(0, 0.4)
    #    plt.ylim(0, 0.4)
    #    ax.set_aspect('equal', adjustable='box')
    fig, ax = plt.subplots()
    ax.errorbar(r_p, r_new, yerr=r_p_err , c='k', marker='^', ls='none', ms=4, capsize=3, alpha=0.4)
    ax.set_xlabel(r'NASA R$\rm{_p}$')
    ax.set_ylabel(r'Recalculated R$\rm{_p}$')
    #plt.xlim(0,0.4)
    #plt.ylim(0,0.4)
    ax.set_aspect('equal', adjustable='box')


    plt.show()

    return r_new

def new_rad_err(df, new_r_s_colname):

    rat = df['pl_ratror']

    r_p = df['pl_rade']
    r_p_errp = df['pl_radeerr2']
    r_p_errm = df['pl_radeerr1']
    rp_err_mean = (abs(r_p_errp) + abs(r_p_errm))/2

    r_s = df['st_rad']
    r_s_errp = df['st_raderr2']
    r_s_errm = df['st_raderr1']
    rs_err_mean = (abs(r_s_errp) + abs(r_s_errm)) / 2

    r_s_new = df[new_r_s_colname]
    r_s_new_errp = df[new_r_s_colname + '_errp']
    r_s_new_errm = df[new_r_s_colname + '_errm']
    r_s_new_err_mean = (abs(r_s_new_errp) + abs(r_s_new_errm)) / 2

    error  =  ( (rp_err_mean / r_p)**2 + (rs_err_mean / r_s)**2 + (r_s_new_err_mean/r_s_new)**2 )**(1/2)
    #error = error*r_p
    #error = error * 100
    df['r_new_error'] = error

    return error

def median(errors):
    median = errors.median()
    return median

def planet_mass(df,m_star_colname, Msini_units='earth'):
    #stole this from radvel - maybe recalculate
    """Calculate Msini

    Calculate Msini for a given K, P, stellar mass, and e

    Args:
        K (float or array: Doppler semi-amplitude [m/s]
        P (float or array): Orbital period [days]
        Mstar (float or array): Mass of star [Msun]
        e (float or array): eccentricity
        Msini_units (Optional[str]): Units of Msini {'earth','jupiter'}
            default: 'earth'

    Returns:
        float or array: Msini [units = Msini_units]

    """
    m_star = np.array(df[m_star_colname])
    # convert inputs to array so they work with units
    P = np.array(df['pl_orbper'])
    Mstar = np.array(m_star)
    K = np.array(df['pl_rvamp'])
    e = np.array(df['pl_orbeccen'])
    G = c.G.value                # added gravitational constant
    Mjup = c.M_jup.value         # added Jupiter's mass
    Msun = c.M_sun.value         # added sun's mass
    Mstar = Mstar*Msun
    Mstar = np.array(Mstar)
    K_0 = 28.4329

    incl = df['pl_orbincl']
    old_mass = df['pl_bmasse']


    P_year = (P * u.d).to(u.year).value
    P = (P * u.d).to(u.second).value

    # First assume that Mp << Mstar
    Msini = K / K_0 * np.sqrt(1.0 - e ** 2.0) * (Mstar/Msun) ** (2.0 / 3.0) * P_year ** (1 / 3.0)

    # Use correct calculation if any elements are >10% of the stellar mass
    if (np.array(((Msini * u.Mjup).to(u.M_sun) / (Mstar/Msun)).value > 0.10)).any():
        warnings.warn("Mpsini << Mstar assumption broken, correcting Msini calculation.")

        a = K*(((2*(np.pi)*G)/P)**(-1/3.))*np.sqrt(1-(e**2))
        Msini = []
        if isinstance(P, float):
            n_elements = 1
        else:
            assert type(K) == type(P) == type(Mstar) == type(e), "All input data types must match."
            assert K.size == P.size == Mstar.size == e.size, "All input arrays must have the same length."
            n_elements = len(P)
        for i in range(n_elements):
            def func(x):
                try:
                    return x - a[i]*((Mstar[i]+x)**(2/3.))
                except IndexError:
                    return x - a * ((Mstar + x) ** (2 / 3.))

            sol = root(func, Mjup)
            Msini.append(sol.x[0])

        Msini = np.array(Msini)
        Msini = Msini/Mjup

    if Msini_units.lower() == 'jupiter':
        pass
    elif Msini_units.lower() == 'earth':
        Msini = (Msini * u.M_jup).to(u.M_earth).value
    else:
        raise Exception("Msini_units must be 'earth', or 'jupiter'")


    mp_recalc = Msini / np.sin(np.deg2rad(incl))

    #fig, ax = plt.subplots()

    #ax.errorbar(old_mass * np.sin(incl), Msini, yerr=m_p_err, c='k', marker='^', ms=4, capsize=3, alpha=0.4)
    #ax.set_xlabel(r'NASA $M_p$')
    #ax.set_ylabel(r'Recalculated $M_p$')
    # plt.xlim()
    # plt.ylim()
    # ax.set_aspect('equal', adjustable='box')
    #plt.show()
    # WROKS ON THIS _ CALCULATE ERROES

    return mp_recalc

def new_mass_err(df, new_m_s_colname):

    mp = df['pl_bmasse']

    K = np.array(df['pl_rvamp'])
    K_errp = df['pl_rvamperr2']
    K_errm = df['pl_rvamperr1']
    K_err_mean = (abs(K_errp) + abs(K_errm)) / 2



    e = np.array(df['pl_orbeccen'])
    df['pl_orbeccenerr2'] = np.where(df.pl_orbeccen == 'nan', 0, df.pl_orbeccenerr2)
    df['pl_orbeccenerr1'] = np.where(df.pl_orbeccen == 'nan', 0, df.pl_orbeccenerr1)
    e_errp = df['pl_orbeccenerr2']
    e_errm = df['pl_orbeccenerr1']

    e_err_mean = (abs(e_errp) + abs(e_errm)) / 2
    rat = e_err_mean / e
    e_err_abs = rat[rat == 'nan'] = 0


    incl = np.array(df['pl_orbincl'])
    incl_errp = df['pl_orbinclerr2']
    incl_errm = df['pl_orbinclerr1']
    incl_err_mean = (abs(incl_errp) + abs(incl_errm)) / 2


    m_p = df['pl_bmasse']
    m_p_errp = df['pl_bmasseerr2']
    m_p_errm = df['pl_bmasseerr1']
    mp_err_mean = (abs(m_p_errp) + abs(m_p_errm))/2

    m_s = df['st_mass']
    m_s_errp = df['st_masserr2']
    m_s_errm = df['st_masserr1']
    ms_err_mean = (abs(m_s_errp) + abs(m_s_errm)) / 2

    m_s_new = df[new_m_s_colname]
    m_s_new_errp = df[new_m_s_colname + '_errp']
    m_s_new_errm = df[new_m_s_colname + '_errm']
    m_s_new_err_mean = (abs(m_s_new_errp) + abs(m_s_new_errm)) / 2

    error  =  ((ms_err_mean/m_s)**2 + (m_s_new_err_mean/m_s_new)**2 + ( e_err_abs )**2 + (K_err_mean / K)**2)**0.5

    #print(e_err_abs)



    #error = error * 100

    df['m_new_error'] = error

    #print(e_err_mean/e, K_err_mean / K)
    return error

def join_input_output(input, output, nasa):
    joined_for_analysis = pd.merge(input, output, on='#starid', suffixes=('_input', '_bg'))
    #bg stands for basta gaia here- this could be changed
    joined_for_analysis = pd.merge(nasa, joined_for_analysis, left_on='gaia_id', right_on='#starid')
    final_idx = (joined_for_analysis['rv_flag'] == 1) * (joined_for_analysis['tran_flag'] == 1) * (
                joined_for_analysis['pl_bmassprov'] == 'Mass') * (joined_for_analysis['pl_rade'] <= 4)
    intermediate = joined_for_analysis[final_idx]
    error_cut_idx = (intermediate['pl_bmasseerr1'] / intermediate['pl_bmasse'] <= 0.25) * (
                intermediate['pl_radeerr1'] / intermediate['pl_rade'] <= 0.25)
    full_analysis = intermediate[error_cut_idx]
    # For now, keeping the shorter planet dataframe in order to analyse all planet parameters
    # however, ;last three steps could be skipped in order to investigate companions / populations
    len(full_analysis)

    return full_analysis


def plot_comparison(ax, var1, var2, output):
    plt.rc('font', size=22)
    plt.rc('axes', titlesize=22)
    plt.rc('legend', fontsize=14)

    ax.scatter(output[var1], output[var2], s=2)
    input_range = max(output[var1])-min(output[var1])
    output_range = max(output[var1])-min(output[var1])
    ax.set_ylabel('NASA')

    ax.axes.xaxis.set_visible(False)
    if 'Teff' in var1:
        ax.set_title('Teff [K]')
    if 'FeH' in var1:
        ax.set_title('FeH [dex]')
    if 'mass' in var1:
        ax.set_title(r'Mass [M$_{\odot}$]')
    if 'rad' in var1:
        ax.set_title(r'Radius [R$_{\odot}$]')
    if 'age' in var1:
        ax.set_title(r'Age [Myr]')

def plot_residual(ax, var1, var2, output):
    residual = output[var1] - output[var2]
    ax.scatter(output[var1], residual, s=2)
    ax.set_xlabel('Gaia')
    ax.set_ylabel('Gaia - NASA')
    range = max(output[var1])-min(output[var1])
    residual_range = max(residual) - min(residual)
    ax.set_aspect(range/(residual_range*4))

def plot(var1, var2, output):

    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols = 1, sharex='all', figsize=(5,10))

    plt.subplots_adjust(hspace=0.1)
    plot_comparison(ax1, var1, var2, output)
    plot_residual(ax2, var1, var2, output)

    plt.show()

def plot_valley_1D(df, rp_new, mp_new, spec_source, mass_split=False):
    plt.rcParams['axes.formatter.min_exponent'] = 4
    #plt.style.use('dark_background')

    plt.rc('font', size=22)
    plt.rc('axes', titlesize=22)
    plt.rc('legend', fontsize=14)

    rp_old = np.array(df['pl_rade'])
    mp_old = np.array(df['pl_bmasse'])

    re = c.R_earth.value
    rj = c.R_jup.value

    rp_new= np.array(rp_new)
    #rp_new = np.array((rp_new * rj) / re)
    mp_new = np.array(mp_new)
    print('For the 1D histogram, there are {} old and {} new radii'.format(len(rp_old), len(rp_new)))

    plt.hist(rp_old, label='Original NASA value', alpha=0.9, bins=10, color = 'red', histtype='step',linewidth=4)
    #plt.hist(rp_old + df['pl_radeerr1'], alpha=0.5, bins=10,color = 'red', histtype='step', linewidth=2, ls='--')
    #plt.hist(rp_old + df['pl_radeerr2'], alpha=0.5, bins=10, color = 'red',histtype='step', linewidth=2, ls='--')


    plt.hist(rp_new, label='BASTA+Gaia', alpha=0.7, bins=10, color='goldenrod', histtype='step', linewidth=4)
    #plt.hist(rp_old - df['r_new_error'], alpha=0.5, bins=10, color='goldenrod', histtype='step', linewidth=2, ls='--')
    #plt.hist(rp_old + df['r_new_error'], alpha=0.5, bins=10, color='goldenrod', histtype='step', linewidth=2, ls='--')
    # plt.xlim(1,2.5)
    plt.xlabel(r'R$_p$ [R$_{\oplus}$]')
    #plt.title('BASTA planet radii')

    plt.legend(loc='lower left')
    plt.tight_layout()

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/radius_valley_1D.png'.format(spec_source))
    plt.show()

    plt.hist(mp_new, label='BASTA+GAIA', alpha=0.9, bins=15, histtype='step')
    plt.hist(mp_old, label='Original NASA value', alpha=0.9, bins=15, histtype='step')
    plt.title('BASTA planet masses')
    plt.xlabel(r'M$_p$ [M$_J$]')
    plt.legend()



    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/mpl_mass_1D.png'.format(spec_source))
    plt.show()

def plot_valley_2D(df, rp_new, spec_source):
    plt.rcParams['axes.formatter.min_exponent'] = 4
    #plt.style.use('dark_background')
    import matplotlib as mpl
    #mpl.rcParams['figure.dpi'] = dpi
    plt.rc('font', size=30)
    plt.rc('axes', titlesize=30)
    plt.rc('legend', fontsize=30)


    rp_old = np.array(df['pl_rade'])
    re = c.R_earth.value
    rj = c.R_jup.value

    #rp_gaia_earth = np.array((rp_new * rj) / re)
    rp_gaia_earth = np.array(rp_new)
    rp_nasa_earth = np.array(rp_old)

    st_mass_old = np.array(df['st_mass'])
    st_mass_new = np.array(df['massfin'])

    period  = np.array(df['pl_orbper'])
    age_old = np.array(df['st_age'])
    age_new = np.array(df['age']/1000)
    # now plotting st_mass valley

    fig, ax = plt.subplots(figsize=(15, 10))

    mass_grad = 0.23
    mass_intercept = 0.27

    # could change this to be bigger errors - have not included different gradients that cynthia used
    x = np.linspace(0.1, 1.4, 1000)
    y = 10 ** (mass_grad * np.log10(x) + mass_intercept)
    y_up = 10 ** ((mass_grad) * np.log10(x) + mass_intercept + 0.01)
    y_lo = 10 ** ((mass_grad) * np.log10(x) + mass_intercept - 0.01)

    plt.grid(True, which="both", alpha=0.3, color='lightgray')
    for i in range(0,len(st_mass_old)):
        plt.arrow(x=st_mass_old[i], y=rp_nasa_earth[i], dx=(st_mass_new[i]-st_mass_old[i]), dy=(rp_gaia_earth[i] - rp_nasa_earth[i]), width=.000001,
                  length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.1,
                  head_starts_at_zero=False)

    plt.scatter((st_mass_new), (rp_gaia_earth), s=40, c='gold', marker='^', label='BASTA+GAIA')
    plt.scatter((st_mass_old), (rp_nasa_earth), s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')
    plt.plot(x, y, label='Ho 22')
    plt.plot(x, y_up, c='lightblue', linestyle='--')
    plt.plot(x, y_lo, c='lightblue', linestyle='--')
    plt.xlabel(r'$log_{10}(M_{star}}$) [M$_\odot$]')
    plt.ylabel(r'$log_{10}$(R$_p$ [R$_\oplus$])')
    plt.xscale('log')
    plt.yscale('log')

    plt.legend()

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/radius_valley_stellar_mass.png'.format(spec_source))

    plt.show()

    # now plotting period valley
    fig, ax = plt.subplots(figsize=(15, 10))

    period_grad = -0.11
    period_intercept = 0.37

    x = np.linspace(0.1, 100, 1000)
    y = 10 ** (period_grad * np.log10(x) + period_intercept)
    y_up = 10 ** ((period_grad) * np.log10(x) + period_intercept + 0.02)
    y_lo = 10 ** ((period_grad) * np.log10(x) + period_intercept - 0.02)

    plt.grid(True, which="both", alpha=0.3, color='lightgray')
    for i in range(len(period)):
        plt.arrow(x=period[i], y=rp_nasa_earth[i], dx=(0), dy=(rp_gaia_earth[i] - rp_nasa_earth[i]), width=.000001,
                  length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.1,
                  head_starts_at_zero=False)

    plt.scatter((period), (rp_gaia_earth), s=40, c='gold', marker='^', label='BASTA+GAIA')
    plt.scatter((period), (rp_nasa_earth), s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')
    plt.plot(x, y, label='Ho 22 valley')
    plt.plot(x, y_up, c='lightblue', linestyle='--')
    plt.plot(x, y_lo, c='lightblue', linestyle='--')
    plt.xlabel(r'$log_{10}$P [d]')
    plt.ylabel(r'$log_{10}$(R$_p$ [R$_E$])')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.legend()

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/radius_valley_period.png'.format(spec_source))
    plt.show()

    #now plotting age valley
    fig, ax = plt.subplots(figsize=(15, 10))

    age_grad = 0.02
    age_intercept = 0.26

    x = np.linspace(1, 20, 1000)
    y = 10 ** (age_grad * np.log10(x) + age_intercept)
    y_up = 10 ** ((age_grad) * np.log10(x) + age_intercept + 0.01)
    y_lo = 10 ** ((age_grad) * np.log10(x) + age_intercept - 0.01)

    plt.grid(True, which="both", alpha=0.3, color='lightgray')
    for i in range(len(age_old)):
        plt.arrow(x=age_old[i], y=rp_nasa_earth[i], dx=(age_new[i]-age_old[i]), dy=(rp_gaia_earth[i] - rp_nasa_earth[i]), width=.000001,
                  length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.1,
                  head_starts_at_zero=False)

    plt.scatter(age_new, rp_gaia_earth, s=40, c='gold', marker='^', label='BASTA+GAIA')
    plt.scatter((age_old), (rp_nasa_earth), s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')
    plt.plot(x, y, label='Ho 22 valley')
    plt.plot(x, y_up, c='lightblue', linestyle='--')
    plt.plot(x, y_lo, c='lightblue', linestyle='--')
    plt.xlabel(r'$log_{10}(Age_{star}}$) [Gyr]')
    plt.ylabel(r'$log_{10}$(R$_p$ [R$_E$])')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1,max(age_new))

    plt.legend()

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/radius_valley_stellar_age_cut.png'.format(spec_source))
    plt.show()

def new_param_df(df, rp_new, mp_new):
    df['rp_new']= rp_new
    df['mp_new'] = mp_new

    return df

def plot_mass_radius(df, rp_new, mp_new, spec_source):
    plt.rcParams['axes.formatter.min_exponent'] = 4
    import matplotlib as mpl
    #plt.style.use('dark_background')
    #mpl.rcParams['figure.dpi'] = dpi
    plt.rc('font', size=22)
    plt.rc('axes', titlesize=22)
    plt.rc('legend', fontsize=14)

    #plt.style.use('dark_background')
    rp_old = np.array(df['pl_rade'])
    re = c.R_earth.value
    rj = c.R_jup.value
    mj = c.M_jup.value

    #rp_new = np.array((rp_new * rj) / re)
    rp_new = np.array(rp_new)
    rp_old = np.array(rp_old)

    mp_old = np.array(df['pl_bmasse'])
    mp_new = np.array(mp_new)

    ms = np.array(df['massfin'])


    fig, ax = plt.subplots(figsize=(15, 10))

    #PLotting planet composition thoretical tracks

    #importing the values from files

    fe_100 = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_fe.ascii', delim_whitespace=True)
    earth = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/earth.ascii', delim_whitespace=True)
    rock = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_rock.ascii', delim_whitespace=True)
    h_cold = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/H_cold.ascii', delim_whitespace=True)
    max_coll = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/max_coll_strip.ascii',
                             delim_whitespace=True)
    h2o_100_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/h20_100_500k.ascii',
                                 delim_whitespace=True)
    h2_2_earth_49_h2o_49_500k = pd.read_table(
        '/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_49_earth_49_h2O_500k.ascii', delim_whitespace=True)
    h2_2_earth_98_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_98_earth_500k.ascii',
                                       delim_whitespace=True)

    plt.plot(fe_100['Mass'], fe_100['Radius'], c='orange', label='100% Fe', alpha=0.4)
    plt.plot(earth['Mass'], earth['Radius'], c='green', label=f'32.5% Fe + 67.5% MgSiO$_3$', alpha=0.4)
    plt.plot(rock['Mass'], rock['Radius'], c='brown', label=f'100% MgSiO$_3$', alpha=0.4)
    plt.plot(h2_2_earth_49_h2o_49_500k['Mass'], h2_2_earth_49_h2o_49_500k['Radius'], linestyle='-',
             label=f'49% Earth-like rocky core  \n+ 49% H$_2$O + \n 2% H$_2$ atmosphere', alpha=0.4)
    plt.plot(h2_2_earth_98_500k['Mass'], h2_2_earth_98_500k['Radius'], linestyle='-',
             label=f'98% Earth-like rocky core \n+ 2% H$_2$ atmosphere', alpha=0.4)
    # plt.plot(h_cold['Mass'], h_cold['Radius'], label  = 'Cold Hydrogen')
    plt.plot(max_coll['Mass'], max_coll['Radius'], linestyle='--', c='lightgrey', label='Maximum collisional\n stripping',
             alpha=0.5)
    plt.plot(h2o_100_500k['Mass'], h2o_100_500k['Radius'], linestyle='-', c='cyan', label=f'100% H$_2$O', alpha=0.5)


    plt.grid(True, which="both", alpha=0.3, color='lightgray')


    for i in range(len(ms)):
        plt.arrow(x=mp_old[i], y=rp_old[i], dx=(mp_new[i] - mp_old[i]),
                  dy=(rp_new[i] - rp_old[i]), width=.000001,
                  length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.14,
                  head_starts_at_zero=False)

    plt.scatter(mp_new, rp_new, s=50, marker='^', label='BASTA+GAIA', c=ms, cmap='viridis')
    plt.scatter(mp_old, rp_old, s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')

    plt.legend()
    plt.clim(0, 1.5)
    plt.xlim(0.1, 60)
    plt.ylim(0, 5)
    plt.colorbar(label=r'Stellar Mass M$_{\odot}$')

    plt.xlabel(r'M$_p$ [M${\oplus}$]')
    plt.ylabel(r'R$_p$ [R$_{\oplus}$]')
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend(loc = 'upper left')

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/mass_radius_comp.png'.format(spec_source))

    plt.show()

def plot_mass_radius_errors(df, rp_new, mp_new, spec_source):
    plt.rcParams['axes.formatter.min_exponent'] = 4
    import matplotlib as mpl
    #plt.style.use('dark_background')
    #mpl.rcParams['figure.dpi'] = dpi
    plt.rc('font', size=22)
    plt.rc('axes', titlesize=22)
    plt.rc('legend', fontsize=14)

    #plt.style.use('dark_background')
    rp_old = np.array(df['pl_rade'])
    re = c.R_earth.value
    rj = c.R_jup.value
    mj = c.M_jup.value

    #rp_new = np.array((rp_new * rj) / re)
    rp_new = np.array(rp_new)
    rp_old = np.array(rp_old)

    mp_old = np.array(df['pl_bmasse'])
    mp_new = np.array(mp_new)

    ms = np.array(df['massfin'])


    fig, ax = plt.subplots(figsize=(15, 10))

    #PLotting planet composition thoretical tracks

    #importing the values from files

    fe_100 = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_fe.ascii', delim_whitespace=True)
    earth = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/earth.ascii', delim_whitespace=True)
    rock = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_rock.ascii', delim_whitespace=True)
    h_cold = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/H_cold.ascii', delim_whitespace=True)
    max_coll = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/max_coll_strip.ascii',
                             delim_whitespace=True)
    h2o_100_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/h20_100_500k.ascii',
                                 delim_whitespace=True)
    h2_2_earth_49_h2o_49_500k = pd.read_table(
        '/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_49_earth_49_h2O_500k.ascii', delim_whitespace=True)
    h2_2_earth_98_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_98_earth_500k.ascii',
                                       delim_whitespace=True)

    plt.plot(fe_100['Mass'], fe_100['Radius'], c='orange', label='100% Fe', alpha=0.4)
    plt.plot(earth['Mass'], earth['Radius'], c='green', label=f'32.5% Fe + 67.5% MgSiO$_3$', alpha=0.4)
    plt.plot(rock['Mass'], rock['Radius'], c='brown', label=f'100% MgSiO$_3$', alpha=0.4)
    plt.plot(h2_2_earth_49_h2o_49_500k['Mass'], h2_2_earth_49_h2o_49_500k['Radius'], linestyle='-',
             label=f'49% Earth-like rocky core  \n+ 49% H$_2$O + \n 2% H$_2$ atmosphere', alpha=0.4)
    plt.plot(h2_2_earth_98_500k['Mass'], h2_2_earth_98_500k['Radius'], linestyle='-',
             label=f'98% Earth-like rocky core \n+ 2% H$_2$ atmosphere', alpha=0.4)
    # plt.plot(h_cold['Mass'], h_cold['Radius'], label  = 'Cold Hydrogen')
    plt.plot(max_coll['Mass'], max_coll['Radius'], linestyle='--', c='lightgrey', label='Maximum collisional\n stripping',
             alpha=0.5)
    plt.plot(h2o_100_500k['Mass'], h2o_100_500k['Radius'], linestyle='-', c='cyan', label=f'100% H$_2$O', alpha=0.5)


    plt.grid(True, which="both", alpha=0.3, color='lightgray')

    yerr = np.array(df['r_new_error'])
    xerr = np.array(df['m_new_error'])


    for i in range(len(ms)):
        plt.arrow(x=mp_old[i], y=rp_old[i], dx=(mp_new[i] - mp_old[i]),
                  dy=(rp_new[i] - rp_old[i]), width=.000001,
                  length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.14,
                  head_starts_at_zero=False)

    plt.errorbar(mp_new, rp_new, xerr=xerr, yerr=yerr, ms=5, marker='^', label='BASTA+GAIA', c='blue', ecolor = 'gray', ls='none')
    plt.scatter(mp_old, rp_old, s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')

    plt.legend()
    #plt.clim(0, 1.5)
    plt.xlim(0.1, 60)
    plt.ylim(0, 5)
    #plt.colorbar(label=r'Stellar Mass M$_{\odot}$')

    plt.xlabel(r'M$_p$ [M${\oplus}$]')
    plt.ylabel(r'R$_p$ [R$_{\oplus}$]')
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend(loc = 'upper left')

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/mass_radius_errors_comp.png'.format(spec_source))

    plt.show()
def mass_split(df, colname, bins=[0,0.3,0.7,0.8,0.9,1,1.1]):

    #splits the data into mass bins of equal widdth
    #this function is bad lol


    masses = plt.hist(df[colname], bins=bins)
    bin_edges = masses[1]
    plt.xlabel(r'BASTA + GAIA stellar mass [M$_{\odot}$]')
    plt.show()

    masscuts = dict()

    for bin in range(len(bin_edges) - 1):
        masscuts['m_{}'.format(bin)] = df[
            (df[colname] >= bin_edges[bin]) & (df[colname] < bin_edges[bin + 1])]
    #print((masscuts.keys()))

    for i in range(0, len(masscuts)):
        plt.hist(masscuts['m_{}'.format(i)]['pl_rade'], histtype='step',
                    label=r'{} < M$_s$ < {} [M$_\odot$]'.format(bin_edges[i], bin_edges[i + 1]))
        plt.legend()
    plt.show()

    #print(bin_edges)
    #print(max(df[colname]))
    return masscuts

def m_r_plot_masses(masses_split, mp_new, rp_new, spec_source):

    ######################MASS RADIUS DIAGRAM##########################################################################

    plt.rcParams.update({'font.size': 25})
    plt.rcParams.update({'xtick.major.size': 25})
    plt.rcParams.update({'xtick.labelsize': 25})
    plt.rcParams.update({'ytick.major.size': 25})
    plt.rcParams.update({'ytick.labelsize': 25})
    plt.rcParams['axes.formatter.min_exponent'] = 4

    fe_100 = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_fe.ascii', delim_whitespace=True)
    earth = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/earth.ascii', delim_whitespace=True)
    rock = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/100_rock.ascii', delim_whitespace=True)
    h_cold = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/H_cold.ascii', delim_whitespace=True)
    max_coll = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/max_coll_strip.ascii',
                             delim_whitespace=True)
    h2o_100_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/h20_100_500k.ascii',
                                 delim_whitespace=True)
    h2_2_earth_49_h2o_49_500k = pd.read_table(
        '/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_49_earth_49_h2O_500k.ascii', delim_whitespace=True)
    h2_2_earth_98_500k = pd.read_table('/Users/afw2/BASTA/Paper1/data/composition_tracks/2_h2_98_earth_500k.ascii',
                                       delim_whitespace=True)

    fig, ax = plt.subplots(2,3, figsize=(25, 17))
    plt.subplots_adjust(hspace=0.4)
    plt.subplots_adjust(wspace=0.4)

    for i in masses_split:
        m_i = '{}'.format(i)
        intm = int(i.strip('m_'))
        if intm <= 2:
            i_x = 0
            i_y = intm
        else:
            i_x = 1
            i_y = intm-3

        #print([i_x, i_y])

        mp_old = np.array(masses_split[m_i]['pl_bmasse'])
        rp_old = np.array(masses_split[m_i]['pl_rade'])
        st_mass_new = np.array(masses_split[m_i]['massfin'])

        rp_new = np.array(masses_split[m_i]['rp_new'])
        mp_new = np.array(masses_split[m_i]['mp_new'])

        #print(len(st_mass_new))

        for j in range(0,len(st_mass_new)):

            #ax[i_x, i_y].legend()
            ax[i_x, i_y].set_xlabel(r'M$_p$ [M${\oplus}$]')
            ax[i_x, i_y].set_ylabel(r'R$_p$ [R$_{\oplus}$]')  #plt.xscale('log')
            #ax[i_x, i_y].set_xscale('log')
            ax[i_x, i_y].grid(True, which="both", alpha=0.3, color='lightgray')
            ax[i_x, i_y].set_title(r'{} < M$_s$ [M$_\odot$] < {}'.format(min(st_mass_new.round(2)), max(st_mass_new.round(2))))


            ax[i_x, i_y].plot(fe_100['Mass'], fe_100['Radius'], c = 'orange', label  = '100% Fe', alpha = 0.1)
            ax[i_x, i_y].plot(earth['Mass'], earth['Radius'], c = 'green',  label  = f'32.5% Fe + 67.5% MgSiO$_3$', alpha = 0.1)
            ax[i_x, i_y].plot(rock['Mass'], rock['Radius'],c = 'brown', label  = f'100% MgSiO$_3$', alpha = 0.4)
            ax[i_x, i_y].plot(h2_2_earth_49_h2o_49_500k['Mass'], h2_2_earth_49_h2o_49_500k['Radius'], linestyle='-', label  = f'49% Earth-like rocky core + 49% H$_2$O + 2% H$_2$ atmosphere', alpha = 0.4)
            ax[i_x, i_y].plot(h2_2_earth_98_500k['Mass'], h2_2_earth_98_500k['Radius'], linestyle='-', label  = f'98% Earth-like rocky core + 2% H$_2$ atmosphere', alpha = 0.1)
            #plt.plot(h_cold['Mass'], h_cold['Radius'], label  = 'Cold Hydrogen')
            ax[i_x, i_y].plot(max_coll['Mass'], max_coll['Radius'], linestyle='--', c='lightgrey', label  = 'Maximum collisional stripping', alpha = 0.1)
            ax[i_x, i_y].plot(h2o_100_500k['Mass'], h2o_100_500k['Radius'], linestyle='-', c='cyan', label  = f'100% H$_2$O', alpha = 0.1)

            ax[i_x, i_y].arrow(x=mp_old[j], y=rp_old[j], dx=(mp_new[j] - mp_old[j]), dy=(rp_new[j] - rp_old[j]),
                               width=.000001,
                               length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.1,
                               head_starts_at_zero=False)

            ax[i_x, i_y].scatter(mp_new[j], rp_new[j], s=70, marker='^', label='BASTA+GAIA', c=st_mass_new[j],
                                 cmap='viridis')
            ax[i_x, i_y].scatter(mp_old[j], rp_old[j], s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')




            ax[i_x, i_y].set_xscale('log')
            ax[i_x, i_y].set_xlim(0.15,60)
            ax[i_x, i_y].set_ylim(0.5,4)

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/mass_split_mass_radius_comp.png'.format(spec_source))
    plt.show()

    ######################2DHIST DIAGRAM##########################################################################

    plt.rcParams.update({'font.size': 25})
    # plt.style.use('dark_background')
    plt.rcParams.update({'xtick.major.size': 25})
    plt.rcParams.update({'xtick.labelsize': 25})
    plt.rcParams.update({'ytick.major.size': 25})
    plt.rcParams.update({'ytick.labelsize': 25})
    plt.rcParams['axes.formatter.min_exponent'] = 4

    fig, ax = plt.subplots(2, 3, figsize=(25, 17))
    plt.subplots_adjust(hspace=0.4)
    plt.subplots_adjust(wspace=0.4)
    for i in masses_split:
        m_i = '{}'.format(i)
        intm = int(i.strip('m_'))
        if intm <= 2:
            i_x = 0
            i_y = intm
        else:
            i_x = 1
            i_y = intm - 3

        #print([i_x, i_y])

        mp_old = np.array(masses_split[m_i]['pl_bmasse'])
        rp_old = np.array(masses_split[m_i]['pl_rade'])
        period = np.array(masses_split[m_i]['pl_orbper'])
        rp_new = np.array(masses_split[m_i]['rp_new'])
        st_mass_new = np.array(masses_split[m_i]['massfin'])

        for j in range(0, len(st_mass_new)):
            # ax[i_x, i_y].legend()
            ax[i_x, i_y].set_xlabel(r'P [d]]')
            ax[i_x, i_y].set_ylabel(r'R$_p$ [R$_{\oplus}$]')

            ax[i_x, i_y].grid(True, which="both", alpha=0.3, color='lightgray')
            ax[i_x, i_y].set_title(r'{} < M$_s$ [M$_\odot$] < {}'.format(min(st_mass_new.round(2)), max(st_mass_new.round(2))))

            period_grad = -0.11
            period_intercept = 0.37

            x = np.linspace(0.1, 100, 1000)
            y = 10 ** (period_grad * np.log10(x) + period_intercept)
            y_up = 10 ** ((period_grad) * np.log10(x) + period_intercept + 0.02)
            y_lo = 10 ** ((period_grad) * np.log10(x) + period_intercept - 0.02)

            plt.grid(True, which="both", alpha=0.3, color='lightgray')
            for i in range(len(period)):
                ax[i_x, i_y].arrow(x=period[i], y=rp_old[i], dx=(0), dy=(rp_new[i] - rp_old[i]),
                          width=.000001,
                          length_includes_head=True, head_width=.0009, head_length=.0009, facecolor='k', alpha=0.1,
                          head_starts_at_zero=False)

            ax[i_x, i_y].scatter((period), (rp_new), s=40, c='gold', marker='^', label='BASTA+GAIA')
            ax[i_x, i_y].scatter((period), (rp_new), s=30, c='lightgrey', alpha=0.4, marker='o', label='NASA')
            ax[i_x, i_y].plot(x, y, label='Ho 22 valley', c = 'lightblue')
            ax[i_x, i_y].plot(x, y_up, c='lightblue', linestyle='--')
            ax[i_x, i_y].plot(x, y_lo, c='lightblue', linestyle='--')
            ax[i_x, i_y].set_xscale('log')
            ax[i_x, i_y].set_yscale('log')
            ax[i_x, i_y].set_xlim(0.7,100)

            ax[i_x, i_y].set_ylim(0.7,4)



    #plt.legend()
    plt.savefig('/Users/afw2/BASTA/Paper1/figures/{}/mass_split_period_radius_comp.png'.format(spec_source))
    plt.show()


def sample_visualisation(df):

    plt.rcParams.update({'font.size': 25})
    # plt.style.use('dark_background')
    plt.rcParams.update({'xtick.major.size': 25})
    plt.rcParams.update({'xtick.labelsize': 25})
    plt.rcParams.update({'ytick.major.size': 25})
    plt.rcParams.update({'ytick.labelsize': 25})
    plt.rcParams['axes.formatter.min_exponent'] = 4



    df = df.drop_duplicates(subset='pl_name')

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize = (23,7))
    plt.subplots_adjust(hspace=0.4)
    plt.subplots_adjust(wspace=0.4)


    xerr = [-df['pl_bmasseerr2'], df['pl_bmasseerr1']]
    yerr = [-df['pl_radeerr2'], df['pl_radeerr1']]

    ax[0].errorbar(df['pl_bmasse'], df['pl_rade'], fmt='o', c='green', xerr=xerr, yerr=yerr,
                   alpha=0.5, elinewidth=0.4)
    ax[0].loglog()
    ax[0].set_ylabel(r'R$_{\rm{p}}$ [R$_{\oplus}$]')
    ax[0].set_xlabel(r'M$_{\rm{p}}$ [M$_{\oplus}$]')

    xerr = [-df['st_tefferr2'], df['st_tefferr1']]
    yerr = [-df['st_loggerr2'], df['st_loggerr1']]

    ax[1].errorbar(df['st_teff'], df['st_logg'], fmt='o', c='green', xerr=xerr, yerr=yerr,
                   alpha=0.5, elinewidth=0.4)

    ax[1].set_xlabel(r'T$_{\rm{eff}}$ [K]')
    ax[1].set_ylabel(r'log(g) [dex]')

    ax[1].grid(alpha=0.5)
    ax[1].set_xlim(7000, 3000)
    ax[1].set_ylim(5.5, 3.5)

    ax[2].hist(df['ks_m'], label='Ks Mag', alpha=0.6, histtype='step')
    ax[2].hist(df['h_m'], label='H Mag', alpha=0.6, histtype='step')
    ax[2].hist(df['j_m'], label='J Mag', alpha=0.6, histtype='step')

    ax[2].set_xlabel('Magnitude')
    ax[2].legend(fontsize=16, loc = 'upper right')

    plt.tight_layout()

    plt.savefig('/Users/afw2/BASTA/Paper1/figures/sample_full.pdf')
    plt.show()