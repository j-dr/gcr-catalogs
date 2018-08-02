"""
CANDELS SAM lightcone galaxy catalog class.
"""
from __future__ import division
import os
import numpy as np
import h5py
from astropy.cosmology import FlatLambdaCDM
from GCR import BaseGenericCatalog

__all__ = ['CANDELSGalaxyCatalog']

class CANDELSGalaxyCatalog(BaseGenericCatalog):
    """
    CANDELS galaxy catalog class. Uses generic quantity and filter mechanisms
    defined by BaseGenericCatalog class.
    """

    def _subclass_init(self, filename, **kwargs):

        assert os.path.isfile(filename), 'Catalog file {} does not exist'.format(filename)
        self._file = filename
        self.lightcone = kwargs.get('lightcone')
        self.area      = float(kwargs.get('area'))

        cosmology      = kwargs.get('cosmology')
        metadata       = kwargs.get('metadata')
        self.candelized     = bool(kwargs.get('candelized', False))

        with h5py.File(self._file, 'r') as fh:
            self.cosmology = FlatLambdaCDM(
                H0=cosmology['H0'],
                Om0=cosmology['Om0'],
                Ob0=cosmology['Ob0']
            )
            try:
                catalog_version = '{}.{}'.format(
                    metadata['versionMajor'],
                    metadata['versionMinor'],
                )
            except KeyError:
                #If no version is specified, it's version 1.0
                catalog_version = '1.0'

        self._quantity_modifiers = {
            'gal_id' :            'truth/gal_id',
            'ra':                 'truth/ra',
            'dec':                'truth/dec',
            'ra_true':            'truth/ra',
            'dec_true':           'truth/dec',
            'redshift':           'truth/redshift',
            'redshift_true':      'truth/z_nopec',
            'halo_id':            'truth/halo_id_nbody',
            'halo_mass':          'truth/m_vir',
            'is_central':         (lambda x : ~x.astype(np.bool), 'truth/gal_type'),
            'stellar_mass':       'truth/mstar',
            'stellar_mass_bulge': 'truth/mbulge',
            'star_formation_rate': 'truth/sfr_ave',
            'size_disk_true':     'truth/r_disk',
            'size_bulge_true':    'truth/r_bulge',
        }

        #apparent mags
        for band in 'ugriz':
            self._quantity_modifiers['mag_true_dust_{}_sdss'.format(band)] = 'truth/sdss_{}_dust'.format(band)
            self._quantity_modifiers['mag_dust_nodust_{}_sdss'.format(band)] = 'truth/sdss_{}'.format(band)

        for band in ['435w', '606w', '775w', '814w', '850lp']:
            self._quantity_modifiers['mag_true_dust_f{}_acs'.format(band)] = 'truth/acsf{}_dust'.format(band)
            self._quantity_modifiers['mag_true_nodust_f{}_acs'.format(band)] = 'truth/acsf{}'.format(band)
            if self.candelized:
                self._quantity_modifiers['mag_f{}_acs'.format(band)] = 'candelized/acsf{}_dust_mag_noisy'.format(band) 
                self._quantity_modifiers['mag_nodust_f{}_acs'.format(band)] = 'candelized/acsf{}_mag_noisy'.format(band)

                self._quantity_modifiers['mag_err_f{}_acs'.format(band)] = 'candelized/acsf{}_dust_magerr'.format(band) 
                self._quantity_modifiers['mag_err_nodust_f{}_acs'.format(band)] = 'candelized/acsf{}_magerr'.format(band)

                self._quantity_modifiers['flux_f{}_acs'.format(band)] = 'candelized/acsf{}_dust_flux_noisy'.format(band) 
                self._quantity_modifiers['flux_nodust_f{}_acs'.format(band)] = 'candelized/acsf{}_flux_noisy'.format(band)

                self._quantity_modifiers['flux_err_f{}_acs'.format(band)] = 'candelized/acsf{}_dust_fluxerr'.format(band) 
                self._quantity_modifiers['flux_err_nodust_f{}_acs'.format(band)] = 'candelized/acsf{}_fluxerr'.format(band)


        for band in ['275w', '336w', '105w', '125w', '160w']:
            self._quantity_modifiers['mag_true_dust_f{}_wfc3'.format(band)] = 'truth/wfc3f{}_dust'.format(band)
            self._quantity_modifiers['mag_true_nodust_f{}_wfc3'.format(band)] = 'truth/wfc3f{}'.format(band)
            self._quantity_modifiers['mag_true_nodust_bulge_f{}_wfc3'.format(band)] = 'truth/wfc3f{}_bulge'.format(band)

            if (self.candelized) & (band != '336w'):
                self._quantity_modifiers['mag_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_dust_mag_noisy'.format(band) 
                self._quantity_modifiers['mag_nodust_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_mag_noisy'.format(band)

                self._quantity_modifiers['mag_err_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_dust_magerr'.format(band) 
                self._quantity_modifiers['mag_err_nodust_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_magerr'.format(band)

                self._quantity_modifiers['flux_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_dust_flux_noisy'.format(band) 
                self._quantity_modifiers['flux_nodust_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_flux_noisy'.format(band)

                self._quantity_modifiers['flux_err_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_dust_fluxerr'.format(band)
                self._quantity_modifiers['flux_err_nodust_f{}_wfc3'.format(band)] = 'candelized/wfc3f{}_fluxerr'.format(band)


        for band in ['J', 'H', 'K']:
            self._quantity_modifiers['mag_true_dust_{}_ukirt'.format(band)]  = 'truth/UKIRT_{}_dust'.format(band)
            self._quantity_modifiers['mag_true_nodust_{}_ukirt'.format(band)]  = 'truth/UKIRT_{}'.format(band)

        for band in ['ch1', 'ch2']:
            self._quantity_modifiers['mag_true_dust_{}_irac'.format(band)]  = 'truth/irac_{}_dust'.format(band)
            self._quantity_modifiers['mag_true_nodust_{}_irac'.format(band)]  = 'truth/irac_{}'.format(band)

            if self.candelized:
                self._quantity_modifiers['mag_{}_irac'.format(band)] = 'candelized/irac_{}_dust_mag_noisy'.format(band) 
                self._quantity_modifiers['mag_nodust_{}_irac'.format(band)] = 'candelized/irac_{}_mag_noisy'.format(band)

                self._quantity_modifiers['mag_err_{}_irac'.format(band)] = 'candelized/irac_{}_dust_magerr'.format(band) 
                self._quantity_modifiers['mag_err_nodust_{}_irac'.format(band)] = 'candelized/irac_{}_magerr'.format(band)

                self._quantity_modifiers['flux_{}_irac'.format(band)] = 'candelized/irac_{}_dust_flux_noisy'.format(band) 
                self._quantity_modifiers['flux_nodust_{}_irac'.format(band)] = 'candelized/irac_{}_flux_noisy'.format(band)

                self._quantity_modifiers['flux_err_{}_irac'.format(band)] = 'candelized/irac_{}_dust_fluxerr'.format(band) 
                self._quantity_modifiers['flux_err_nodust_{}_irac'.format(band)] = 'candelized/irac_{}_fluxerr'.format(band)


        for band in ['FUV', 'NUV']:
            self._quantity_modifiers['mag_true_dust_{}_galex'.format(band)]  = 'truth/galex_{}_dust'.format(band)
            self._quantity_modifiers['mag_true_nodust_{}_galex'.format(band)]  = 'truth/galex_{}'.format(band)

        self._quantity_modifiers['mag_true_dust_u_ctio'] = 'truth/ctio_U_dust'
        self._quantity_modifiers['mag_true_nodust_u_ctio'] = 'truth/ctio_U'

        self._quantity_modifiers['mag_true_dust_u_cfhtls'] = 'truth/CFHTLS_u_dust'
        self._quantity_modifiers['mag_true_nodust_u_cfhtls'] = 'truth/CFHTLS_u'

        self._quantity_modifiers['mag_true_dust_u38_musyc'] = 'truth/musyc_u38_dust'
        self._quantity_modifiers['mag_true_nodust_u38_musyc'] = 'truth/musyc_u38'

        #rest frame mags
        for band in ['UV1500', 'UV2300', 'UV2800', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']:
            self._quantity_modifiers['Mag_true_dust_{}_z0'.format(band)] = 'truth/{}_rest_dust'.format(band)
            self._quantity_modifiers['Mag_true_nodust_{}_z0'.format(band)] = 'truth/{}_rest'.format(band)

        if self.candelized:
            self._quantity_modifiers['p_detection'] = 'candelized/dp'
            self._quantity_modifiers['redshift_phot'] = 'candelized/z_est'

    def _generate_native_quantity_list(self):
        with h5py.File(self._file, 'r') as fh:
            hgroup = fh['galaxyProperties']
            hobjects = []
            #get all the names of objects in this tree
            hgroup.visit(hobjects.append)
            #filter out the group objects and keep the dataste objects
            hdatasets = [hobject for hobject in hobjects if type(hgroup[hobject]) == h5py.Dataset]
            native_quantities = set(hdatasets)
        return native_quantities


    def _iter_native_dataset(self, native_filters=None):
        assert not native_filters, '*native_filters* is not supported'
        with h5py.File(self._file, 'r') as fh:
            def native_quantity_getter(native_quantity):
                return fh['galaxyProperties/{}'.format(native_quantity)].value
            yield native_quantity_getter

    def _get_native_quantity_info_dict(self, quantity, default=None):
        raise(NotImplementedError)


    def _get_quantity_info_dict(self, quantity, default=None):
        raise(NotImplementedError)
        #TODO needs some fixing
        # print "in get quantity"
        # native_name = None
        # if quantity in self._quantity_modifiers:
        #     print "in quant modifers"
        #     q_mod = self._quantity_modifiers[quantity]
        #     if isinstance(q_mod,(tuple,list)):
        #         print "it's a list object, len:",len(length)

        #         if(len(length) > 2):
        #             return default #This value is composed of a function on
        #             #native quantities. So we have no idea what the units are
        #         else:
        #             #Note: This is just a renamed column.
        #             return self._get_native_quantity_info_dict(q_mod[1],default)
        #     else:
        #         print "it's a string: ",q_mod
        #         return self._get_native_quantity_info_dict(q_mod,default)
        # elif quantity in self._native_quantities:
        #     print "in get native quant"
        #     return self._get_native_quantity_info_dict(quantity,default)




