"""
Alpha Q galaxy catalog class.
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
    Lu galaxy catalog class. Uses generic quantity and filter mechanisms
    defined by BaseGenericCatalog class.
    """

    def _subclass_init(self, filename, **kwargs):

        assert os.path.isfile(filename), 'Catalog file {} does not exist'.format(filename)
        self._file = filename
        self.lightcone = kwargs.get('lightcone')
        cosmology      = kwargs.get('cosmology')
        metadata       = kwargs.get('metadata')

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
            'gal_id' :            'gal_id',
            'ra':                 'ra',
            'dec':                'dec',
            'ra_true':            'ra',
            'dec_true':           'dec',
            'redshift':           'redshift',
            'redshift_true':      'z_nopec',
            'halo_id':            'halo_id_nbody',
            'halo_mass':          'm_vir',
            'is_central':         (lambda x : ~x.astype(np.bool), 'gal_type'),
            'stellar_mass':       'mstar',
            'size_disk_true':     'r_disk',
            'size_bulge_true':    'r_bulge',
        }

        #apparent mags
        for band in 'ugriz':
            self._quantity_modifiers['mag_{}_sdss'.format(band)] = 'sdss_{}'.format(band)

        for band in ['435w', '606w', '775w', '814w', '850lp']:
            self._quantity_modifiers['mag_{}_acs'.format(band)] = 'acsf{}'.format(band)

        for band in ['275w', '336w', '105w', '125w', '160w']:
            self._quantity_modifiers['mag_{}_wfc3'.format(band)] = 'wfc3f{}'.format(band)

        for band in ['J', 'H', 'K']:
            self._quantity_modifiers['mag_{}_ukirt'.format(band)]  = 'UKIRT_{}'.format(band)

        for band in ['ch1', 'ch2']:
            self._quantity_modifiers['mag_{}_irac'.format(band)]  = 'irac_{}'.format(band)

        for band in ['FUV', 'NUV']:
            self._quantity_modifiers['mag_{}_galex'.format(band)]  = 'galex_{}'.format(band)

        self._quantity_modifiers['mag_u_ctio'] = 'ctio_U'
        self._quantity_modifiers['mag_u_cfhtls'] = 'CFHTLS_u'
        self._quantity_modifiers['mag_u38_musyc'] = 'musyc_u38'

        #rest frame mags
        for band in ['UV1500', 'UV2300', 'UV2800', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']:
            self._quantity_modifiers['Mag_{}_z0'.format(band)] = '{}_rest'.format(band)


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




