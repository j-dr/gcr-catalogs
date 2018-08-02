"""
UniverseMachine lightcone galaxy catalog class.
"""
from __future__ import division
import os
import numpy as np
import h5py
from astropy.cosmology import FlatLambdaCDM
from GCR import BaseGenericCatalog

__all__ = ['CANDELSGalaxyCatalog']

class UniverseMachineGalaxyCatalog(BaseGenericCatalog):
    """
    UniverseMaching galaxy catalog class. Uses generic quantity and filter mechanisms
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
            'gal_id' :                     'truth/ID',
            'position_x':                  'truth/X',
            'position_y':                  'truth/Y',
            'position_z':                  'truth/Z',
            'redshift':                    'truth/Zlos',
            'redshift_true':               'truth/Zcosmo',
            'halo_id':                     'truth/ID',
            'halo_mass':                   'truth/M',
            'is_central':                  (lambda x : x==-1, 'truth/UPID'),
            'stellar_mass':                (lambda x : x / 10**10, 'truth/obs_SM'),
            'star_formation_rate':         'truth/obs_SFR',
            'stellar_mass_true':           (lambda x : x / 10**10, 'truth/SM'),
            'star_formation_rate_true':    'truth/SFR',
            'specific_star_formation_rate':'truth/SSFR',
            'Mag_true_nodust_UV1500_z0' :  (lambda x, y : x - y, 'truth/obs_UV', 'truth/A_UV'),
            'Mag_true_dust_UV1500_z0' :    'truth/obs_UV',
        }

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
