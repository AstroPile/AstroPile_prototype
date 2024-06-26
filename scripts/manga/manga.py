
import itertools

import datasets
import h5py

from datasets import Features, Value, Array2D
from datasets.data_files import DataFilesPatternsDict


# TODO: Add BibTeX citation
# Find for instance the citation on arxiv or on the dataset repo/website
_CITATION = """\
@ARTICLE{2015ApJ...798....7B,
       author = {{Bundy}, Kevin and {Bershady}, Matthew A. and {Law}, David R. and {Yan}, Renbin and {Drory}, Niv and {MacDonald}, Nicholas and {Wake}, David A. and {Cherinka}, Brian and {S{\'a}nchez-Gallego}, Jos{\'e} R. and {Weijmans}, Anne-Marie and {Thomas}, Daniel and {Tremonti}, Christy and {Masters}, Karen and {Coccato}, Lodovico and {Diamond-Stanic}, Aleksandar M. and {Arag{\'o}n-Salamanca}, Alfonso and {Avila-Reese}, Vladimir and {Badenes}, Carles and {Falc{\'o}n-Barroso}, J{\'e}sus and {Belfiore}, Francesco and {Bizyaev}, Dmitry and {Blanc}, Guillermo A. and {Bland-Hawthorn}, Joss and {Blanton}, Michael R. and {Brownstein}, Joel R. and {Byler}, Nell and {Cappellari}, Michele and {Conroy}, Charlie and {Dutton}, Aaron A. and {Emsellem}, Eric and {Etherington}, James and {Frinchaboy}, Peter M. and {Fu}, Hai and {Gunn}, James E. and {Harding}, Paul and {Johnston}, Evelyn J. and {Kauffmann}, Guinevere and {Kinemuchi}, Karen and {Klaene}, Mark A. and {Knapen}, Johan H. and {Leauthaud}, Alexie and {Li}, Cheng and {Lin}, Lihwai and {Maiolino}, Roberto and {Malanushenko}, Viktor and {Malanushenko}, Elena and {Mao}, Shude and {Maraston}, Claudia and {McDermid}, Richard M. and {Merrifield}, Michael R. and {Nichol}, Robert C. and {Oravetz}, Daniel and {Pan}, Kaike and {Parejko}, John K. and {Sanchez}, Sebastian F. and {Schlegel}, David and {Simmons}, Audrey and {Steele}, Oliver and {Steinmetz}, Matthias and {Thanjavur}, Karun and {Thompson}, Benjamin A. and {Tinker}, Jeremy L. and {van den Bosch}, Remco C.~E. and {Westfall}, Kyle B. and {Wilkinson}, David and {Wright}, Shelley and {Xiao}, Ting and {Zhang}, Kai},
        title = "{Overview of the SDSS-IV MaNGA Survey: Mapping nearby Galaxies at Apache Point Observatory}",
      journal = {\apj},
     keywords = {galaxies: evolution, galaxies: general, surveys, techniques: imaging spectroscopy, Astrophysics - Astrophysics of Galaxies},
         year = 2015,
        month = jan,
       volume = {798},
       number = {1},
          eid = {7},
        pages = {7},
          doi = {10.1088/0004-637X/798/1/7},
archivePrefix = {arXiv},
       eprint = {1412.1482},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015ApJ...798....7B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

# TODO: Add description of the dataset here
# You can copy an official description
_DESCRIPTION = """
An IFU dataset from the SDSS-IV MaNGA survey, a wide-field, optical, IFU survey of ~10,000
nearby galaxies. This dataset contains the following data products for each galaxy: the 3D data cubes,
and reconstructed griz images from the MaNGA Data Reduction Pipeline (DRP), and all the derived
analsysis maps from the MaNGA Data Analyis Pipeline (DAP).
"""

# TODO: Add a link to an official homepage for the dataset here
_HOMEPAGE = "https://www.sdss4.org/dr17/manga/"

# TODO: Add the licence for the dataset here if you can find it
_LICENSE = "BSD-3-Clause"

_VERSION = "1.0.0"


class MaNGA(datasets.GeneratorBasedBuilder):

    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(name="manga",
                               version=VERSION,
                               data_files=DataFilesPatternsDict.from_patterns({'train': ['out/manga/healpix=*/*.hdf5']}),
                               description="SDSS-IV MaNGA IFU datacubes and maps"),
    ]

    DEFAULT_CONFIG_NAME = "manga"

    _image_size = 96
    _image_filters = ['G', 'R', 'I', 'Z']
    _spectrum_size = 4563

    @classmethod
    def _info(cls):
        """ Defines features within the dataset """

        features = {}

        # add metadata features
        features['object_id'] = Value("string")
        features['ra'] = Value("float32")
        features['dec'] = Value("float32")
        features['healpix'] = Value("int16")
        features['z'] = Value("float32")
        features['spaxel_size'] = Value("float32")
        features['spaxel_size_units'] = Value("string")

        # add the spaxel features
        features['spaxels'] = [{
            "flux": Array2D(shape=(1, cls._spectrum_size), dtype='float32'),
            "ivar": Array2D(shape=(1, cls._spectrum_size), dtype='float32'),
            "mask": Array2D(shape=(1, cls._spectrum_size), dtype='int64'),
            "lsf": Array2D(shape=(1, cls._spectrum_size), dtype='float32'),
            "lambda": Array2D(shape=(1, cls._spectrum_size), dtype='float32'),
            "x": Value('int8'),
            "y": Value('int8'),
            'spaxel_idx': Value('int16'),
            "flux_units": Value('string'),
            "lambda_units": Value('string'),
            "skycoo_x": Value('float32'),
            "skycoo_y": Value('float32'),
            "ellcoo_r": Value('float32'),
            "ellcoo_rre": Value('float32'),
            "ellcoo_rkpc": Value('float32'),
            "ellcoo_theta": Value('float32'),
            "skycoo_units": Value('string'),
            "ellcoo_r_units": Value('string'),
            "ellcoo_rre_units": Value('string'),
            "ellcoo_rkpc_units": Value('string'),
            "ellcoo_theta_units": Value('string')
            }]

        # add the reconstructed image features
        features['images'] = [{
            'filter': Value('string'),
            'array': Array2D(shape=(cls._image_size, cls._image_size), dtype='float32'),
            'array_units': Value('string'),
            'psf': Array2D(shape=(cls._image_size, cls._image_size), dtype='float32'),
            'psf_units': Value('string'),
            'scale': Value('float32'),
            'scale_units': Value('string')
        }]

        # add the dap map features
        features['maps'] = [{
            "group": Value('string'),
            "label": Value('string'),
            "array": Array2D(shape=(cls._image_size, cls._image_size), dtype='float32'),
            "ivar": Array2D(shape=(cls._image_size, cls._image_size), dtype='float32'),
            "mask": Array2D(shape=(cls._image_size, cls._image_size), dtype='float32'),
            'array_units': Value('string')
        }]

        return datasets.DatasetInfo(
            # This is the description that will appear on the datasets page.
            description=_DESCRIPTION,
            # This defines the different columns of the dataset and their types
            features=Features(features),
            # Homepage of the dataset for documentation
            homepage=_HOMEPAGE,
            # License for the dataset if available
            license=_LICENSE,
            # Citation for the dataset
            citation=_CITATION,
        )

    def _split_generators(self, dl_manager):
        if not self.config.data_files:
            raise ValueError(f"At least one data file must be specified, but got data_files={self.config.data_files}")

        data_files = dl_manager.download_and_extract(self.config.data_files)

        if isinstance(data_files, (str, list, tuple)):
            files = data_files
            if isinstance(files, str):
                files = [files]

            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]

            return [
                datasets.SplitGenerator(
                    name=datasets.Split.TRAIN, gen_kwargs={"files": files}
                )
            ]

        splits = []
        for split_name, files in data_files.items():
            if isinstance(files, str):
                files = [files]

            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            splits.append(
                datasets.SplitGenerator(name=split_name, gen_kwargs={"files": files})
            )
        return splits

    def _generate_examples(self, files, object_ids=None):
        """ Yields examples as (key, example) tuples.
        """

        def reshape_if_needed(col_name, data):
            """ reshape the array data if it is in the list of columns """
            # define the columns that need reshaping - 1d array columns
            if col_name in {'flux', 'ivar', 'mask', 'lsf', 'lambda'}:
                return data.reshape(1, -1)
            return data

        for file in itertools.chain.from_iterable(files):
            with h5py.File(file, "r") as data:

                # loop over groups in hdf file
                for key in data.keys():
                    if object_ids and key not in object_ids:
                        continue

                    grp = data[key]
                    objid = grp['object_id'].asstr()[()]

                    example = {
                        'object_id': objid,
                        'ra': grp['ra'][()],
                        'dec': grp['dec'][()],
                        'healpix': grp['healpix'][()],
                        'z': grp['z'][()],
                        'spaxel_size': grp['spaxel_size'][()],
                        'spaxel_size_units': grp['spaxel_size_unit'].asstr()[()]

                    }

                    spax_cols = ('flux', 'ivar', 'mask', 'lsf', 'lambda', 'x', 'y', 'spaxel_idx', 'flux_units', 'lambda_units',
                                 'skycoo_x', 'skycoo_y', 'ellcoo_r', 'ellcoo_rre', 'ellcoo_rkpc', 'ellcoo_theta', 'skycoo_units',
                                 'ellcoo_r_units', 'ellcoo_rre_units', 'ellcoo_rkpc_units', 'ellcoo_theta_units')

                    example['spaxels'] = [
                        {col_name: reshape_if_needed(col_name, i[col_idx]) for col_idx, col_name in enumerate(spax_cols)}
                        for i in grp['spaxels']
                    ]

                    im_cols = ('filter', 'array', 'array_units', 'psf', 'psf_units', 'scale', 'scale_units')
                    example['images'] = [dict(zip(im_cols, i)) for i in grp['images']]

                    map_cols = ('group', 'label', 'array', 'ivar', 'mask', 'array_units')
                    example['maps'] = [dict(zip(map_cols, i)) for i in grp['maps']]

                    yield objid, example
