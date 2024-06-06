from astropy.table import Table
from snmachine import example_data, sndata, snfeatures, tsne_plot, snclassifier
from snmachine.utils.plasticc_pipeline import create_folder_structure, get_directories
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pickle

def translate_dataset_to_snmachine():
    from tqdm import tqdm
    from datasets import load_dataset
    from snmachine.sndata import EmptyDataset
    dset = load_dataset('../../scripts/yse/yse.py', trust_remote_code=True, split='train')
    dataset = EmptyDataset()
    for i in tqdm(range(len(dset['lightcurve']))):
        if dset['redshift'][i] == -99:
            continue
        N = len(dset['lightcurve'][i]['time'])
        t = Table(
                {
                    'mjd': dset['lightcurve'][i]['time'],
                    'filter': dset['lightcurve'][i]['band'],
                    'flux': dset['lightcurve'][i]['flux'],
                    'flux_error': dset['lightcurve'][i]['flux_err'],
                    'zp': [27.5 for _ in range(N)],
                    'zpsys': ['ab' for _ in range(N)],
                },
                meta = {
                    'name': dset['object_id'][i],
                    'type': dset['obj_type'][i],
                    'z': dset['redshift'][i]
                    }
                )
        dataset.insert_lightcurve(t)
    with open('snmachine_dset.pickle', 'wb') as f:
        pickle.dump(dataset, f)
    return dataset

def initialize():
    sns.set(font_scale=1., style="ticks")
    dataset_name = 'yse'
    analysis_name = os.path.join(f'output_{dataset_name}')
    folder_path = '.'
    if not os.path.exists(os.path.join(folder_path, analysis_name)):
        create_folder_structure(folder_path, analysis_name)
    directories = get_directories(folder_path, analysis_name)
    path_saved_features = directories['features_directory']
    path_saved_interm = directories['intermediate_files_directory']
    path_saved_classifier = directories['classifications_directory']
    run_name = os.path.join(path_saved_features, f'{dataset_name}_all')
    return dataset_name, analysis_name, directories, path_saved_features, path_saved_interm, path_saved_classifier, run_name

def load_and_clean_dataset(clip_mjd=True, pre_clip=30, post_clip=70, min_num_obs=5):
    # in lieu of loading example data
    if os.path.exists('snmachine_dset.pickle'):
        with open('snmachine_dset.pickle', 'rb') as f:
            dataset = pickle.load(f)
    else:
        dataset = translate_dataset_to_snmachine()
    weight_power = 4
    for name in dataset.data:
        # Scrubbing padding
        dataset.data[name].remove_rows(np.where(dataset.data[name]['mjd'] == 0))
        # clip around peak_mjd
        if clip_mjd:
            peak_mjd_guess = np.average(dataset.data[name]['mjd'], weights=dataset.data[name]['flux']**weight_power/dataset.data[name]['flux_error']**weight_power)
            dataset.data[name].remove_rows(np.where(dataset.data[name]['mjd'] < peak_mjd_guess - pre_clip))
            dataset.data[name].remove_rows(np.where(dataset.data[name]['mjd'] > peak_mjd_guess + post_clip))
        # normalizing to min(mjd) 0
        if len(dataset.data[name]['mjd']):
            dataset.data[name]['mjd'] -= min(dataset.data[name]['mjd'])
        # renaming filters
        dataset.data[name]['filter'] = dataset.data[name]['filter'].astype('<U6')
        for filt_before, filt_after in zip(('X', 'Y', 'g', 'r', 'i', 'z'), ('ztfg', 'ztfr', 'ps1::g', 'ps1::r', 'ps1::i', 'ps1::z')):
            idx = np.where(dataset.data[name]['filter'] == filt_before)
            dataset.data[name]['filter'][idx] = filt_after
    dataset.filter_set = ['ztfg', 'ztfr', 'ps1::g', 'ps1::r', 'ps1::i', 'ps1::z']
    # clip bad redshifts
    z = dataset.get_redshift()
    z_mask = z != -99
    for name in dataset.get_object_names()[~z_mask]:
        dataset.data.pop(name)
    dataset.object_names = dataset.object_names[z_mask]
    # clip num obs
    lengths_mask, bad_length_names = [], []
    for name in dataset.data:
        if len(dataset.data[name]['mjd']) >= min_num_obs:
            lengths_mask.append(True)
        else:
            lengths_mask.append(False)
            bad_length_names.append(name)
    for name in bad_length_names:
        dataset.data.pop(name)
    dataset.object_names = dataset.object_names[lengths_mask]
    # Fixing types to Ia, II, Ibc, other
    types = dataset.get_types()
    for i, typ_group in enumerate([
    ('SN-Ia-norm', 'SN-Ia-91bg-like', 'SN-Ia-91T-like', 'SN-Ia-SC', 'SN-Ia-CSM'),
    ('SN-IIb', 'SN-IIn', 'SN-II'),
    ('SN-Ib', 'SN-Ibn', 'SN-Ib-pec', 'SN-Ic', 'SN-Ic-BL'),
    ('SN-Iax', 'NA', 'LBV', 'SN', 'Other', 'LRN', 'TDE', 'SLSN-I', 'SLSN-II'),]):
        for typ in typ_group:
            types['Type'][types['Type']==typ] = i+1
    return dataset, types


def extract_salt_features(dataset, path_saved_interm, run_name):
    salt2_features_instance = snfeatures.TemplateFeatures(sampler='leastsq')
    if os.path.exists(f'{run_name}_templates.dat'):
        salt2_features = Table.read(f'{run_name}_templates.dat',
                                    format='ascii')
    else:
        salt2_features = salt2_features_instance.extract_features(
            dataset, use_redshift=True, number_processes=4,
            chain_directory=path_saved_interm)
        salt2_features.write(f'{run_name}_templates.dat', format='ascii')
    dataset.set_model(salt2_features_instance.fit_sn, salt2_features)
    for c in salt2_features.colnames[1:]:
        salt2_features[c][np.isnan(salt2_features[c])] = 0
    salt2_features_pd = salt2_features.to_pandas()
    salt2_features_pd.set_index('Object', inplace=True)
    return salt2_features, salt2_features_pd

def extract_wavelet_features(dataset, path_saved_interm, run_name):
    # Wavelet based feature extraction
    wavelet_features_instance = snfeatures.WaveletFeatures(output_root=path_saved_interm)
    if os.path.exists(f'{run_name}_wavelets.dat'):
        from snmachine import gps
        gps.read_gp_files_into_models(dataset, path_saved_interm)
        wavelet_features = Table.read(f'{run_name}_wavelets.dat', format='ascii')
        better_name = wavelet_features['Object'].astype(str)
        wavelet_features.replace_column('Object', better_name)
        pca_vals = np.loadtxt(f'{run_name}_wavelets_PCA_vals.dat')
        pca_vec = np.loadtxt(f'{run_name}_wavelets_PCA_vec.dat')
        pca_mean = np.loadtxt(f'{run_name}_wavelets_PCA_mean.dat')
    else:
        wavelet_features_pd = wavelet_features_instance.extract_features(
            dataset, number_gp=100, t_min=0, t_max=dataset.get_max_length(),
            number_processes=1, gp_dim=1, number_comps=20,
            output_root=path_saved_interm)

        # Change to astropy table
        wavelet_features = Table(wavelet_features_pd.to_numpy())
        column_names = wavelet_features.colnames
        column_names.insert(0, 'Object')
        objs_name = list(wavelet_features_pd.index)
        wavelet_features.add_column(objs_name, name='Object')
        wavelet_features = wavelet_features[column_names]
        wavelet_features.write(f'{run_name}_wavelets.dat', format='ascii')

        # Load from where they were saved
        pca_vals = np.load(os.path.join(path_saved_interm,
                           'eigenvalues.npy'))
        pca_vec = np.load(os.path.join(path_saved_interm,
                          'eigenvectors.npy'))
        pca_mean = np.load(os.path.join(path_saved_interm, 'means.npy'))

        # Save in a different way
        np.savetxt(f'{run_name}_wavelets_PCA_vals.dat', pca_vals)
        np.savetxt(f'{run_name}_wavelets_PCA_vec.dat', pca_vec)
        np.savetxt(f'{run_name}_wavelets_PCA_mean.dat', pca_mean)

    wavelet_features_pd = wavelet_features.to_pandas()
    wavelet_features_pd.set_index('Object', inplace=True)
    dataset.set_model(wavelet_features_instance.fit_sn,
                      wavelet_features_pd, dataset, 'sym2',
                      path_saved_interm, dataset.filter_set)
    return wavelet_features, wavelet_features_pd

def classify_features(features_pd, types, path_saved_classifier, which_column=0, train_set=0.3, number_processes=4, classifier_list=['svm', 'random_forest', 'boost_dt', 'boost_rf'], figname='default'):
    data_labels = types.to_pandas()
    data_labels.set_index('Object', inplace=True)
    data_labels = data_labels['Type']
    print('The available classifiers are:', snclassifier.choice_of_classifiers)
    fig = plt.figure()
    classifier_instances, cms = snclassifier.run_several_classifiers(
        classifier_list=classifier_list,
        features=features_pd, labels=data_labels,  scoring='accuracy',
        train_set=train_set, scale_features=True, which_column=which_column,
        random_seed=42, output_root=path_saved_classifier,
        **{'plot_roc_curve': True, 'number_processes': number_processes})
    if figname == 'default':
        figname = 'classification.png')
    fig.savefig(figname)
    plt.show()

if __name__=='__main__':
    salt_classification = False
    wavelet_classification = True
    dataset, types = load_and_clean_dataset()
    dataset_name, analysis_name, directories, path_saved_features, path_saved_interm, path_saved_classifier, run_name = initialize()
    if salt_classification:
        salt_features, salt2_features_pd = extract_salt_features(dataset, path_saved_interm, run_name)
        fig = plt.figure(figsize=(5,4))
        tsne_plot.plot(salt_features, types['Type'], type_dict={'1':'SN Ia', '2':'SN II', '3': 'SN Ibc', '4': 'Other'})
        fig.savefig('salt_tsne.png')
        classify_features(salt_features_pd, types, path_saved_classifier, figname='salt_classification')
    if wavelet_classification:
        wavelet_features, wavelet_features_pd = extract_wavelet_features(dataset, path_saved_interm, run_name)
        fig = plt.figure(figsize=(5,4))
        tsne_plot.plot(wavelet_features, types['Type'], type_dict={'1':'SN Ia', '2':'SN II', '3': 'SN Ibc', '4': 'Other'})
        fig.savefig('wavelet_tsne.png')
        classify_features(wavelet_features_pd, types, path_saved_classifier, figname='wavelet_classification.png')
